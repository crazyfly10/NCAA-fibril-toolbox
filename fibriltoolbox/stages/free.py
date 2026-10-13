"""Run 2: the fibril with free copies of the capped residue around it. Build, run (33-34),
analyse (35)."""
import json
import re
import shutil

import numpy as np

from .. import amber, gmxmd
from ..box import advice, atoms_in, auto_gap_buffer, built
from ..core import StageError, pdb_atoms, pdb_chains, stage, temp_changes
from .linked import build_gromacs, equilibrate, fibril_pdb, ncaa_info, ncaa_reach, reach_caveat, prep_for_tleap, production, restraint


@stage("free-build", "free", "28-32")
def free_build(ctx):
    """relaxed fibril + N free capped residues -> solvated GROMACS system"""
    F, d = ctx.cfg.get("free", {}), ctx.dir
    name, files, relax, desc = ncaa_info(ctx)
    out, group = F.get("name", "free"), F.get("group", "Solute")
    n_link = int(F.get("n_linkers", 40))
    frames = ctx.cfg["ncaa"].get("relax", {}).get("frames_ps", [250, 500, 750, 1000])

    if ctx.begin("s28", "28 · The fibril in a box, ready for tleap"):
        shutil.copy(fibril_pdb(ctx, "free"), d / "fibril.pdb")
        ctx.sh(f"{ctx.gmx} editconf -f fibril.pdb -o fibril_box.pdb -princ -c -d {float(F.get('insert_d', 1.5))} -bt triclinic",
               stdin="Protein\nSystem\n", log="editconf")
        prep_for_tleap(ctx, d, "fibril_box.pdb", "fib_noh.pdb")
        ctx.end()

    if ctx.begin("s29", f"29 · Insert {n_link} free capped residues around the fibril"):
        if n_link % len(frames):
            raise StageError(f"n_linkers = {n_link} is not a multiple of the {len(frames)} starting frames")
        per = n_link // len(frames)
        shutil.copy(d / "fibril_box.pdb", d / "start_0.pdb")
        for k, t in enumerate(frames, 1):
            shutil.copy(relax / f"lnk_{t}.pdb", d / f"lnk_{t}.pdb")
            o = ctx.sh(f"{ctx.gmx} insert-molecules -f start_{k - 1}.pdb -ci lnk_{t}.pdb -nmol {per} -try {int(F.get('insert_tries', 1000))} -o start_{k}.pdb", log=f"insert_{t}")
            m = re.search(r"Added (\d+) molecules \(out of (\d+) requested\)", o)
            ctx.note(f"frame {t} ps: " + (m.group(0) if m else "no 'Added' line in the output"))
            if not m or m.group(1) != m.group(2):
                raise StageError(f"insert-molecules placed fewer than {per} molecules: enlarge [free] insert_d, or raise insert_tries")
        shutil.copy(d / f"start_{len(frames)}.pdb", d / "start.pdb")
        seq = ("ACE", name, "NME")
        o = ctx.script("split_linkers.py", ["fibril_box.pdb", "start.pdb", "linkers_noh.pdb"],
                       dict(SEQ=f"py:{seq!r}", HEAVY=desc["heavy_atoms"]))
        ctx.grab(o, re.compile(r"^linkers|^wrote"))
        ctx.end()

    if ctx.begin("s30", "30 · Build, solvate and add salt in tleap"):
        amber.copy_ncaa_files(ctx, files, d, name)
        body = f"fib = loadPdb fib_noh.pdb\nlnk = loadPdb linkers_noh.pdb\n{out} = combine {{ fib lnk }}\n"
        box = F.get("box", "auto")
        try:
            r_nm = ncaa_reach(F, relax)
        except FileNotFoundError:
            r_nm = None
        if r_nm is None and box == "auto":
            raise StageError("box = \"auto\" needs the reach measured by ncaa-relax (reach.json): run it, or give [free] reach or box as a number")
        if r_nm is not None:
            margin = float(F.get("image_margin_nm", 2.0))
            abox, lines, ext, cext = auto_gap_buffer(d / "fib_noh.pdb", [d / "linkers_noh.pdb"], r_nm + margin)
            fixed = 12.0 if box == "auto" else float(box)
            if box == "auto":
                ctx.note(f"auto box: a free {name} reaches {r_nm:.2f} nm; margin {margin:.1f} nm")
                for l in lines:
                    ctx.note("auto box: " + l)
                reach_caveat(ctx, F, relax)
            for l in advice("[free] box", box, atoms_in((ext + 2 * np.array(abox)) / 10), atoms_in((ext + 2 * fixed) / 10), f"a {fixed:g} A buffer",
                            f"A free {name} reaches {r_nm:.2f} nm, and one lying on the fibril could then touch the fibril's image too.",
                            "gmx mindist -pi on the Fibril group after the run"):
                ctx.note(l)
            if box == "auto":
                box = abox
        charge, waters, pairs = amber.tleap_build(ctx, d, name, out, body, out, box, "Box", float(F.get("salt", 0.15)))
        shape = dict(fibril_nm=[float(x) / 10 for x in cext]) if r_nm is not None else {}
        (d / "build_info.json").write_text(json.dumps(dict(charge=charge, waters=waters, pairs=pairs, box=box, **shape)) + "\n")
        ctx.end()

    info = json.loads((d / "build_info.json").read_text())
    if ctx.begin("s31", "31 · Convert with acpype, then its repairs"):
        g0 = build_gromacs(ctx, d, out, info["waters"])
        if "fibril_nm" in info:
            for l in built([float(x) for x in gmxmd.box_of(g0 / f"{out}_GMX.gro")[:3]], info["fibril_nm"], "the fibril"):
                ctx.note(l)
            ctx.note(f"the free {name} molecules fill the water around the fibril and cross the box faces; the distances above are the fibril's own")
        ctx.end()

    g = d / f"{out}.amb2gmx"
    if ctx.begin("s32", "32 · Restraint file, include lines and index groups"):
        chains = pdb_chains(d / "fib_noh.pdb")
        nfib = sum(len(r) for _, r in chains)
        rs = restraint(ctx, "free", chains)
        _, natoms, nres = amber.restraint_files(ctx, g, out, f"POSRES_{rs[0]}" if rs else "POSRES_NONE")
        if nres != nfib + 3 * n_link:
            raise StageError(f"the solute has {nres} residues, expected {nfib} + 3 x {n_link}")
        groups = {"System": "all", "Fibril": f"resindex 1 to {nfib}", "Linkers": f"resindex {nfib + 1} to {nres}",
                  group: f"resindex 1 to {nres}", "Water_and_ions": f"not resindex 1 to {nres}"}
        if rs:
            groups[rs[1]] = rs[2]
        amber.make_index(ctx, g, f"{out}_GMX.gro", groups)
        ctx.end()


@stage("free-run", "free", "33-34")
def free_run(ctx):
    """EM, NVT, NPT (4 fs, repartitioned), production with the chosen chains restrained"""
    F = ctx.cfg.get("free", {})
    out, group = F.get("name", "free"), F.get("group", "Solute")
    g = ctx.dir / f"{out}.amb2gmx"
    rs = restraint(ctx, "free", pdb_chains(ctx.dir / "fib_noh.pdb"))
    ctx.say("   (step 33: on a system this size minimisation is expected to end NOT converged; it does not hold the stage, NVT is the judge)")
    equilibrate(ctx, g, out, group, rs, "free_npt.mdp", F.get("density_range"), float(F.get("temperature", 298.0)), F,
                where="[free]", gate_em=False)
    if ctx.begin("s34", "34 · Production with the chosen chains restrained"):
        production(ctx, g, group, rs, "free_prod.mdp", float(F.get("production_ns", 200)), {"compressed-x-grps": group, **temp_changes("free_prod.mdp", F.get("temperature"))})
        ctx.end()


@stage("free-analyse", "free", "35")
def free_analyse(ctx):
    """has the free-molecule run settled, do the molecules exchange, which part makes the contacts (free_linker_analysis.py)"""
    F = ctx.cfg.get("free", {})
    A = F.get("analysis", {})
    LA = ctx.cfg.get("linked", {}).get("analysis", {})          # the residue is the same one: share its description
    name, _, _, desc = ncaa_info(ctx)
    out, group = F.get("name", "free"), F.get("group", "Solute")
    g = ctx.dir / f"{out}.amb2gmx"
    if ctx.begin("s35", "35 · Analysis: has it settled, and do the free molecules exchange?"):
        ctx.sh(f"{ctx.gmx} trjconv -s prod.tpr -f prod.cpt -n index.ndx -o solute.gro", cwd=g, stdin=group + "\n", log="solute_gro")
        args = ["solute.gro", "prod.xtc", "--every", str(A.get("every_ps", 10))]
        if "from_ns" in A:
            args += ["--from", str(A["from_ns"])]
        pick = lambda key, default=None: A.get(key, LA.get(key, default))
        ov = dict(RESNAME=name, CORE=desc["backbone_names"] + desc["parent_side_names"],
                  ANCHOR=pick("anchor", (desc["parent_side_names"] or ["CA"])[-1]), TIP=pick("tip"),
                  NEAR_BONDS=int(pick("near_bonds", 7)), TIP_BONDS=int(pick("tip_bonds", 7)))
        if "part_names" in A:
            if len(A["part_names"]) != 4:
                raise StageError("[free.analysis] part_names needs four names: backbone and caps, near the anchor, tip, the rest")
            ov["PARTS"] = list(A["part_names"])
        for key, const in (("contact_A", "CONTACT"), ("min_visit_ps", "MIN_VISIT_PS"), ("block_ps", "BLOCK_PS"), ("smooth_ns", "SMOOTH_NS")):
            if key in A:
                ov[const] = A[key]
        o = ctx.script(A.get("script", "free_linker_analysis.py"), args, ov, cwd=g)
        ctx.grab(o, re.compile(r"\[(PASS|FAIL|WARN|info)\]|^linker parts"))
        gmxmd.draw(ctx, g, A.get("plot_script", "plot_free_linker.py"), ["--separate"] if A.get("separate_figures") else [])
        ctx.note(f"full report: {g / 'free_linker_report.txt'}; figure: free_linker_overview.png; four CSV tables beside it")
        ctx.end()
