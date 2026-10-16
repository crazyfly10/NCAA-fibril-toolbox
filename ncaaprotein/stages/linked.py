"""Run 1: the protein with the residue attached. Build (steps 13-17), run (18-20), analyse (21)."""
import json
import re
import shutil
from pathlib import Path

import numpy as np

from .. import amber, gmxmd
from ..box import advice, atoms_in, auto_buffer, built
from ..core import StageError, nsteps_for, remove_spec, restraint_atoms, restraint_tag, temp_changes, pdb_atoms, pdb_chains, resindex_ranges, stage, write_mdp


def ncaa_info(ctx):
    """Where the residue's files are, and what ncaa-param recorded about it."""
    n = ctx.cfg["ncaa"]
    files = ctx.input(n["files_from"]) if "files_from" in n else ctx.other("ncaa")
    desc = json.loads((files / "ncaa_map.json").read_text())
    relax = (ctx.input(n["relax_from"]) if "relax_from" in n else ctx.other("ncaa_relax")) / f"ace_{n['name'].lower()}_nme.amb2gmx"
    return n["name"], files, relax, desc


def ncaa_reach(c, relax):
    """The residue's reach in nm: `reach` is "median", "p95" (default) or "max" from its run in water, "extended" for
    its fully stretched length (no simulation behind it, so the same every time), or a number."""
    q = c.get("reach", "p95")
    if isinstance(q, (int, float)):
        return float(q)
    keys = {"max": "max_nm", "p95": "p95_nm", "median": "median_nm", "extended": "extended_nm"}
    if q not in keys:
        raise StageError(f'reach has to be "median", "p95", "max", "extended" or a number in nm, not {q!r}')
    data = json.loads((relax / "reach.json").read_text())
    if keys[q] not in data:
        raise StageError(f'reach.json has no {keys[q]}: run ncaa-relax again to have it worked out')
    return data[keys[q]]


def reach_caveat(ctx, c, relax):
    """A line for the report when the automatic box rests on the median or 95% reach of a short residue run."""
    f = relax / "reach.json"
    if isinstance(c.get("reach", "p95"), (int, float)) or c.get("reach", "p95") in ("max", "extended") or not f.exists():
        return
    ns = json.loads(f.read_text()).get("run_ns")
    if ns is not None and ns < 5:
        ctx.note(f"auto box: the reach comes from a {ns:g} ns run of the residue, too short for this value to be settled: the box "
                 "will differ a little from one such run to the next. Run ncaa-relax for 10 ns or more, or set reach = \"extended\"")


def protein_pdb(ctx, section):
    c = ctx.cfg.get(section, {})
    p = ctx.input(c["protein"]) if "protein" in c else ctx.other("protein_relax") / "protein.pdb"
    if not p.exists():
        raise StageError(f"{p} not found: run protein-prepare and protein-run, or give [{section}] protein = \"your_relaxed_protein.pdb\"")
    return p


def restraint(ctx, section, chains):
    """The production restraint as an index group -> (tag, group name, selection, fc). Either `chains` (each a whole
    chain "G" or a stretch of one "G:226-282", residue numbers as in the structure) with `atoms`, or `select`: any
    gmx select expression, on the built system (residues numbered 1..N through the solute)."""
    r = ctx.cfg.get(section, {}).get("restrain", {})
    fc = float(r.get("fc", 200))
    if r.get("select"):
        tag = re.sub(r"[^A-Za-z0-9]+", "_", str(r.get("name", "SEL")))
        return tag, f"{tag}_posre", str(r["select"]), fc
    if not r.get("chains"):
        return None
    tag = restraint_tag(r["chains"])
    atoms = restraint_atoms(r.get("atoms", "N CA C O"))
    rng = [resindex_ranges(chains, ch) for ch in r["chains"]]
    sel = "(" + " or ".join(f"resindex {a} to {b}" for a, b in rng) + ")"
    if atoms == "heavy":
        sel += ' and not name "H*"'
    elif atoms != "all":
        sel += " and name " + " ".join(atoms)
    return tag, f"{tag}_posre", sel, fc


def prep_for_tleap(ctx, d, src, dst):
    rm = remove_spec(ctx.cfg)          # [remove], whichever `when`: a relaxed structure given by hand may still have it
    out = ctx.script("prep_protein_tleap.py", [src, dst], dict(DROP_RESI=rm[1] if rm else -1, DROP_RESN=rm[0] if rm else "none"))
    ctx.grab(out, re.compile(r"^wrote|may not know"))


def build_gromacs(ctx, d, out, waters):
    """acpype and its repairs, shared by the linked and the free-NCAA builds."""
    name = ctx.cfg["ncaa"]["name"]
    g = amber.acpype(ctx, d, out)
    top = (g / f"{out}_GMX.top").read_text()
    mol, natoms, nres = amber.first_molecule(g / f"{out}_GMX.top")
    if f"N-{name}" not in top:
        raise StageError(f"{out}_GMX.top has no CMAP type for {name}: is acpype 2026.9.3 or newer, and was frcmod.ff19SB_{name} loaded?")
    ncmap = len(re.findall(r"^\s*\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+1\s*$", top.split("[ cmap ]")[-1].split("[ ")[0], re.M)) if "[ cmap ]" in top else 0
    ctx.note(f"{out}_GMX.top: first molecule '{mol}'; CMAP type N-{name} defined; {ncmap} CMAP terms in the solute")
    if waters + nres > 99999:
        ctx.grab(ctx.script("fix_gro_columns.py", [f"{out}_GMX.gro"], cwd=g), re.compile(r"re-aligned|aligned"))
    ctx.note(f"{out}_GMX.gro: {int(open(g / f'{out}_GMX.gro').readlines()[1])} atoms, box {' x '.join(gmxmd.box_of(g / f'{out}_GMX.gro')[:3])} nm")
    return g


def equilibrate(ctx, g, out, group, rs, npt_mdp, density, ref_t, sec=None, where="[linked]", gate_em=True, num=18):
    """EM, NVT, NPT with every solute heavy atom restrained.
    [linked] / [free] nvt_ps and npt_ps change the lengths; the defaults are 100 ps and 1 ns."""
    tc = {"tc-grps": f"{group} Water_and_ions"}
    sec = sec or {}
    T = sec.get("temperature")                                 # None = the mdp files' 298 K
    nvt_n = {"nsteps": nsteps_for("linked_nvt.mdp", float(sec["nvt_ps"]) / 1000)} if "nvt_ps" in sec else {}
    npt_n = {"nsteps": nsteps_for(npt_mdp, float(sec["npt_ps"]) / 1000)} if "npt_ps" in sec else {}
    if ctx.begin("em", f"{num} · Energy minimisation"):
        write_mdp("linked_em.mdp", g / "linked_em.mdp")
        gmxmd.run_em(ctx, g, "linked_em.mdp", f"{out}_GMX.gro", "topol.top", n="index.ndx", maxwarn=1, second_pass=False)
        ctx.end()
    gsec = (where, sec)
    if gate_em:
        ctx.check_gate("em", gsec, "the minimisation", "em.png", "NVT")
    if ctx.begin("nvt", f"{num} · NVT with the solute restrained"):
        write_mdp("linked_nvt.mdp", g / "linked_nvt.mdp", {**tc, **nvt_n, **temp_changes("linked_nvt.mdp", T)})
        gmxmd.run_md(ctx, g, "nvt", "linked_nvt.mdp", "em.gro", "topol.top", r="em.gro", n="index.ndx", maxwarn=1)
        gmxmd.judge(ctx, g, "nvt", [group, "Water_and_ions"], ref_t, table=sec)
        ctx.end()
    ctx.check_gate("nvt", gsec, "the NVT equilibration", "nvt.png", "NPT")
    if ctx.begin("npt", f"{num} · NPT with the solute restrained"):
        write_mdp(npt_mdp, g / npt_mdp, {**tc, **npt_n, **temp_changes(npt_mdp, T)})
        gmxmd.run_md(ctx, g, "npt", npt_mdp, "nvt.gro", "topol.top", r="nvt.gro", t="nvt.cpt", n="index.ndx", maxwarn=1)
        gmxmd.judge(ctx, g, "npt", [group, "Water_and_ions"], ref_t, gmxmd.density_range(ctx, density, g / "nvt.gro"), table=sec)
        ctx.note(f"box before NPT {' x '.join(gmxmd.box_of(g / 'nvt.gro')[:3])}, after {' x '.join(gmxmd.box_of(g / 'npt.gro')[:3])} nm")
        ctx.end()
    ctx.check_gate("npt", gsec, "the NPT equilibration", "npt.png", "production")
    ctx.pause(gsec, "production")


def production(ctx, g, group, rs, mdp, ns, extra=None):
    """Production with only the chosen chains' backbone restrained."""
    changes = {"tc-grps": f"{group} Water_and_ions", "nsteps": nsteps_for(mdp, ns), "define": None, **(extra or {})}
    if rs:
        tag, rgroup, _, fc = rs
        ctx.sh(f"{ctx.gmx} genrestr -f npt.gro -n index.ndx -o posre_{tag}.itp -fc {fc:g} {fc:g} {fc:g}", cwd=g, stdin=rgroup + "\n", log="genrestr")
        k = sum(1 for l in open(g / f"posre_{tag}.itp") if re.match(r"\s*\d", l))
        ctx.note(f"posre_{tag}.itp: {k} atoms ({rgroup}) at {fc:g} kJ/mol/nm2")
        changes["define"] = f"-DPOSRES_{tag}"
    write_mdp(mdp, g / mdp, changes)
    out = gmxmd.grompp(ctx, g, mdp, "npt.gro", "topol.top", "prod.tpr", r="npt.gro", t="npt.cpt", n="index.ndx", maxwarn=1)
    ctx.grab(out, re.compile(r"minimum mass"))
    gmxmd.mdrun(ctx, g, "prod")
    perf = gmxmd.performance(g, "prod")
    ctx.note(f"production: {ns:g} ns finished, {gmxmd.lincs(g, 'prod')} LINCS warnings" + (f", {perf:.0f} ns/day" if perf else ""))
    en = gmxmd.energy(ctx, g, "prod.edr", ["Temperature", f"T-{group}", "Pressure", "Density"], "prod_energy.xvg")
    ctx.grab(en, re.compile(r"^(Temperature|T-|Pressure|Density)"))


@stage("linked-build", "linked", "13-17")
def linked_build(ctx):
    """relaxed protein + the residue at chosen sites -> solvated GROMACS system"""
    L, d = ctx.cfg.get("linked", {}), ctx.dir
    name, files, relax, desc = ncaa_info(ctx)
    out = L.get("name", "linked")
    group = L.get("group", "Solute")
    sites = [[s["chain"], int(s["resi"])] for s in ctx.cfg["attach"]["sites"]]
    graft = f"{name.lower()}.pdb"

    if ctx.begin("s13", "13 · Place the residue on the protein (PyMOL)"):
        shutil.copy(protein_pdb(ctx, "linked"), d / "protein.pdb")
        shutil.copy(relax / graft, d / graft)
        o = ctx.script("place_ncaa.py", pymol=True, overrides=dict(
            PROTEIN="protein.pdb", GRAFT=graft, OUT="protein_ncaa.pdb", RESN=name, SITES=sites, CHI=desc["chi"],
            CORE=desc["backbone_names"] + desc["parent_side_names"], MAX_FIT=float(ctx.cfg["attach"].get("max_fit", 0.6))))
        ctx.grab(o, re.compile(r"fit RMSD|best|wrote"))
        # clash test: the backbone fit can be perfect while the rest of the residue runs into the protein
        limit = float(ctx.cfg["attach"].get("clash_within", 1.5))
        if limit not in (1.5, 2.0, 2.5):
            raise StageError("[attach] clash_within has to be 1.5, 2.0 or 2.5 (A)")
        buried = []
        for (ch, resi), m in zip(sites, re.finditer(r"best \((\d+), (\d+), (\d+)\)", o)):
            n = dict(zip((1.5, 2.0, 2.5), (int(x) for x in m.groups())))
            bad = n[limit]
            verdict = "FAIL, the site is buried" if bad else ("PASS, but the site is tight" if n[2.5] else "PASS")
            ctx.note(f"clash test, site {ch}{resi}: at the best rotamer {n[1.5]} / {n[2.0]} / {n[2.5]} heavy atoms of the {name} lie "
                     f"within 1.5 / 2.0 / 2.5 A of the protein; limit: none within {limit:g} A  ->  {verdict}")
            if bad:
                buried.append(f"{ch}{resi}")
        if buried and not ctx.cfg["attach"].get("allow_clash", False):
            raise StageError(f"the {name} clashes with the protein at {', '.join(buried)}: its backbone fits, but there is no room for "
                             "the rest of it, so the residue there is buried and cannot be used as an attachment point. Choose a "
                             "residue whose side chain is exposed (check its SASA in the relaxed protein) or redesign the protein. "
                             "To build it anyway, set [attach] allow_clash = true")
        ca = [(l[21], int(l[22:26]), l[17:21].strip()) for l in pdb_atoms(d / "protein_ncaa.pdb") if l[12:16].strip() == "CA"]
        for ch, resi in sites:
            k = [i for i, x in enumerate(ca) if x[:2] == (ch, resi)]
            ok = len(k) == 1 and ca[k[0]][2] == name and (k[0] == 0 or ca[k[0] - 1][0] != ch or ca[k[0] - 1][1] < resi)
            ctx.note(f"{name} at {ch}{resi}: " + ("in its own residue position" if ok else "NOT in sequence -- tleap would bond it wrongly"))
            if not ok:
                raise StageError(f"{name} is not in place at {ch}{resi}")
        ctx.end()

    if ctx.begin("s14", "14 · AMBER names and TER records"):
        prep_for_tleap(ctx, d, "protein_ncaa.pdb", "protein_ncaa_noh.pdb")
        k = sum(1 for l in pdb_atoms(d / "protein_ncaa_noh.pdb") if l[17:20] == name)
        ctx.note(f"{name} heavy atoms: {k} (expected {desc['heavy_atoms'][name]} per residue, {desc['heavy_atoms'][name] * len(sites)} for {len(sites)})")
        ctx.end()

    if ctx.begin("s15", "15 · Build, solvate and add salt in tleap"):
        amber.copy_ncaa_files(ctx, files, d, name)
        box = L.get("box", "auto")
        try:
            r_nm = ncaa_reach(L, relax)
        except FileNotFoundError:
            r_nm = None                      # no ncaa-relax run to read the reach from: a number has to be given
        if r_nm is None and box == "auto":
            raise StageError("box = \"auto\" needs the reach measured by ncaa-relax (reach.json): run it, or give [linked] reach or box as a number")
        if r_nm is not None:
            parts = [(ch, resi, r_nm) for ch, resi in sites]
            parts += [(f["chain"], int(f["anchor"]), float(f["reach_nm"])) for f in L.get("flexible", [])]
            abox, lines, ext, stretched = auto_buffer(d / "protein_ncaa_noh.pdb", parts, float(L.get("image_margin_nm", 2.0)))
            fixed = 12.0 if box == "auto" else float(box)
            if box == "auto":
                for l in lines:
                    ctx.note("auto box: " + l)
                reach_caveat(ctx, L, relax)
            for l in advice("[linked] box", box, atoms_in((ext + 2 * np.array(abox)) / 10), atoms_in((ext + 2 * fixed) / 10), f"a {fixed:g} A buffer",
                            f"The attached {name} reaches {r_nm:.2f} nm from its backbone.", "the periodic image check of step 20"):
                ctx.note(l)
            if box == "auto":
                box = abox
        charge, waters, pairs = amber.tleap_build(ctx, d, name, out, "prot = loadPdb protein_ncaa_noh.pdb\n", "prot", box, "Box",
                                                  float(L.get("salt", 0.15)), extra_save="savePdb prot build_check.pdb\n")
        k = sum(1 for l in open(d / "build_check.pdb") if l[17:20] == name and l.startswith(("ATOM", "HETATM")))
        ctx.note(f"build_check.pdb: {k} {name} atoms (expected {desc['n_residue_atoms'] * len(sites)} = {len(sites)} x {desc['n_residue_atoms']})")
        shape = dict(extent_nm=[float(x) / 10 for x in ext], stretched_nm=[float(x) / 10 for x in stretched]) if r_nm is not None else {}
        (d / "build_info.json").write_text(json.dumps(dict(charge=charge, waters=waters, pairs=pairs, box=box, **shape)) + "\n")
        ctx.end()

    info = json.loads((d / "build_info.json").read_text())
    if ctx.begin("s16", "16 · Convert to GROMACS with acpype"):
        g = build_gromacs(ctx, d, out, info["waters"])
        if "extent_nm" in info:
            for l in built([float(x) for x in gmxmd.box_of(g / f"{out}_GMX.gro")[:3]], info["extent_nm"],
                           f"the protein with its {name}", info["stretched_nm"]):
                ctx.note(l)
        ctx.end()

    g = d / f"{out}.amb2gmx"
    if ctx.begin("s17", "17 · Thermostat group, restraint files and topology includes"):
        chains = pdb_chains(d / "protein_ncaa_noh.pdb")
        rs = restraint(ctx, "linked", chains)
        _, natoms, nres = amber.restraint_files(ctx, g, out, f"POSRES_{rs[0]}" if rs else "POSRES_NONE")
        groups = {"System": "all", group: f"resindex 1 to {nres}", "Water_and_ions": f"not resindex 1 to {nres}"}
        if rs:
            groups[rs[1]] = rs[2]
        amber.make_index(ctx, g, f"{out}_GMX.gro", groups)
        ctx.end()


@stage("linked-run", "linked", "18-20")
def linked_run(ctx):
    """EM, NVT, NPT, production with the chosen chains restrained, then the periodic-image check"""
    L = ctx.cfg.get("linked", {})
    out, group = L.get("name", "linked"), L.get("group", "Solute")
    g = ctx.dir / f"{out}.amb2gmx"
    rs = restraint(ctx, "linked", pdb_chains(ctx.dir / "protein_ncaa_noh.pdb"))
    equilibrate(ctx, g, out, group, rs, "linked_npt.mdp", L.get("density_range"), float(L.get("temperature", 298.0)), L)
    if ctx.begin("s19", "19 · Production with the chosen chains restrained"):
        production(ctx, g, group, rs, "linked_prod.mdp", float(L.get("production_ns", 200)),
                   temp_changes("linked_prod.mdp", L.get("temperature")))
        ctx.end()
    if ctx.begin("s20", "20 · Periodic image distance of the protein with its NCAAs"):
        dmin = gmxmd.image_check_single(ctx, g, "prod.tpr", "prod.xtc", "index.ndx", group, int(L.get("image_check_dt_ps", 100)))
        if dmin < 2.0:
            b = json.loads((ctx.dir / "build_info.json").read_text())["box"]
            add = (2.0 - dmin) * 10 / 2 / 0.95
            ctx.note(f"to clear 2.0 nm, rebuild with at least {add:.0f} A more buffer per side (the box was built with {b})" +
                     ("; the two parts touched, so they could have reached further than this: box = \"auto\" sizes the "
                      "buffer from the residue's reach instead" if dmin < 0.5 else ""))
        ctx.end()


@stage("linked-analyse", "linked", "21")
def linked_analyse(ctx):
    """conformations and adsorption of each attached residue; from when the run is usable (linked_eq_check.py)"""
    L = ctx.cfg.get("linked", {})
    A = L.get("analysis", {})
    name, _, _, desc = ncaa_info(ctx)
    out, group = L.get("name", "linked"), L.get("group", "Solute")
    g = ctx.dir / f"{out}.amb2gmx"
    if ctx.begin("s21", "21 · Conformations and adsorption of the attached residues"):
        ctx.sh(f"{ctx.gmx} trjconv -s prod.tpr -f prod.cpt -n index.ndx -o solute.gro", cwd=g, stdin=group + "\n", log="solute_gro")
        ctx.sh(f"{ctx.gmx} trjconv -s prod.tpr -f prod.xtc -n index.ndx -o solute.xtc", cwd=g, stdin=group + "\n", log="solute_xtc")
        chains = pdb_chains(ctx.dir / "protein_ncaa_noh.pdb")
        names, short = {}, {}
        for s in ctx.cfg["attach"]["sites"]:
            k = resindex_ranges(chains, f"{s['chain']}:{s['resi']}-{s['resi']}")[0]
            names[f"{name}{k}"] = s.get("label", f"chain {s['chain']}") + f", residue {s['resi']}"
            short[f"{name}{k}"] = s.get("label", f"chain {s['chain']}")
        # the residue as the analysis needs it: what it shares with its parent, where the attached part starts,
        # and where it ends (the user's choice; otherwise the atom most bonds from the anchor)
        core = desc["backbone_names"] + desc["parent_side_names"]
        ov = dict(RESNAME=name, CORE=core, ANCHOR=A.get("anchor", (desc["parent_side_names"] or ["CA"])[-1]), TIP=A.get("tip"),
                  NEAR_BONDS=int(A.get("near_bonds", 7)), TIP_BONDS=int(A.get("tip_bonds", 7)), SITE_NAMES=names, SITE_SHORT=short)
        for key, const in (("contact_A", "CONTACT"), ("on_min", "ON_MIN"), ("start_step_ns", "START_STEP_NS"), ("block_ns", "BLOCK_NS"),
                           ("tail_ns", "TAIL_NS"), ("smooth_ns", "SMOOTH_NS")):
            if key in A:
                ov[const] = A[key]
        args = ["solute.gro", "solute.xtc"] + ([str(A["every_ps"])] if "every_ps" in A else [])
        o = ctx.script(A.get("script", "linked_eq_check.py"), args, ov, cwd=g)
        (g / "linked_eq_report.txt").write_text(o)
        gmxmd.draw(ctx, g, A.get("plot_script", "plot_linked_eq.py"))
        ctx.note(f"sites: {names}")
        ctx.grab(o, re.compile(r"^each |lying on it|latest of these|for the bonded distributions|tip from|too short to judge|halves still differ"))
        ctx.note(f"full report: {g / 'linked_eq_report.txt'}; figure and table: linked_eq.png, linked_eq.csv")
        ctx.end()
