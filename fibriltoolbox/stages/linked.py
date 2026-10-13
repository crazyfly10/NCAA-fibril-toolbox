"""Run 1: the fibril with the residue attached. Build (protocol steps 18-21, 23), run (22, 24-26), analyse (27)."""
import json
import re
import shutil
from pathlib import Path

import numpy as np

from .. import amber, gmxmd
from ..box import advice, atoms_in, auto_buffer, built
from ..core import StageError, nsteps_for, pdb_atoms, pdb_chains, resindex_ranges, stage, write_mdp


def ncaa_info(ctx):
    """Where the residue's files are, and what ncaa-param recorded about it."""
    n = ctx.cfg["ncaa"]
    files = ctx.input(n["files_from"]) if "files_from" in n else ctx.other("ncaa")
    desc = json.loads((files / "ncaa_map.json").read_text())
    relax = (ctx.input(n["relax_from"]) if "relax_from" in n else ctx.other("ncaa_relax")) / f"ace_{n['name'].lower()}_nme.amb2gmx"
    return n["name"], files, relax, desc


def ncaa_reach(c, relax):
    """The residue's reach in nm: `reach` is "median", "p95" (default), "max", or a number."""
    q = c.get("reach", "p95")
    if isinstance(q, (int, float)):
        return float(q)
    return json.loads((relax / "reach.json").read_text())[{"max": "max_nm", "p95": "p95_nm", "median": "median_nm"}[q]]


def fibril_pdb(ctx, section):
    c = ctx.cfg.get(section, {})
    p = ctx.input(c["fibril"]) if "fibril" in c else ctx.other("fibril_relax") / "fibril.pdb"
    if not p.exists():
        raise StageError(f"{p} not found: run fibril-prepare and fibril-run, or give [{section}] fibril = \"your_relaxed_fibril.pdb\"")
    return p


def restraint(ctx, section, chains):
    """The production restraint as an index group: backbone of the chosen chains -> (tag, group name, selection, fc)."""
    r = ctx.cfg.get(section, {}).get("restrain", {})
    if not r.get("chains"):
        return None
    tag = "".join(r["chains"])
    atoms = r.get("atoms", "N CA C O")
    rng = [resindex_ranges(chains, ch) for ch in r["chains"]]
    sel = "(" + " or ".join(f"resindex {a} to {b}" for a, b in rng) + f") and name {atoms}"
    return tag, f"{tag}_BB", sel, float(r.get("fc", 200))


def prep_for_tleap(ctx, d, src, dst):
    drop = ctx.cfg.get("tleap", {}).get("drop", {})
    out = ctx.script("prep_fibril_tleap.py", [src, dst], dict(DROP_RESI=int(drop.get("resi", -1)), DROP_RESN=drop.get("resn", "none")))
    ctx.grab(out, re.compile(r"^wrote|may not know"))


def build_gromacs(ctx, d, out, waters):
    """acpype and its repairs, shared by the linked and the free-linker builds (protocol steps 21 and 31)."""
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


def equilibrate(ctx, g, out, group, rs, npt_mdp, density, ref_t, sec=None, where="[linked]", gate_em=True):
    """EM, NVT, NPT with every solute heavy atom restrained (protocol steps 22, 24 and 33).
    [linked] / [free] nvt_ps and npt_ps change the lengths; the defaults are the protocol's 100 ps and 1 ns."""
    tc = {"tc-grps": f"{group} Water_and_ions"}
    sec = sec or {}
    nvt_n = {"nsteps": nsteps_for("linked_nvt.mdp", float(sec["nvt_ps"]) / 1000)} if "nvt_ps" in sec else {}
    npt_n = {"nsteps": nsteps_for(npt_mdp, float(sec["npt_ps"]) / 1000)} if "npt_ps" in sec else {}
    if ctx.begin("em", "Energy minimisation"):
        write_mdp("linked_em.mdp", g / "linked_em.mdp")
        gmxmd.run_em(ctx, g, "linked_em.mdp", f"{out}_GMX.gro", "topol.top", n="index.ndx", maxwarn=1, second_pass=False)
        ctx.end()
    gsec = (where, sec)
    if gate_em:
        ctx.check_gate("em", gsec, "the minimisation", "em.png", "NVT")
    if ctx.begin("nvt", "NVT with the solute restrained"):
        write_mdp("linked_nvt.mdp", g / "linked_nvt.mdp", {**tc, **nvt_n})
        gmxmd.run_md(ctx, g, "nvt", "linked_nvt.mdp", "em.gro", "topol.top", r="em.gro", n="index.ndx", maxwarn=1)
        gmxmd.judge(ctx, g, "nvt", [group, "Water_and_ions"], ref_t)
        ctx.end()
    ctx.check_gate("nvt", gsec, "the NVT equilibration", "nvt.png", "NPT")
    if ctx.begin("npt", "NPT with the solute restrained"):
        write_mdp(npt_mdp, g / npt_mdp, {**tc, **npt_n})
        gmxmd.run_md(ctx, g, "npt", npt_mdp, "nvt.gro", "topol.top", r="nvt.gro", t="nvt.cpt", n="index.ndx", maxwarn=1)
        gmxmd.judge(ctx, g, "npt", [group, "Water_and_ions"], ref_t, gmxmd.density_range(ctx, density, g / "nvt.gro"))
        ctx.note(f"box before NPT {' x '.join(gmxmd.box_of(g / 'nvt.gro')[:3])}, after {' x '.join(gmxmd.box_of(g / 'npt.gro')[:3])} nm")
        ctx.end()
    ctx.check_gate("npt", gsec, "the NPT equilibration", "npt.png", "production")
    ctx.pause(gsec, "production")


def production(ctx, g, group, rs, mdp, ns, extra=None):
    """Production with only the chosen chains' backbone restrained (protocol steps 25 and 34)."""
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


@stage("linked-build", "linked", "18-21, 23")
def linked_build(ctx):
    """relaxed fibril + the residue at chosen sites -> solvated GROMACS system"""
    L, d = ctx.cfg.get("linked", {}), ctx.dir
    name, files, relax, desc = ncaa_info(ctx)
    out = L.get("name", "linked")
    group = L.get("group", "Solute")
    sites = [[s["chain"], int(s["resi"])] for s in ctx.cfg["attach"]["sites"]]
    graft = f"{name.lower()}.pdb"

    if ctx.begin("s18", "18 · Place the residue on the fibril (PyMOL)"):
        shutil.copy(fibril_pdb(ctx, "linked"), d / "fibril.pdb")
        shutil.copy(relax / graft, d / graft)
        o = ctx.script("place_ncaa.py", pymol=True, overrides=dict(
            FIBRIL="fibril.pdb", LINKER=graft, OUT="full_fib.pdb", RESN=name, SITES=sites, CHI=desc["chi"],
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
                             "residue whose side chain is exposed (check its SASA in the relaxed fibril) or redesign the protein. "
                             "To build it anyway, set [attach] allow_clash = true")
        ca = [(l[21], int(l[22:26]), l[17:21].strip()) for l in pdb_atoms(d / "full_fib.pdb") if l[12:16].strip() == "CA"]
        for ch, resi in sites:
            k = [i for i, x in enumerate(ca) if x[:2] == (ch, resi)]
            ok = len(k) == 1 and ca[k[0]][2] == name and (k[0] == 0 or ca[k[0] - 1][0] != ch or ca[k[0] - 1][1] < resi)
            ctx.note(f"{name} at {ch}{resi}: " + ("in its own residue position" if ok else "NOT in sequence -- tleap would bond it wrongly"))
            if not ok:
                raise StageError(f"{name} is not in place at {ch}{resi}")
        ctx.end()

    if ctx.begin("s19", "19 · AMBER names and TER records"):
        prep_for_tleap(ctx, d, "full_fib.pdb", "full_fib_noh.pdb")
        k = sum(1 for l in pdb_atoms(d / "full_fib_noh.pdb") if l[17:20] == name)
        ctx.note(f"{name} heavy atoms: {k} (expected {desc['heavy_atoms'][name]} per residue, {desc['heavy_atoms'][name] * len(sites)} for {len(sites)})")
        ctx.end()

    if ctx.begin("s20", "20 · Build, solvate and add salt in tleap"):
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
            abox, lines, ext, stretched = auto_buffer(d / "full_fib_noh.pdb", parts, float(L.get("image_margin_nm", 2.0)))
            fixed = 12.0 if box == "auto" else float(box)
            if box == "auto":
                for l in lines:
                    ctx.note("auto box: " + l)
            for l in advice("[linked] box", box, atoms_in((ext + 2 * np.array(abox)) / 10), atoms_in((ext + 2 * fixed) / 10), f"a {fixed:g} A buffer",
                            f"The attached {name} reaches {r_nm:.2f} nm from its backbone.", "the periodic image check of step 26"):
                ctx.note(l)
            if box == "auto":
                box = abox
        charge, waters, pairs = amber.tleap_build(ctx, d, name, out, "fib = loadPdb full_fib_noh.pdb\n", "fib", box, "Box",
                                                  float(L.get("salt", 0.15)), extra_save="savePdb fib fib_check.pdb\n")
        k = sum(1 for l in open(d / "fib_check.pdb") if l[17:20] == name and l.startswith(("ATOM", "HETATM")))
        ctx.note(f"fib_check.pdb: {k} {name} atoms (expected {desc['n_residue_atoms'] * len(sites)} = {len(sites)} x {desc['n_residue_atoms']})")
        shape = dict(extent_nm=[float(x) / 10 for x in ext], stretched_nm=[float(x) / 10 for x in stretched]) if r_nm is not None else {}
        (d / "build_info.json").write_text(json.dumps(dict(charge=charge, waters=waters, pairs=pairs, box=box, **shape)) + "\n")
        ctx.end()

    info = json.loads((d / "build_info.json").read_text())
    if ctx.begin("s21", "21 · Convert to GROMACS with acpype"):
        g = build_gromacs(ctx, d, out, info["waters"])
        if "extent_nm" in info:
            for l in built([float(x) for x in gmxmd.box_of(g / f"{out}_GMX.gro")[:3]], info["extent_nm"],
                           f"the fibril with its {name}", info["stretched_nm"]):
                ctx.note(l)
        ctx.end()

    g = d / f"{out}.amb2gmx"
    if ctx.begin("s23", "23 · Thermostat group, restraint files and topology includes"):
        chains = pdb_chains(d / "full_fib_noh.pdb")
        rs = restraint(ctx, "linked", chains)
        _, natoms, nres = amber.restraint_files(ctx, g, out, f"POSRES_{rs[0]}" if rs else "POSRES_NONE")
        groups = {"System": "all", group: f"resindex 1 to {nres}", "Water_and_ions": f"not resindex 1 to {nres}"}
        if rs:
            groups[rs[1]] = rs[2]
        amber.make_index(ctx, g, f"{out}_GMX.gro", groups)
        ctx.end()


@stage("linked-run", "linked", "22, 24-26")
def linked_run(ctx):
    """EM, NVT, NPT, production with the chosen chains restrained, then the periodic-image check"""
    L = ctx.cfg.get("linked", {})
    out, group = L.get("name", "linked"), L.get("group", "Solute")
    g = ctx.dir / f"{out}.amb2gmx"
    rs = restraint(ctx, "linked", pdb_chains(ctx.dir / "full_fib_noh.pdb"))
    equilibrate(ctx, g, out, group, rs, "linked_npt.mdp", L.get("density_range"), float(L.get("temperature", 298.0)), L)
    if ctx.begin("s25", "25 · Production with the chosen chains restrained"):
        production(ctx, g, group, rs, "linked_prod.mdp", float(L.get("production_ns", 200)))
        ctx.end()
    if ctx.begin("s26", "26 · Periodic image distance of the fibril with its linkers"):
        dmin = gmxmd.image_check_single(ctx, g, "prod.tpr", "prod.xtc", "index.ndx", group, int(L.get("image_check_dt_ps", 100)))
        if dmin < 2.0:
            b = json.loads((ctx.dir / "build_info.json").read_text())["box"]
            add = (2.0 - dmin) * 10 / 2 / 0.95
            ctx.note(f"to clear 2.0 nm, rebuild with at least {add:.0f} A more buffer per side (the box was built with {b})" +
                     ("; the two parts touched, so they could have reached further than this: box = \"auto\" sizes the "
                      "buffer from the residue's reach instead" if dmin < 0.5 else ""))
        ctx.end()


@stage("linked-analyse", "linked", "27")
def linked_analyse(ctx):
    """linker conformations and adsorption (the study's linked_eq_check.py; its atom groups describe CYL)"""
    L = ctx.cfg.get("linked", {})
    A = L.get("analysis", {})
    name = ctx.cfg["ncaa"]["name"]
    out, group = L.get("name", "linked"), L.get("group", "Solute")
    g = ctx.dir / f"{out}.amb2gmx"
    if ctx.begin("s27", "27 · Linker conformations and adsorption in the linked run"):
        ctx.sh(f"{ctx.gmx} trjconv -s prod.tpr -f prod.cpt -n index.ndx -o solute.gro", cwd=g, stdin=group + "\n", log="solute_gro")
        ctx.sh(f"{ctx.gmx} trjconv -s prod.tpr -f prod.xtc -n index.ndx -o solute.xtc", cwd=g, stdin=group + "\n", log="solute_xtc")
        chains = pdb_chains(ctx.dir / "full_fib_noh.pdb")
        names, short = {}, {}
        for s in ctx.cfg["attach"]["sites"]:
            k = resindex_ranges(chains, f"{s['chain']}:{s['resi']}-{s['resi']}")[0]
            names[f"{name}{k}"] = s.get("label", f"chain {s['chain']}") + f", residue {s['resi']}"
            short[f"{name}{k}"] = s.get("label", f"chain {s['chain']}")
        skip = {name, "WAT", "SOL", "HOH", "Na+", "Cl-", "NA", "CL", "K+", "K"}
        args = ["solute.gro", "solute.xtc"] + ([str(A["every_ps"])] if "every_ps" in A else [])
        o = ctx.script(A.get("script", "linked_eq_check.py"), args, dict(SITE_NAMES=names, SITE_SHORT=short, SKIP_RES=f"py:{skip!r}"), cwd=g)
        (g / "linked_eq_report.txt").write_text(o)
        gmxmd.draw(ctx, g, A.get("plot_script", "plot_linked_eq.py"))
        ctx.note(f"sites: {names}")
        ctx.grab(o, re.compile(r"lying on it|latest of these|for the bonded distributions|biotin from"))
        ctx.note(f"full report: {g / 'linked_eq_report.txt'}; figure and table: linked_eq.png, linked_eq.csv")
        ctx.end()
