"""The non-canonical residue: parameters (steps 5-11) and its run in water (step 12)."""
import json
import os
import re
import shutil
from pathlib import Path

from .. import amber, gmxmd
from ..core import StageError, nsteps_for, stage, temp_changes, write_mdp
from ..ncaa_map import MapError, retype, write_mc

SQM = "qm_theory='AM1', grms_tol=0.0005, scfconv=1.d-10, ndiis_attempts=700, maxcyc=0"
CHARGE_METHODS = {"abcg2": "ABCG2", "bcc": "AM1-BCC", "mul": "Mulliken (AM1)", "gas": "Gasteiger",
                  "resp": "RESP", "esp": "ESP (Merz-Kollman)", "rc": "read-in"}      # the antechamber -c methods
SQM_METHODS = ("abcg2", "bcc", "mul")                   # antechamber runs sqm for these, with sqm_keywords
QM_FILE = ("resp", "esp")                               # these fit charges to a QM run done beforehand
QM_FORMATS = {".log": "gout", ".out": "gout", ".gesp": "gesp", ".dat": "gamess"}
NOT_METHODS = {
    "cm1": "antechamber lists it, but sqm cannot compute CM1 charges and antechamber stops",
    "cm2": "antechamber lists it, but sqm cannot compute CM2 charges and antechamber stops",
    "wc": "wc writes charges to a file, it does not compute them",
    "dc": "dc deletes the charges (all zero)",
}


def charge_method(n):
    """The [ncaa] charge_method, checked, with the settings it needs."""
    cm = str(n.get("charge_method", "abcg2")).lower()
    if cm in NOT_METHODS:
        raise StageError(f"[ncaa] charge_method = \"{cm}\": {NOT_METHODS[cm]}. Choose one of {', '.join(CHARGE_METHODS)}")
    if cm not in CHARGE_METHODS:
        raise StageError(f"[ncaa] charge_method = \"{cm}\" is not an antechamber charge method. "
                         f"Choose one of {', '.join(CHARGE_METHODS)}")
    if cm in QM_FILE and not n.get("qm_output"):
        raise StageError(f"[ncaa] charge_method = \"{cm}\" fits the charges to the electrostatic potential of a QM run: "
                         "set qm_output to its Gaussian output (.log, .out or .gesp) or GAMESS .dat file, "
                         "run on the capped molecule with its atoms in the order of the mol2")
    if cm == "rc" and not n.get("charge_file"):
        raise StageError("[ncaa] charge_method = \"rc\" reads the charges from a file: set charge_file "
                         "(one charge per atom in the order of the mol2, as antechamber -c wc writes them)")
    return cm


def read_ac(path):
    """[(element, charge), ...] of an antechamber .ac file, in atom order (the element read from the GAFF type)."""
    rows = []
    for l in open(path):
        if l.startswith("ATOM"):
            t = l.split()[-1].lower()
            rows.append((t[:2] if t[:2] in ("cl", "br") else t[0], float(l[54:64])))
    return rows


def write_charges(path, q):
    """A charge file antechamber -c rc reads: eight charges a line."""
    path.write_text("\n".join("".join(f"{c:10.6f}" for c in q[i:i + 8]) for i in range(0, len(q), 8)) + "\n")


def need(d, *files, step):
    """Stop if a step left one of its files missing or empty: a program can exit quietly without writing it."""
    bad = [f for f in files if not (d / f).exists() or (d / f).stat().st_size == 0]
    if bad:
        raise StageError(f"step {step} did not write {', '.join(bad)}; see the newest files in {d / 'logs'}")


def read_prepin(path):
    """[(name, type, charge), ...] of a prepin's real atoms ."""
    rows = []
    for l in open(path).read().split("\n")[7:]:
        f = l.split()
        if len(f) > 3 and f[0].isdigit() and f[1] != "DUMM":
            rows.append((f[1], f[2], float(f[-1])))
    return rows


def parent_cmap(frcmod, parent):
    """The CMAP grid ff19SB assigns to a residue, as text (from the CMAP section of frcmod.ff19SB)."""
    text = open(frcmod).read()
    for block in text[text.index("\nCMAP"):].split("%FLAG CMAP_COUNT")[1:]:
        names = block.split("%FLAG CMAP_RESLIST")[1].split("\n")[1].split()
        if parent in names:
            return block.split("%FLAG CMAP_PARAMETER")[1].lstrip("\n").rstrip() + "\n"
    raise StageError(f"ff19SB has no CMAP listed for {parent}")


@stage("ncaa-param", "ncaa", "5-11")
def ncaa_param(ctx):
    """capped NCAA .mol2 -> XXX.prepin, frcmod.XXX, frcmod.ff19SB_XXX, checked against the parent residue"""
    n = ctx.cfg["ncaa"]
    name, parent, nc, d = n["name"], n.get("parent", "CYX"), int(n.get("net_charge", 0)), ctx.dir
    if len(name) != 3 or name in ("CYS", "CYX", "CYM"):
        raise StageError("[ncaa] name must be three characters and not an existing residue name")

    cq_file = d / "charge_choice.txt"
    cq = str(n.get("charge_method", "abcg2")).lower()
    old = cq_file.read_text().strip() if cq_file.exists() else "bcc"     # runs from before the setting was saved used bcc
    if "s5" in ctx.done and old != cq:
        for k in ("s5", "s6", "s7", "s8", "s9", "s10", "s11"):          # the charges changed: every step after depends on them
            ctx.done.pop(k, None)
        ctx.say(f"   [ncaa] charge_method is now \"{cq}\" (was \"{old}\"): steps 5 to 11 are redone; "
                "rerun ncaa-relax and the stages after it too")
    if ctx.begin("s5", "5 · antechamber with GAFF2 types"):
        src = ctx.input(n.get("mol2") or n["structure"])
        fmt = src.suffix.lstrip(".").lower()                  # mol2, sdf or pdb: antechamber reads all three
        shutil.copy(src, d / f"ncaa.{fmt}")
        cm = charge_method(n)
        if cm in QM_FILE:                                     # RESP/ESP: charges from the user's QM run, names from the mol2
            qm = ctx.input(n["qm_output"])
            if not qm.exists():
                raise StageError(f"[ncaa] qm_output: {qm} does not exist")
            qfmt = n.get("qm_format") or QM_FORMATS.get(qm.suffix.lower())
            if qfmt not in QM_FORMATS.values():
                raise StageError(f"[ncaa] qm_format: give one of {', '.join(sorted(set(QM_FORMATS.values())))} for {qm.name}")
            shutil.copy(qm, d / f"qm{qm.suffix}")
            ctx.sh(f"antechamber -i qm{qm.suffix} -fi {qfmt} -o qm_charges.ac -fo ac -c {cm} -nc {nc} -pf y", log="antechamber_qm")
            need(d, "qm_charges.ac", step=5)
            ctx.sh(f"antechamber -i ncaa.{fmt} -fi {fmt} -o ncaa_order.ac -fo ac -c dc -pf y", log="antechamber_order")
            need(d, "ncaa_order.ac", step=5)
            q, el = read_ac(d / "qm_charges.ac"), read_ac(d / "ncaa_order.ac")
            if [e for e, _ in q] != [e for e, _ in el]:
                raise StageError(f"{qm.name} does not list the atoms of {src.name} in the same order "
                                 f"({len(q)} vs {len(el)} atoms): run the QM job on the molecule as drawn in {src.name}")
            write_charges(d / "charges.crg", [c for _, c in q])
            opts = "-c rc -cf charges.crg"
        elif cm == "rc":
            cf = ctx.input(n["charge_file"])
            if not cf.exists():
                raise StageError(f"[ncaa] charge_file: {cf} does not exist")
            shutil.copy(cf, d / "charges.crg")
            opts = "-c rc -cf charges.crg"
        else:
            opts = f"-c {cm}" + (f' -ek "{n.get("sqm_keywords", SQM)}"' if cm in SQM_METHODS else "")
        ctx.sh(f'antechamber -i ncaa.{fmt} -fi {fmt} -o {name}_capped.ac -fo ac -rn {name} '
               f'-nc {nc} -at gaff2 -pf y {opts}', log="antechamber")
        need(d, f"{name}_capped.ac", step=5)
        k = sum(1 for l in open(d / f"{name}_capped.ac") if l.startswith("ATOM"))
        ctx.note(f"{name}_capped.ac: {k} atoms, {CHARGE_METHODS[cm]} charges, net charge {nc}")
        cq_file.write_text(cm + "\n")
        ctx.end()

    if ctx.begin("s6", "6 · Names, types and charges (found in the bond graph, no atom numbers)"):
        try:
            out, desc = retype(d / f"{name}_capped.ac", d / f"{name}_fixed.ac", parent, nc,
                               n.get("set_backbone_charges", True), n.get("parent_atoms"), d / "ncaa_map.json", name)
        except MapError as e:
            raise StageError(str(e))
        (d / "logs" / "s6_atom_map.txt").write_text("\n".join(out) + "\n")
        ctx.note(f"parent {parent}: backbone N H CA HA C O and side-chain atoms {' '.join(desc['parent_side_names'])} put on ff19SB types")
        ctx.note("grafting torsions: " + "; ".join("-".join(c) for c in desc["chi"]))
        ctx.note(out[-1])
        ctx.end()
    desc = json.loads((d / "ncaa_map.json").read_text())

    if ctx.begin("s7", "7 · The main-chain file and prepgen"):
        write_mc(d / f"{name}.mc", desc)
        ctx.sh(f"prepgen -i {name}_fixed.ac -o {name}.prepin -m {name}.mc -rn {name}", log="prepgen")
        need(d, f"{name}.prepin", step=7)
        rows = read_prepin(d / f"{name}.prepin")
        q = sum(r[2] for r in rows)
        ctx.note(f"{name}.prepin: {len(rows)} atoms, total charge {q:.4f}")
        ctx.note("types: " + "  ".join(f"{a} {t}" for a, t, _ in rows if a in desc["parent_types"]))
        if len(rows) != desc["n_residue_atoms"] or abs(q - nc) > 1e-3:
            raise StageError(f"prepin check failed: expected {desc['n_residue_atoms']} atoms and charge {nc}")
        ctx.end()

    if ctx.begin("s8", "8 · parmchk2 against GAFF2 and ff19SB, then merge_frcmod.py"):
        P = Path(os.environ["AMBERHOME"]) / "dat" / "leap" / "parm"
        (d / "ff19SB_XC.atc").write_text("PARM    XC      0               0       12.01   0       6\nEQUA    CX\nEQUA    CT\nEQUA    c3\n")
        keep = []
        for l in (P / "frcmod.ff19SB").read_text().split("\n"):
            if l.startswith("CMAP"):
                break
            keep.append("IMPROPER" if re.fullmatch(r"IMPR *", l) else "NONBON" if re.fullmatch(r"NONB *", l) else l)
        (d / "frcmod.ff19SB.parmchk2").write_text("\n".join(keep) + "\n")
        ctx.sh(f"parmchk2 -i {name}_fixed.ac -f ac -o frcmod.{name}.gaff2 -s 2 -atc ff19SB_XC.atc", log="parmchk2_gaff2")
        ctx.sh(f"parmchk2 -i {name}_fixed.ac -f ac -o frcmod.{name}.ff19SB -p {P / 'parm19.dat'} "
               f"-afrc frcmod.ff19SB.parmchk2 -atc ff19SB_XC.atc", log="parmchk2_ff19SB")
        need(d, f"frcmod.{name}.gaff2", f"frcmod.{name}.ff19SB", step=8)
        out = ctx.script("merge_frcmod.py", ["--primary", f"frcmod.{name}.gaff2", "--fallback", f"frcmod.{name}.ff19SB",
                                             "-o", f"frcmod.{name}", "--flag-above", str(n.get("penalty_flag_above", 10.0))])
        need(d, f"frcmod.{name}", step=8)
        ctx.note(out.strip().split("\n")[-1])
        ctx.end()

    cm_file = d / "cmap_choice.txt"
    if "s9" in ctx.done and cm_file.exists() and cm_file.read_text().strip() != n.get("cmap", "generic"):
        for k in ("s9", "s10", "s11"):          # the cmap setting changed: these three steps depend on it
            ctx.done.pop(k, None)
        ctx.say(f"   [ncaa] cmap is now \"{n.get('cmap', 'generic')}\": steps 9 to 11 are redone")
    if ctx.begin("s9", "9 · CMAP for the new residue"):
        P = Path(os.environ["AMBERHOME"]) / "dat" / "leap" / "parm"
        tpl = (P / "frcmod.ff19SB_XXX").read_text().replace("XXX", name)
        if n.get("cmap", "generic") == "parent":
            grid = parent_cmap(P / "frcmod.ff19SB", parent)
            tpl = tpl[:tpl.index("%FLAG CMAP_PARAMETER")] + "%FLAG CMAP_PARAMETER\n" + grid
            ctx.note(f"frcmod.ff19SB_{name}: the CMAP ff19SB gives {parent}, under the name {name}")
        else:
            ctx.note(f"frcmod.ff19SB_{name}: the ff19SB generic (Leu) CMAP under the name {name}")
        (d / f"frcmod.ff19SB_{name}").write_text(tpl)
        cm_file.write_text(n.get("cmap", "generic") + "\n")
        ctx.end()

    if ctx.begin("s10", "10 · tleap: test peptides for checking"):
        (d / "check.leap").write_text(amber.LEAP_HEAD.format(name=name) +
                                      f"p = sequence {{ ACE ALA {name} ALA NME }}\ncheck p\ncharge p\n"
                                      f"saveAmberParm p test_{name}.parm7 test_{name}.rst7\n"
                                      f"q = sequence {{ ACE ALA {parent} ALA NME }}\n"
                                      f"saveAmberParm q test_{parent}.parm7 test_{parent}.rst7\nquit\n")
        log = ctx.sh("tleap -f check.leap", log="tleap_check")
        (d / "check_leap.log").write_text(log)
        bad = [l for l in log.split("\n") if amber.LEAP_BAD.search(l)]
        m = re.search(r"Total unperturbed charge:\s*([-\d.]+)", log)
        if bad or not m:
            raise StageError("tleap could not build the test peptide:\n" + "\n".join(bad[:8]))
        need(d, f"test_{name}.parm7", f"test_{parent}.parm7", step=10)
        ctx.note(f"ACE-ALA-{name}-ALA-NME: total charge {float(m.group(1)):.6f}; no FATAL, no missing type or parameter")
        ctx.end()

    if ctx.begin("s11", "11 · Check the prmtop"):
        out = ctx.script("check_ff19sb_prmtop.py", [f"test_{name}.parm7", "--res", name, "--ref", f"test_{parent}.parm7",
                                                    "--ref-res", parent], check=False)
        (d / "check_prmtop.txt").write_text(out)
        last = ([l for l in out.split("\n") if "PROBLEM" in l or "PASSED" in l] or [out.strip().split("\n")[-1]])[-1].strip()
        ctx.note(last)
        own_map = f"CMAP grid identical to {parent}: False" in out and "1 PROBLEM(S) FOUND" in out and n.get("cmap", "generic") == "generic"
        if own_map and "cmap" in n.get("accept", []):
            ctx.note(f"the one difference is the CMAP itself: ff19SB gives {parent} its own map and {name} carries the generic one; "
                     "accepted in the config ([ncaa] accept)")
        elif own_map:
            raise StageError(
                f"step 11 found one difference from {parent}, the CMAP: ff19SB gives {parent} its own backbone map, and {name} was "
                f"given the generic one (cmap = \"generic\").\n  Choose one under [ncaa] in the config, then run the same command again:\n"
                f"    cmap = \"parent\"      {name} gets {parent}'s map; steps 9 to 11 are redone\n"
                f"    accept = [\"cmap\"]    keep the generic map on purpose\n  Full report: {d / 'check_prmtop.txt'}")
        elif "ALL CHECKS PASSED" not in out:
            raise StageError(f"step 11 (check_ff19sb_prmtop.py) did not pass: {last}\n  Full report: {d / 'check_prmtop.txt'}")
        ctx.end()


def extended_reach(prmtop, resname):
    """The farthest any heavy atom of the residue can be from its CA, in nm, with the chain pulled straight.
    Along the bonded path from CA to each atom, every pair of bonds spans at most the distance its equilibrium
    angle allows; the sum over the path is the all-trans length, which no conformation exceeds (bond and angle
    vibration aside). It needs no simulation, so it is the same every time."""
    import math
    import networkx as nx
    import parmed as pmd
    parm = pmd.load_file(str(prmtop))
    res = [r for r in parm.residues if r.name == resname][0]
    idx = {a.idx for a in res.atoms}
    g = nx.Graph()
    for b in parm.bonds:
        if b.atom1.idx in idx and b.atom2.idx in idx:
            g.add_edge(b.atom1.idx, b.atom2.idx, w=b.type.req)
    theta = {}
    for an in parm.angles:
        theta[(an.atom1.idx, an.atom2.idx, an.atom3.idx)] = theta[(an.atom3.idx, an.atom2.idx, an.atom1.idx)] = an.type.theteq
    ca = [a.idx for a in res.atoms if a.name == "CA"][0]
    paths = nx.single_source_dijkstra_path(g, ca, weight="w")
    best = 0.0
    for a in res.atoms:
        if a.atomic_number == 1 or a.idx not in paths:
            continue
        p, d, i = paths[a.idx], 0.0, 0
        while i + 2 < len(p):
            b1, b2 = g[p[i]][p[i + 1]]["w"], g[p[i + 1]][p[i + 2]]["w"]
            t = math.radians(theta.get((p[i], p[i + 1], p[i + 2]), 180.0))
            d += math.sqrt(b1 * b1 + b2 * b2 - 2 * b1 * b2 * math.cos(t))
            i += 2
        if i + 1 < len(p):
            d += g[p[i]][p[i + 1]]["w"]
        best = max(best, d)
    return best / 10.0


@stage("ncaa-relax", "ncaa_relax", "12")
def ncaa_relax(ctx):
    """the three files -> capped NCAA in water: a frame to graft, frames to scatter, and its reach"""
    n = ctx.cfg["ncaa"]
    r = n.get("relax", {})
    name, d = n["name"], ctx.dir
    out = f"ace_{name.lower()}_nme"
    g = d / f"{out}.amb2gmx"
    group = f"ACE_{name}_NME"
    flags = r.get("mdrun_flags", "")
    T = r.get("temperature", n.get("temperature"))            # None = the mdp files' 298 K
    sec = ("[ncaa.relax]", r)
    shape = {"oct": "Oct", "box": "Box"}.get(str(r.get("box_shape", "oct")).lower())
    if shape is None:
        raise StageError('[ncaa.relax] box_shape has to be "oct" (truncated octahedron) or "box" (rectangular)')
    ns = float(r.get("production_ns", 10.0))                   # long enough for the 95% reach to settle (see s12g)
    wanted = [float(r.get("graft_frame_ps", 1000))] + [float(t) for t in r.get("frames_ps", [250, 500, 750, 1000])]
    if max(wanted) > ns * 1000 + 1e-6:
        raise StageError(f"[ncaa.relax] asks for a frame at {max(wanted):g} ps, but production_ns = {ns:g} ends at {ns * 1000:g} ps: "
                         "lengthen the run or take earlier frames (graft_frame_ps, frames_ps)")

    if ctx.begin("s12a", "12 · NCAA coordinates from an ACE-X-NME run: tleap build"):
        amber.copy_ncaa_files(ctx, ctx.input(n["files_from"]) if "files_from" in n else ctx.other("ncaa"), d, name)
        amber.tleap_build(ctx, d, name, out, f"m = sequence {{ ACE {name} NME }}\n", "m", r.get("box", 12.0), shape,
                          r.get("salt", n.get("salt", 0.15)), r.get("ion_pairs"))
        ctx.end()

    if ctx.begin("s12b", "12 · Convert to GROMACS with acpype"):
        amber.acpype(ctx, d, out)
        ctx.end()

    top, gro = f"{out}_GMX.top", f"{out}_GMX.gro"
    if ctx.begin("s12c", "12 · minimisation"):
        write_mdp("ncaa_em.mdp", g / "em.mdp")
        em = gmxmd.run_em(ctx, g, "em.mdp", gro, top)
        (g / "em_final.txt").write_text(em + "\n")
        ctx.end()
    em = (g / "em_final.txt").read_text().strip()
    ctx.check_gate("em", sec, "the minimisation (after its second pass, if one was needed)", "em.png", "NPT")

    if ctx.begin("s12d", "12 · NPT equilibration"):
        ch = temp_changes("ncaa_eq.mdp", T)
        if "npt_ps" in r:
            ch["nsteps"] = nsteps_for("ncaa_eq.mdp", float(r["npt_ps"]) / 1000)
        write_mdp("ncaa_eq.mdp", g / "eq.mdp", ch)
        ctx.note(f"NPT: {float(r.get('npt_ps', 200)):g} ps at {float(T if T is not None else 298):g} K")
        gmxmd.run_md(ctx, g, "npt", "eq.mdp", em, top, flags=flags)
        gmxmd.judge(ctx, g, "npt", ["System"], float(T if T is not None else 298.0), table=r, flags=flags)
        ctx.end()

    ctx.check_gate("npt", sec, "the NPT equilibration of the residue in water", "npt.png", "production")
    ctx.pause(sec, "production")
    if ctx.begin("s12e", "12 · production"):
        write_mdp("ncaa_prod.mdp", g / "prod.mdp", {"nsteps": nsteps_for("ncaa_prod.mdp", ns), **temp_changes("ncaa_prod.mdp", T)})
        gmxmd.run_md(ctx, g, "prod", "prod.mdp", "npt.gro", top, t="npt.cpt", flags=flags)
        ctx.end()

    if ctx.begin("s12f", "12 · export the frame to graft and the frames to scatter"):
        ctx.sh(f"""{ctx.gmx} select -s prod.tpr -on ncaa.ndx -select '"{group}" resname ACE {name} NME'""", cwd=g, log="select")
        graft = f"{name.lower()}.pdb"
        t0 = r.get("graft_frame_ps", 1000)
        ctx.sh(f"{ctx.gmx} trjconv -s prod.tpr -f prod.xtc -n ncaa.ndx -o {graft} -dump {t0} -pbc mol -conect", cwd=g, stdin=group + "\n", log="graft_frame")
        frames = r.get("frames_ps", [250, 500, 750, 1000])
        for t in frames:
            ctx.sh(f"{ctx.gmx} trjconv -s prod.tpr -f prod.xtc -n ncaa.ndx -o ncaa_{t}.pdb -dump {t} -pbc mol", cwd=g, stdin=group + "\n", log=f"frame_{t}")
        ctx.note(f"{graft}: the frame at {t0} ps, whole molecule; ncaa_<t>.pdb at {', '.join(str(t) for t in frames)} ps")
        ctx.end()

    if ctx.begin("s12g", "12 · reach of the residue in water (for box = \"auto\")"):
        import MDAnalysis as mda
        import numpy as np
        ctx.sh(f"{ctx.gmx} trjconv -s prod.tpr -f prod.xtc -n ncaa.ndx -o mol.xtc -pbc mol", cwd=g, stdin=group + "\n", log="mol_xtc")
        ctx.sh(f"{ctx.gmx} trjconv -s prod.tpr -f prod.xtc -n ncaa.ndx -o mol.gro -dump 0 -pbc mol", cwd=g, stdin=group + "\n", log="mol_gro")
        u = mda.Universe(str(g / "mol.gro"), str(g / "mol.xtc"))
        ca = u.select_atoms(f"resname {name} and name CA")
        heavy = u.select_atoms(f"resname {name} and not name H*")
        reach = np.array([np.linalg.norm(heavy.positions - ca.positions[0], axis=1).max() for _ in u.trajectory]) / 10.0
        res = dict(frames=int(len(reach)), run_ns=ns, max_nm=float(reach.max()), p95_nm=float(np.percentile(reach, 95)),
                   median_nm=float(np.median(reach)), note="largest CA-to-heavy-atom distance of the residue, per frame")
        res["extended_nm"] = extended_reach(d / f"{out}.prmtop", name)
        (g / "reach.json").write_text(json.dumps(res, indent=1) + "\n")
        ctx.note(f"reach (CA to the farthest heavy atom): median {res['median_nm']:.2f} nm, 95% {res['p95_nm']:.2f} nm, "
                 f"largest {res['max_nm']:.2f} nm over {res['frames']} frames")
        ctx.note(f"fully stretched (all-trans path from CA, no simulation needed): {res['extended_nm']:.2f} nm; reach = \"extended\" "
                 "sizes the automatic box from this, the same every time")
        ctx.note("the reach is used only by the automatic box of linked-build and free-build. With a box set by hand it is not used, "
                 "and this run then only has to be long enough for the frames exported above")
        if ns < 5:
            ctx.note(f"this run is {ns:g} ns: the median and 95% reach of a flexible residue are not settled in so short a run, "
                     "and the largest reach varies more still. Run 10 ns or more if the automatic box will use it")
        ctx.end()
