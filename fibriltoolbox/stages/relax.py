"""The protein-only fibril in water: pdb2gmx, box, equilibration, production and its checks."""
import re
import shutil
from pathlib import Path

from .. import gmxmd
from ..box import advice, atoms_in, auto_lengths, built
from ..core import StageError, chain_breaks, nsteps_for, restraint_tag, temp_changes, pdb_atoms, pdb_chains, resindex_ranges, stage, write_mdp

GROUPS = ["Protein", "Water_and_ions"]


def _structure(ctx):
    c = ctx.cfg.get("relax", {})
    if "structure" not in c:
        raise StageError("no structure to relax: give [relax] structure = \"your_protein.pdb\" (a protein-only model of the whole assembly)")
    p = ctx.input(c["structure"])
    if not p.exists():
        raise StageError(f"[relax] structure: {p} not found")
    return p


def _chain_posre(d, spec, atoms, fc):
    """Position restraints for one chain, or a stretch of it ("G:226-282"), numbered inside its own moleculetype
    (pdb2gmx writes one .itp per chain). Returns the chain's .itp, the restraint file's name and the atom count."""
    chain, _, rng = str(spec).partition(":")
    lo, hi = (int(x) for x in rng.split("-")) if rng else (None, None)
    itp = d / f"topol_Protein_chain_{chain}.itp"
    if not itp.exists():
        raise StageError(f"{itp.name} not found: pdb2gmx writes one per chain when the structure has several chains")
    idx, sect = [], None
    for l in itp.read_text().split("\n"):
        s = l.strip()
        if s.startswith("["):
            sect = s.strip("[] ").lower()
        elif sect == "atoms" and s and not s.startswith(";") and s.split()[0].isdigit() and s.split()[4] in atoms:
            if lo is None or lo <= int(s.split()[2]) <= hi:
                idx.append(int(s.split()[0]))
    if not idx:
        raise StageError(f"[relax] restrain: no atoms {sorted(atoms)} found for {spec} (chain {chain}" + (f", residues {lo}-{hi})" if rng else ")"))
    name = f"posre_{restraint_tag([spec])}.itp"
    (d / name).write_text("[ position_restraints ]\n; atom  type      fx      fy      fz\n" +
                          "".join(f"{i:6d}     1  {fc:g}  {fc:g}  {fc:g}\n" for i in idx))
    return itp, name, len(idx)


def _aligned(gro, chains, flexible):
    """Coordinates (nm) of the aligned structure, and the C-alpha of each flexible part's anchor residue.
    flexible: [{chain = "L", anchor = 20, reach_nm = 2.5}, ...]; residues are counted as pdb2gmx kept them."""
    L = open(gro).readlines()
    rows = L[2:2 + int(L[1])]
    xyz = [[float(l[20 + 8 * i:28 + 8 * i]) for i in range(3)] for l in rows]
    resindex, prev = [], None
    for l in rows:
        if l[:10] != prev:
            prev = l[:10]
            resindex.append((resindex[-1] if resindex else 0) + 1)
        else:
            resindex.append(resindex[-1])
    parts = []
    for f in flexible:
        k = resindex_ranges(chains, f"{f['chain']}:{int(f['anchor'])}-{int(f['anchor'])}")[0]
        ca = [x for l, x, r in zip(rows, xyz, resindex) if r == k and l[10:15].strip() == "CA"]
        if not ca:
            raise StageError(f"no CA for chain {f['chain']} residue {f['anchor']} in {gro}")
        parts.append((ca[0], float(f["reach_nm"]), f"{f['chain']}{f['anchor']}"))
    return xyz, parts


@stage("fibril-prepare", "fibril_relax", "7 (system)")
def fibril_prepare(ctx):
    """protein-only structure -> topology, box, water and ions: the system ready to minimise"""
    c, d = ctx.cfg.get("relax", {}), ctx.dir
    src = _structure(ctx)

    if ctx.begin("s7a", "7 · Solvate the fibril: pdb2gmx (ff19SB, OPC), box, water, 0.15 M NaCl"):
        shutil.copy(src, d / src.name)
        out = ctx.sh(f"{ctx.gmx} pdb2gmx -f {src.name} -o prot.gro -p topol.top -ignh -ff amber19sb -water opc", log="pdb2gmx")
        ss = [l.strip() for l in out.split("\n") if l.startswith("Linking")]
        ctx.note(f"pdb2gmx: {len(ss)} disulfide bonds made" + (": " + "; ".join(ss[:4]) if ss else " (every cysteine stays reduced)"))
        ctx.sh(f"{ctx.gmx} editconf -f prot.gro -o aligned.gro -princ -c", stdin="Protein\nSystem\n", log="editconf_princ")
        box_d = c.get("box_d", "auto")
        flex = c.get("flexible", [])
        lengths, lines, ext, stretched = auto_lengths(*_aligned(d / "aligned.gro", pdb_chains(d / src.name), flex), float(c.get("image_margin_nm", 3.0)))
        fixed = 1.2 if box_d == "auto" else float(box_d)
        if box_d == "auto":
            for l in lines:
                ctx.note("auto box: " + l)
            ctx.sh(f"{ctx.gmx} editconf -f aligned.gro -o box.gro -c -box {' '.join(f'{x:.3f}' for x in lengths)}", log="editconf_box")
        else:
            ctx.sh(f"{ctx.gmx} editconf -f aligned.gro -o box.gro -c -d {fixed} -bt triclinic", log="editconf_box")
        for l in advice("[relax] box_d", box_d, atoms_in(lengths, 130), atoms_in(ext + 2 * fixed, 130), f"{fixed * 10:g} A of water on each side",
                        "A protein tail or loop that moves is not held back by anything" +
                        ("." if flex else "; none is listed under flexible."), "the periodic image check of step 8"):
            ctx.note(l)
        for l in built([float(x) for x in gmxmd.box_of(d / "box.gro")[:3]], ext, "the protein", stretched if flex else None):
            ctx.note(l)
        ctx.grab(ctx.script("check_box2.py", ["box.gro"]), re.compile(r"PASS|FAIL"))
        tip4p = Path(shutil.which(ctx.gmx)).resolve().parent.parent / "share" / "gromacs" / "top" / "tip4p.gro"
        ctx.sh(f"{ctx.gmx} solvate -cp box.gro -cs {tip4p} -o solv.gro -p topol.top", log="solvate")
        (d / "ions.mdp").write_text("integrator = steep\nnsteps = 0\n")
        ctx.sh(f"{ctx.gmx} grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1", log="grompp_ions")
        out = ctx.sh(f"{ctx.gmx} genion -s ions.tpr -o ionised.gro -p topol.top -pname NA -nname CL -neutral "
                     f"-conc {float(c.get('salt', 0.15))}", stdin="SOL\n", log="genion")
        ctx.grab(out, re.compile(r"Will try to add"))
        ctx.note(f"ionised.gro: {int(open(d / 'ionised.gro').readlines()[1])} atoms, box {' x '.join(gmxmd.box_of(d / 'ionised.gro')[:3])} nm")
        ctx.end()



@stage("fibril-run", "fibril_relax", "7 (runs), 8-9")
def fibril_run(ctx):
    """minimisation, NVT, NPT and production, each going on only if the one before passed its checks; then the
    periodic-image and stability checks and the relaxed fibril.pdb"""
    c, d = ctx.cfg.get("relax", {}), ctx.dir
    src = _structure(ctx)
    T = c.get("temperature")                                   # None = the mdp files' 298 K
    ref_t = float(T if T is not None else 298.0)
    sec = ("[relax]", c)
    if not (d / "ionised.gro").exists():
        raise StageError("no system to run: run fibril-prepare first")

    if ctx.begin("s7b", "7 · Energy minimisation"):
        write_mdp("fibril_em.mdp", d / "em.mdp")
        gmxmd.run_em(ctx, d, "em.mdp", "ionised.gro", "topol.top", second_pass=False)
        ctx.end()

    ctx.check_gate("em", sec, "the minimisation", "em.png", "NVT")
    if ctx.begin("s7c", "7 · NVT with every protein heavy atom restrained"):
        write_mdp("fibril_nvt.mdp", d / "nvt.mdp", temp_changes("fibril_nvt.mdp", T))
        gmxmd.run_md(ctx, d, "nvt", "nvt.mdp", "em.gro", "topol.top", r="em.gro")
        gmxmd.judge(ctx, d, "nvt", GROUPS, ref_t, table=c)
        ctx.end()

    ctx.check_gate("nvt", sec, "the NVT equilibration", "nvt.png", "NPT")
    if ctx.begin("s7d", "7 · NPT with every protein heavy atom restrained"):
        write_mdp("fibril_npt.mdp", d / "npt.mdp", temp_changes("fibril_npt.mdp", T))
        gmxmd.run_md(ctx, d, "npt", "npt.mdp", "nvt.gro", "topol.top", r="nvt.gro", t="nvt.cpt")
        gmxmd.judge(ctx, d, "npt", GROUPS, ref_t, gmxmd.density_range(ctx, c.get("density_range"), d / "nvt.gro"), table=c)
        ctx.end()

    ctx.check_gate("npt", sec, "the NPT equilibration", "npt.png", "production")
    ctx.pause(sec, "production")
    rs = c.get("restrain", {})
    chains = list(rs.get("chains", []))
    tag = restraint_tag(chains)
    if ctx.begin("s7e", "7 · Production, with the chosen chains (or stretches of them) held"):
        for ch in chains:
            itp, name, k = _chain_posre(d, ch, set(rs.get("atoms", ["CA"])), float(rs.get("fc", 1000)))
            with open(itp, "a") as f:
                f.write(f'\n; added by fibriltoolbox: production restraint\n#ifdef POSRES_{tag}\n#include "{name}"\n#endif\n')
            ctx.note(f"{name}: {k} atoms of {ch} at {float(rs.get('fc', 1000)):g} kJ/mol/nm2, behind POSRES_{tag}")
        ns = float(c.get("production_ns", 200))
        write_mdp("fibril_prod.mdp", d / "prod.mdp", {"nsteps": nsteps_for("fibril_prod.mdp", ns),
                                                     "define": f"-DPOSRES_{tag}" if chains else None,
                                                     **temp_changes("fibril_prod.mdp", T)})
        gmxmd.run_md(ctx, d, "prod", "prod.mdp", "npt.gro", "topol.top", r="npt.gro", t="npt.cpt")
        ctx.note(f"production: {ns:g} ns")
        ctx.end()

    if ctx.begin("s8", "8 · Check the run: whole assembly, periodic image distance, stability"):
        ctx.sh(f"{ctx.gmx} trjconv -s prod.tpr -f prod.xtc -o cluster_prod.xtc -pbc cluster", stdin="Protein\nSystem\n", log="trjconv_cluster")
        ctx.sh(f"{ctx.gmx} trjconv -s prod.tpr -f cluster_prod.xtc -o nojump_prod.xtc -pbc nojump", stdin="System\n", log="trjconv_nojump")
        ctx.sh(f"{ctx.gmx} trjconv -s prod.tpr -f prod.cpt -o prod_so_far.gro", stdin="System\n", log="trjconv_gro")
        # one frame every 100 ps: on every frame of a long run this measurement takes far longer than the rest of the stage
        ctx.sh(f"{ctx.gmx} mindist -s prod.tpr -f nojump_prod.xtc -pi -od mindist_pi.xvg -dt {int(c.get('image_check_dt_ps', 100))}",
               stdin="Protein\n", log="mindist")
        _, text, verdict = gmxmd.image_verdict(d / "mindist_pi.xvg")
        ctx.note(text)
        ctx.note(verdict)
        clean = gmxmd.image_clean_stretch(d / "mindist_pi.xvg")
        if clean:
            ctx.note(clean)
        # the general check: did the assembly hold together and stop changing shape? Any protein, no settings needed
        if c.get("assembly_check", True):
            out = ctx.script("check_assembly_stability.py", ["prod_so_far.gro", "nojump_prod.xtc", str(c.get("image_check_dt_ps", 100))],
                             {k2: c[k1] for k1, k2 in (("contacts_kept_min", "Q_MIN"), ("stability_window", "EQ_FRACTION")) if k1 in c})
            (d / "assembly_stability.txt").write_text(out)
            ctx.grab(out, re.compile(r"chains found|interfaces:|\[(PASS|FAIL)\]|held together|did not pass"))
            gmxmd.draw(ctx, d, "plot_assembly_stability.py")
        # a check of your own, written for your system (the worked example has one for its capped rod)
        if "stability" in c:
            st = c["stability"]
            own = lambda name: str(ctx.input(name)) if ctx.input(name).exists() else name      # beside the config, or in the package
            out = ctx.script(own(st["script"]), ["prod_so_far.gro", "nojump_prod.xtc", "prod.edr"], dict(st.get("settings", {})))
            (d / "stability.txt").write_text(out)
            ctx.grab(out)
            if "plot_script" in st:
                gmxmd.draw(ctx, d, own(st["plot_script"]))
        ctx.end()

    if ctx.begin("s9", "9 · The relaxed structure for the next stages" + (", and the hydrogen bonds asked for" if c.get("hbond") else "")):
        ch = pdb_chains(d / src.name)
        half = float(c.get("production_ns", 200)) * 500.0
        for hb in c.get("hbond", []):
            a, b = resindex_ranges(ch, hb["ref"]), resindex_ranges(ch, hb["target"])
            xvg = f"hb_{hb['name'].replace('-', '')}.xvg"
            ctx.sh(f"{ctx.gmx} hbond -s prod.tpr -f nojump_prod.xtc -r 'resindex {a[0]} to {a[1]}' "
                   f"-t 'resindex {b[0]} to {b[1]}' -num {xvg}", log="hbond_" + hb["name"])
            v = [float(l.split()[1]) for l in open(d / xvg) if l.strip() and l[0] not in "#@" and float(l.split()[0]) >= half]
            ctx.note(f"hydrogen bonds {hb['name']} ({hb['ref']} to {hb['target']}): {sum(v) / len(v):.1f} on average over the second half")
        ctx.sh(f"{ctx.gmx} trjconv -s prod.tpr -f prod.cpt -o fibril.pdb -pbc cluster", stdin="Protein\nProtein\n", log="fibril_pdb")
        br = chain_breaks(d / "fibril.pdb")
        cz = pdb_chains(d / "fibril.pdb")
        ctx.note(f"fibril.pdb: {len(pdb_atoms(d / 'fibril.pdb'))} atoms, {len(cz)} chains, {sum(len(r) for _, r in cz)} residues; " +
                 ("no chain breaks" if not br else f"{len(br)} chain break(s): do not use this file"))
        if br:
            raise StageError("fibril.pdb has chain breaks: " + ", ".join(f"chain {a} before residue {b} ({x:.1f} A)" for a, b, x in br[:5]))
        ctx.end()
