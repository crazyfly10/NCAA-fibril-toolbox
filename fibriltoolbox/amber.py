"""AmberTools side: the tleap build with its two passes, acpype and its known repairs, restraint and index files."""
import re
import shutil
from pathlib import Path

from .core import StageError

LEAP_HEAD = """source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.opc
loadAmberPrep {name}.prepin
loadAmberParams frcmod.{name}
loadAmberParams frcmod.ff19SB_{name}
"""
LEAP_BAD = re.compile(r"FATAL|does not have a type|could not find", re.I)


def copy_ncaa_files(ctx, src, dst, name):
    """The three files that carry the residue into tleap."""
    for f in (f"{name}.prepin", f"frcmod.{name}", f"frcmod.ff19SB_{name}"):
        if not (Path(src) / f).exists():
            raise StageError(f"{f} not found in {src}: run ncaa-param first, or point [ncaa] files_from at it")
        if not (Path(dst) / f).exists():
            shutil.copy(Path(src) / f, Path(dst) / f)


def _box_cmd(unit, box, shape="Box"):
    if isinstance(box, (list, tuple)):
        return f"solvate{shape} {unit} OPCBOX {{ {' '.join(f'{b:.1f}' for b in box)} }}"
    return f"solvate{shape} {unit} OPCBOX {float(box):.1f}"


def tleap_build(ctx, d, name, out, body, unit, box, shape="Box", salt=0.15, pairs=None, extra_save=""):
    """Build, solvate and add salt. tleap takes ion counts, not a concentration, so it
    runs twice: once to read the charge and the number of waters, once with the ions. `body` is the tleap text that
    makes `unit`; returns the charge, the waters added, and the ion pairs used."""
    d = Path(d)
    head = LEAP_HEAD.format(name=name) + body + f"check {unit}\ncharge {unit}\n" + extra_save
    (d / f"{out}_pass1.leap").write_text(head + _box_cmd(unit, box, shape) + "\nquit\n")
    log1 = ctx.sh(f"tleap -f {out}_pass1.leap", cwd=d, log="tleap_pass1")
    (d / f"{out}_pass1.log").write_text(log1)
    m = re.search(r"Total unperturbed charge:\s*([-\d.]+)", log1)
    w = re.findall(r"Added (\d+) residues", log1)
    if LEAP_BAD.search(log1) or not m or not w:
        bad = [l for l in log1.split("\n") if LEAP_BAD.search(l)][:8]
        raise StageError("tleap (first pass) failed:\n" + "\n".join(bad or log1.strip().split("\n")[-15:]))
    charge, waters = float(m.group(1)), int(w[-1])
    q = round(charge)
    if pairs is None:
        pairs = round(salt * waters / 55.5)
    ions = (f"addIonsRand {unit} {'Cl-' if q > 0 else 'Na+'} 0\n" if q else "")
    ions += f"addIonsRand {unit} Na+ {pairs} Cl- {pairs}\n" if pairs else ""
    (d / f"{out}.leap").write_text(head + _box_cmd(unit, box, shape) + "\n" + ions +
                                   f"saveAmberParm {unit} {out}.prmtop {out}.inpcrd\nquit\n")
    log2 = ctx.sh(f"tleap -f {out}.leap", cwd=d, log="tleap_build")
    (d / f"{out}_leap.log").write_text(log2)
    if LEAP_BAD.search(log2) or not (d / f"{out}.prmtop").exists():
        bad = [l for l in log2.split("\n") if LEAP_BAD.search(l)][:8]
        raise StageError("tleap (build) failed:\n" + "\n".join(bad or log2.strip().split("\n")[-15:]))
    ctx.note(f"tleap: charge {charge:+.0f}, {waters} waters after solvation, {abs(q)} neutralising "
             f"{'Cl-' if q > 0 else 'Na+'} + {pairs} Na+/Cl- pairs ({salt} M); "
             f"{'Unit is OK' if 'Unit is OK' in log2 else 'check the leap log for warnings'}")
    return charge, waters, pairs


def acpype(ctx, d, out):
    """Amber -> GROMACS, then the ion-order repair acpype still needs."""
    d = Path(d)
    g = d / f"{out}.amb2gmx"
    if not (g / f"{out}_GMX.top").exists():
        ctx.sh(f"acpype -p {out}.prmtop -x {out}.inpcrd -b {out}", cwd=d, log="acpype")
    ctx.grab(ctx.script("fix_ion_order.py", [f"{out}_GMX.top", f"{out}_GMX.gro"], cwd=g), re.compile("ions|molecules"))
    return g


def fix_octahedron_box(gro):
    """acpype writes the truncated-octahedron box mirrored: negate box elements 6 and 9."""
    gro = Path(gro)
    shutil.copy(gro, str(gro) + ".orig")
    lines = gro.read_text().rstrip("\n").split("\n")
    v = [float(x) for x in lines[-1].split()]
    if len(v) != 9:
        raise StageError(f"{gro.name}: the box line has {len(v)} numbers, expected 9 for a truncated octahedron")
    v[5], v[8] = -v[5], -v[8]
    lines[-1] = "".join(f"{x:10.5f}" for x in v)
    gro.write_text("\n".join(lines) + "\n")


def first_molecule(top):
    """Name, atom count and residue count of the first moleculetype in a GROMACS topology (the solute)."""
    name, n, nres, m, sect = None, 0, 0, 0, None
    for line in open(top):
        s = line.strip()
        if s.startswith("["):
            sect = s.strip("[] ").lower()
            if sect == "moleculetype":
                m += 1
                if m == 2:
                    break
            continue
        if not s or s.startswith(";") or m != 1:
            continue
        f = s.split()
        if sect == "moleculetype" and name is None:
            name = f[0]
        elif sect == "atoms" and f[0].isdigit():
            n, nres = int(f[0]), int(f[2])
    return name, n, nres


def restraint_files(ctx, g, out, extra_define):
    """acpype's posre file also lists the ions, numbered past the end of the solute, and its topology includes
    no restraint file at all. Trim the first, write topol.top with both include blocks."""
    g = Path(g)
    name, natoms, nres = first_molecule(g / f"{out}_GMX.top")
    posre = g / f"posre_{out}.itp"
    shutil.copy(posre, str(posre) + ".acpype")
    kept = [l for l in posre.read_text().split("\n") if not re.match(r"\s*\d", l) or int(l.split()[0]) <= natoms]
    posre.write_text("\n".join(kept).rstrip("\n") + "\n")
    block = (f'#ifdef POSRES\n#include "posre_{out}.itp"\n#endif\n\n'
             f'#ifdef {extra_define}\n#include "posre_{extra_define[7:]}.itp"\n#endif\n\n')
    lines, m, res = (g / f"{out}_GMX.top").read_text().split("\n"), 0, []
    for l in lines:
        if re.match(r"\[ *moleculetype *\]", l):
            m += 1
            if m == 2:
                res.append(block.rstrip("\n") + "\n")
        res.append(l)
    (g / "topol.top").write_text("\n".join(res))
    n_posre = sum(1 for l in kept if re.match(r"\s*\d", l))
    ctx.note(f"solute molecule '{name}': {natoms} atoms, {nres} residues; posre_{out}.itp trimmed to {n_posre} heavy atoms; topol.top written")
    return name, natoms, nres


def make_index(ctx, g, gro, groups):
    """Index groups by residue index, with `gmx select`. groups: {name: selection}."""
    sel = "; ".join(f'"{k}" {v}' for k, v in groups.items())
    ctx.sh(f"{ctx.gmx} select -s {gro} -on index.ndx -select '{sel}'", cwd=g, log="select")
    counts, cur = {}, None
    for l in open(Path(g) / "index.ndx"):
        if l.startswith("["):
            cur = l.strip("[] \n")
            counts[cur] = 0
        elif cur:
            counts[cur] += len(l.split())
    ctx.note("index groups: " + ", ".join(f"{k} {v}" for k, v in counts.items()))
    return counts
