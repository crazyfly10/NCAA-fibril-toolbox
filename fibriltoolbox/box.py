"""Box size worked out from the structure: `box = "auto"` (see the README, "Your own system").

The solute may come within reach of its own periodic image when a flexible part stretches out. For every flexible
part a sphere is drawn around its anchor atom, with the part's reach as radius; along each axis the box has to
hold the solute and those spheres, plus a margin, and it shrinks a few percent in NPT. What is left on each side
is the buffer handed to tleap's solvateBox.

Every stage that makes a box takes "auto" (the default) or a number: `auto_buffer` for the fibril with the residue
attached, `auto_gap_buffer` for the fibril among free molecules, `auto_lengths` for the protein-only fibril.
"""
import numpy as np

from .core import pdb_atoms


def auto_buffer(pdb, parts, margin_nm=2.0, shrink=0.95, floor_A=10.0):
    """parts: [(chain, resid, reach in nm), ...]. Returns ([bx, by, bz] in A, text lines, solute extent in A, extent with the parts stretched out in A)."""
    atoms = pdb_atoms(pdb)
    xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in atoms])
    lo, hi = xyz.min(0), xyz.max(0)
    elo, ehi = lo.copy(), hi.copy()
    lines = []
    for chain, resid, reach in parts:
        a = [x for l, x in zip(atoms, xyz) if l[21] == chain and int(l[22:26]) == int(resid) and l[12:16].strip() == "CA"]
        if not a:
            raise ValueError(f"no CA for chain {chain} residue {resid} in {pdb}")
        elo, ehi = np.minimum(elo, a[0] - reach * 10), np.maximum(ehi, a[0] + reach * 10)
        lines.append(f"flexible part at {chain}{resid}: reach {reach:.2f} nm")
    extent = hi - lo
    need = (ehi - elo) + margin_nm * 10.0                 # box length wanted after NPT, A
    buf = np.maximum((need / shrink - extent) / 2.0, floor_A)
    box = (extent + 2 * buf) / 10.0
    lines.append("solute extent %.1f x %.1f x %.1f nm; with the flexible parts stretched out %.1f x %.1f x %.1f nm"
                 % (*(extent / 10), *((ehi - elo) / 10)))
    lines.append("buffer %.1f / %.1f / %.1f A  ->  box about %.1f x %.1f x %.1f nm before NPT (%.0f nm3, roughly %s atoms)"
                 % (*buf, *box, box.prod(), format(round(box.prod() * 0.135) * 1000, ',')))
    return [round(float(b), 1) for b in buf], lines, extent, ehi - elo


def auto_gap_buffer(core_pdb, other_pdbs, gap_nm, shrink=0.95, floor_A=10.0):
    """A core (the fibril) with free molecules around it: the free molecules cross the box faces anyway, so what the
    box has to prevent is one of them touching the core and the core's image at once. Along each axis the core and
    its image are kept `gap_nm` apart after NPT. Returns ([bx, by, bz] in A around the whole solute, text lines, solute extent in A, core extent in A)."""
    def xyz(f):
        return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in pdb_atoms(f)])
    core = xyz(core_pdb)
    every = np.vstack([core] + [xyz(f) for f in other_pdbs])
    cext, ext = core.max(0) - core.min(0), every.max(0) - every.min(0)
    need = cext + gap_nm * 10.0
    buf = np.maximum((need / shrink - ext) / 2.0, floor_A)
    box = (ext + 2 * buf) / 10.0
    lines = ["fibril extent %.1f x %.1f x %.1f nm; with the free molecules %.1f x %.1f x %.1f nm" % (*(cext / 10), *(ext / 10)),
             "fibril kept %.2f nm from its own image after NPT" % gap_nm,
             "buffer %.1f / %.1f / %.1f A  ->  box about %.1f x %.1f x %.1f nm before NPT (%.0f nm3, roughly %s atoms)"
             % (*buf, *box, box.prod(), format(round(box.prod() * 0.135) * 1000, ','))]
    return [round(float(b), 1) for b in buf], lines, ext, cext


def auto_lengths(xyz_nm, parts, margin_nm=3.0, shrink=0.95):
    """Box lengths (nm) for a structure already aligned to the axes. parts: [(anchor xyz in nm, reach in nm, label), ...]
    for the segments known to be flexible. Returns ([lx, ly, lz], text lines, extent in nm)."""
    xyz_nm = np.asarray(xyz_nm)
    lo, hi = xyz_nm.min(0), xyz_nm.max(0)
    elo, ehi = lo.copy(), hi.copy()
    lines = []
    for a, reach, label in parts:
        elo, ehi = np.minimum(elo, np.asarray(a) - reach), np.maximum(ehi, np.asarray(a) + reach)
        lines.append(f"flexible part at {label}: reach {reach:.2f} nm")
    box = ((ehi - elo) + margin_nm) / shrink
    lines.append("extent %.1f x %.1f x %.1f nm" % tuple(hi - lo) + ("; with the flexible parts stretched out %.1f x %.1f x %.1f nm"
                 % tuple(ehi - elo) if parts else "; no flexible parts listed, so the margin has to cover any tail that moves"))
    lines.append("margin %.1f nm  ->  box %.2f x %.2f x %.2f nm before NPT (%.0f nm3, roughly %s atoms)"
                 % (margin_nm, *box, box.prod(), format(round(box.prod() * 0.13) * 1000, ',')))
    return [round(float(b), 3) for b in box], lines, hi - lo, ehi - elo


def atoms_in(lengths_nm, per_nm3=135):
    """Rough atom count of a water box (four-point water): about 135 per nm3 of the box asked of tleap, 130 with gmx solvate."""
    return float(np.prod(lengths_nm)) * per_nm3


def advice(setting, chosen, n_auto, n_fixed, fixed_text, flexible_text, check):
    """What the choice of box costs and risks, said in every stage that makes a box.
    chosen: "auto" or the user's number; n_auto, n_fixed: atom estimates for the automatic box and for the fixed
    one (the user's number, or a 12 A buffer for comparison when the box is automatic)."""
    r = n_auto / n_fixed
    k = lambda n: f"{round(n, -3):,.0f}"
    if chosen == "auto" and r <= 1.05:
        return [f"automatic box: roughly {k(n_auto)} atoms, no larger than with {fixed_text} ({k(n_fixed)}), so here it costs no extra computer time",
                f"to choose the size yourself, set {setting} to a number; {check} says whether the box held"]
    if chosen == "auto":
        return [f"automatic box: this probably removes periodic-image artefacts, but it costs computer time: roughly {k(n_auto)} atoms "
                f"against {k(n_fixed)} with {fixed_text} ({r:.1f} times), and each nanosecond takes about that much longer",
                f"to choose the size yourself, set {setting} to a number; {check} says whether the box held"]
    out = [f"box set by you ({fixed_text}): roughly {k(n_fixed)} atoms; the automatic box would be roughly {k(n_auto)} atoms ({r:.1f} times)"]
    if r > 1.1:
        out.append(f"WARNING: this box is smaller than the automatic one. {flexible_text} A small margin such as 12 A lets a "
                   f"flexible region near the box edge meet its own periodic image: read {check} before using the run")
    return out


def built(box_nm, extent_nm, what="the construct", stretched_nm=None):
    """The box as it was built, and how far the construct is from it: text lines for the report.
    box_nm, extent_nm: lengths along x, y, z; stretched_nm: the extent with the flexible parts stretched out, if known."""
    box, ext = np.asarray(box_nm, float), np.asarray(extent_nm, float)
    side = (box - ext) / 2.0
    lines = ["box as built: %.2f x %.2f x %.2f nm (x, y, z); %s is %.1f x %.1f x %.1f nm" % (*box, what, *ext),
             "water between %s and the box edge: %.2f / %.2f / %.2f nm on each side; from %s to its own periodic image: "
             "%.2f / %.2f / %.2f nm" % (what, *side, what, *(box - ext))]
    if stretched_nm is not None:
        lines.append("with the flexible parts stretched out as far as they reach, the distance to the image would still be "
                     "%.2f / %.2f / %.2f nm" % tuple(box - np.asarray(stretched_nm, float)))
    lines.append("these are distances before NPT, which shrinks the box by a few percent")
    return lines
