#!/usr/bin/env python3
"""
check_box2.py -- check a box before solvation, for ANY box shape
(rectangular, triclinic, dodecahedron, octahedron). Reads only the .gro.

    python3 check_box2.py box.gro [min_image_nm]

Prints three checks:
  1. every atom inside the cell -- what check_box.py tested. This only means
     something for a RECTANGULAR box. A dodecahedron is stored as a slanted
     block spanned by three vectors; editconf centres the protein in the
     12-faced cell, and an elongated protein can poke past the block's slanted
     faces. Those atoms just sit in the neighbouring cell, which periodic
     boundaries handle, so for a slanted box this is reported, not failed.
  2. the closest distance between the protein and its own periodic images,
     which is what actually decides whether the box is big enough. It must be
     at least min_image_nm (default 2.0 nm, the floor the guide uses);
     editconf -d 1.2 gives about 2.4 nm.
  3. that the protein's centre lies within 0.1 nm of the box centre.

All atoms are used, hydrogens included, as gmx mindist -pi does.
"""
import sys
import numpy as np

MIN_IMAGE = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0     # nm


def read_gro(path):
    L = open(path).readlines()
    n = int(L[1])
    xyz = np.array([[float(l[20 + 8 * i:28 + 8 * i]) for i in range(3)] for l in L[2:2 + n]])
    b = [float(v) for v in L[2 + n].split()]
    if len(b) == 3:
        b += [0.0] * 6
    # .gro order: v1x v2y v3z v1y v1z v2x v2z v3x v3y
    box = np.array([[b[0], b[3], b[4]],
                    [b[5], b[1], b[6]],
                    [b[7], b[8], b[2]]])        # rows are the box vectors
    return xyz, box


def min_dist(A, B, chunk=400):
    """closest distance between any point of A and any point of B"""
    try:
        from scipy.spatial import cKDTree
        return float(cKDTree(B).query(A, k=1)[0].min())
    except ImportError:
        best = np.inf
        B32 = B.astype(np.float32)
        for i in range(0, len(A), chunk):
            d = A[i:i + chunk, None, :].astype(np.float32) - B32[None, :, :]
            best = min(best, float(np.sqrt((d * d).sum(-1)).min()))
        return best


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    xyz, box = read_gro(sys.argv[1])
    tric = np.count_nonzero(box - np.diag(np.diag(box))) > 0
    print(f"{len(xyz)} atoms; box vectors (nm):")
    for name, v in zip("abc", box):
        print(f"   {name} = ({v[0]:7.3f} {v[1]:7.3f} {v[2]:7.3f})")
    print(f"   shape: {'triclinic (e.g. dodecahedron)' if tric else 'rectangular'}, "
          f"volume {abs(np.linalg.det(box)):.1f} nm^3")

    # 1. inside the cell: fractional coordinates must all lie in [0, 1)
    frac = xyz @ np.linalg.inv(box)
    out = int(((frac < 0) | (frac >= 1)).any(1).sum())
    centre = box.sum(0) / 2
    off = np.linalg.norm(xyz.mean(0) - centre)
    note = ("" if not (tric and out) else
            "  -- fine for a slanted box: they sit in the neighbouring cell")
    print(f"\n1. atoms outside the stored cell: {out}   "
          f"(protein centre {off:.2f} nm from the box centre){note}")

    # 2. protein to its own periodic images, over all 26 neighbouring cells
    best, shift = np.inf, None
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                if i == j == k == 0:
                    continue
                t = i * box[0] + j * box[1] + k * box[2]
                d = min_dist(xyz, xyz + t)
                if d < best:
                    best, shift = d, (i, j, k)
    span = max(np.ptp(xyz @ u) for u in np.eye(3))
    print(f"2. closest protein-to-image distance: {best:.2f} nm  (image {shift}); "
          f"protein spans up to {span:.2f} nm along x, y or z")

    inside_ok = out == 0 or tric          # only a rectangular box must hold every atom
    centred = off < 0.1
    ok = inside_ok and best >= MIN_IMAGE and centred
    print(f"\n{'PASS' if ok else 'FAIL'}: image distance {best:.2f} nm "
          f"{'>=' if best >= MIN_IMAGE else '<'} {MIN_IMAGE:g} nm"
          + ("" if centred else f", protein {off:.2f} nm off centre")
          + ("" if inside_ok else f", {out} atoms outside a rectangular box"))
    if not inside_ok or not centred:
        print("   re-run editconf with -c to centre the protein")
    if best < MIN_IMAGE:
        print("   too close to its own image: raise -d in editconf")


if __name__ == "__main__":
    main()
