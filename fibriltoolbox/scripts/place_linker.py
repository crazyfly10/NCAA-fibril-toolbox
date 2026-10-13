# PyMOL:  run place_linker.py
# Grafts the equilibrated CYL (from cyl.pdb) onto two residues of the fibril.
from pymol import cmd

FIBRIL = "fibril.pdb"          # equilibrated fibril, protein only, chains A-L
LINKER = "cyl.pdb"             # ACE-CYL-NME frame from §17
OUT    = "full_fib.pdb"
SITE_A, SITE_L = 239, 130      # residues replaced by CYL (numbering as in FIBRIL)
MAX_FIT = 0.6                  # A: a larger fit RMSD means the site residue is broken or wrong

cmd.delete("all")                         # clean start if re-run in the same session
cmd.load(FIBRIL, "fib")
cmd.load(LINKER, "cyl")
cmd.remove("cyl and not resn CYL")        # strips ACE/NME caps and water

def place(tag, chain, resi, others):
    site = f"fib and chain {chain} and resi {resi}"
    cmd.create(tag, "cyl")
    pairs = []
    for a in ("N", "CA", "C", "CB"):      # atoms every non-Gly residue has
        pairs += [f"{tag} and name {a}", f"{site} and name {a}"]
    rmsd = cmd.pair_fit(*pairs)
    print(f"{tag} fit RMSD {rmsd:.3f}")
    if rmsd > MAX_FIT:
        raise SystemExit(f"\n{tag}: fit RMSD {rmsd:.1f} A is above {MAX_FIT} A, so chain {chain} residue {resi} "
                         f"of {FIBRIL} is not a normal residue.\nUsually it is split across the periodic box: "
                         "export fibril.pdb again with whole molecules (section 9). Nothing was written.")

    L = f"{tag} and not hydro and not name N+CA+CB+SG+C+O"
    ref = f"((fib and not ({site}))" + "".join(" or " + o for o in others) + ") and not hydro"
    best, ang = None, None
    for a1 in range(-180, 180, 15):                      # chi1  CB-SG
        cmd.set_dihedral(f"{tag} and name N",  f"{tag} and name CA",
                         f"{tag} and name CB", f"{tag} and name SG", a1)
        for a2 in range(-180, 180, 15):                  # chi2  SG-C4
            cmd.set_dihedral(f"{tag} and name CA", f"{tag} and name CB",
                             f"{tag} and name SG", f"{tag} and name C4", a2)
            c = tuple(cmd.count_atoms(f"({L}) within {d} of ({ref})")
                      for d in (1.5, 2.0, 2.5))
            if best is None or c < best: best, ang = c, (a1, a2)

    cmd.set_dihedral(f"{tag} and name N",  f"{tag} and name CA",
                     f"{tag} and name CB", f"{tag} and name SG", ang[0])
    cmd.set_dihedral(f"{tag} and name CA", f"{tag} and name CB",
                     f"{tag} and name SG", f"{tag} and name C4", ang[1])
    print(f"{tag} best {best} at chi1={ang[0]} chi2={ang[1]}")

place("cylA", "A", SITE_A, [])
place("cylL", "L", SITE_L, ["cylA"])      # second site avoids the first

cmd.alter("cylA", f"chain='A'; resi='{SITE_A}'; resn='CYL'; segi=''")
cmd.alter("cylL", f"chain='L'; resi='{SITE_L}'; resn='CYL'; segi=''")
cmd.remove(f"fib and chain A and resi {SITE_A}")
cmd.remove(f"fib and chain L and resi {SITE_L}")
cmd.create("linked", "fib or cylA or cylL")
cmd.sort("linked")

# save with TER only between chains (PyMOL's own TERs would also go before CYL,
# and tleap reads every TER as a chain break)
cmd.set("pdb_use_ter_records", 0)
cmd.save(OUT, "linked")
atoms = [l for l in open(OUT) if l.startswith(("ATOM", "HETATM"))]
with open(OUT, "w") as f:
    for i, l in enumerate(atoms):
        if i and l[21] != atoms[i - 1][21]:
            f.write("TER\n")
        f.write(l)
    f.write("TER\nEND\n")
print(f"wrote {OUT}")
