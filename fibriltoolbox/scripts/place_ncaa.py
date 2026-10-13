# PyMOL:  run place_ncaa.py
# Grafts the equilibrated non-canonical residue onto residues of the fibril (step 18).
# Takes any residue name, any number of sites and the side-chain torsions of any residue: superimpose on N, CA, C, CB, scan the first side-chain torsions for the
# fewest clashes, replace the residue. Each site avoids the ones placed before it.
import itertools
from pymol import cmd

FIBRIL = "fibril.pdb"          # equilibrated fibril, protein only
LINKER = "cyl.pdb"             # one frame of ACE-X-NME from the run in water
OUT    = "full_fib.pdb"
RESN   = "CYL"
SITES  = [["A", 239], ["L", 130]]          # [chain, residue number] as in FIBRIL
CHI    = [["N", "CA", "CB", "SG"], ["CA", "CB", "SG", "C4"]]    # torsions scanned (from ncaa_map.json)
CORE   = ["N", "CA", "CB", "SG", "C", "O"]  # atoms that sit where the old residue sat: not counted as clashes
MAX_FIT = 0.6                  # A: a larger fit RMSD means the site residue is broken or wrong
STEP   = 15                    # degrees

cmd.delete("all")                         # clean start if re-run in the same session
cmd.load(FIBRIL, "fib")
cmd.load(LINKER, "lnk")
cmd.remove(f"lnk and not resn {RESN}")    # strips ACE/NME caps and water


def turn(tag, angles):
    for names, a in zip(CHI, angles):
        cmd.set_dihedral(*[f"{tag} and name {n}" for n in names], a)


def place(tag, chain, resi, others):
    site = f"fib and chain {chain} and resi {resi}"
    if cmd.count_atoms(f"{site} and name CA") != 1:
        raise SystemExit(f"\n{tag}: {FIBRIL} has no residue {resi} in chain {chain}. Nothing was written.")
    cmd.create(tag, "lnk")
    pairs = []
    for a in ("N", "CA", "C", "CB"):      # atoms every non-Gly residue has
        pairs += [f"{tag} and name {a}", f"{site} and name {a}"]
    rmsd = cmd.pair_fit(*pairs)
    print(f"{tag} fit RMSD {rmsd:.3f}")
    if rmsd > MAX_FIT:
        raise SystemExit(f"\n{tag}: fit RMSD {rmsd:.1f} A is above {MAX_FIT} A, so chain {chain} residue {resi} "
                         f"of {FIBRIL} is not a normal residue.\nUsually it is split across the periodic box: "
                         "export the fibril again with whole molecules. Nothing was written.")

    L = f"{tag} and not hydro and not name {'+'.join(CORE)}"
    ref = f"((fib and not ({site}))" + "".join(" or " + o for o in others) + ") and not hydro"
    best, ang = None, None
    for angles in itertools.product(range(-180, 180, STEP), repeat=len(CHI)):
        turn(tag, angles)
        c = tuple(cmd.count_atoms(f"({L}) within {d} of ({ref})") for d in (1.5, 2.0, 2.5))
        if best is None or c < best:
            best, ang = c, angles
    turn(tag, ang)
    print(f"{tag} best {best} at " + " ".join(f"chi{i + 1}={a}" for i, a in enumerate(ang)))


tags = []
for k, (chain, resi) in enumerate(SITES, 1):
    tag = f"ncaa{k}_{chain}{resi}"
    place(tag, chain, resi, list(tags))
    tags.append(tag)

for tag, (chain, resi) in zip(tags, SITES):
    cmd.alter(tag, f"chain='{chain}'; resi='{resi}'; resn='{RESN}'; segi=''")
    cmd.remove(f"fib and chain {chain} and resi {resi}")
cmd.create("linked", "fib or " + " or ".join(tags))
cmd.sort("linked")

# save with TER only between chains (PyMOL's own TERs would also go before the new residue,
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
