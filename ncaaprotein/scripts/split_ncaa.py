#!/usr/bin/env python3
"""
split_ncaa.py -- take the free NCAAs out of gmx insert-molecules' output and
write them as a PDB tleap reads correctly.

    python3 split_ncaa.py protein_box.pdb start.pdb ncaa_noh.pdb

insert-molecules appends the new molecules after the protein and leaves the
protein's atoms untouched, so everything past the protein's atom count is NCAA.
The protein itself is NOT taken from start.pdb: prep_protein_tleap.py writes it
on its own (protein_noh.pdb), so the protein's renaming never reaches the NCAAs:
renaming NME -> MET (GROMACS cuts NMET to NME) there would turn every NCAA's
N-methyl cap into a methionine.

What this writes:
  - heavy atoms only; tleap rebuilds the hydrogens from the residue's .prepin and the
    ACE/NME templates, as prep_protein_tleap.py does for the protein
  - residues renumbered 1..3N, atoms renumbered 1..
  - one TER after each NME and nowhere else: every ACE-X-NME is its own
    molecule, and no TER falls between ACE and X -- tleap reads every TER as a
    chain break, and X would lose its peptide bond to ACE
"""
import sys

SEQ = ("ACE", "NCA", "NME")                # the toolbox passes [ncaa] name for the middle one
HEAVY = {"ACE": 3, "NCA": 0, "NME": 2}     # heavy atoms per residue; the toolbox passes the count from ncaa-param


def atom_lines(path):
    return [l for l in open(path) if l.startswith(("ATOM", "HETATM"))]


def is_hydrogen(line):
    # Amber names every hydrogen H..., sometimes behind a digit (1HH3).
    # Decided from the name, not the element column, which GROMACS guesses.
    return line[12:16].strip().lstrip("0123456789").upper().startswith("H")


def residues(lines):
    """group consecutive atom lines into residues by (chain, resSeq, iCode, resName)"""
    out, key = [], None
    for l in lines:
        k = (l[21], l[22:26], l[26], l[17:21].strip())
        if k != key:
            out.append([])
            key = k
        out[-1].append(l)
    return out


if len(sys.argv) != 4:
    sys.exit(__doc__)
prot_pdb, ins_pdb, out_pdb = sys.argv[1:]

n_prot = len(atom_lines(prot_pdb))
allat = atom_lines(ins_pdb)
if len(allat) <= n_prot:
    sys.exit(f"{ins_pdb} has {len(allat)} atoms, the protein alone {n_prot}: nothing was inserted")

head = allat[n_prot - 1][17:21].strip()
if head in SEQ:
    sys.exit(f"atom {n_prot} of {ins_pdb} is in residue {head}, not the protein -- "
             "is protein_box.pdb the file insert-molecules was given with -f?")

res = residues(allat[n_prot:])
if len(res) % 3:
    sys.exit(f"{len(res)} NCAA residues is not a multiple of 3")
names = [r[0][17:21].strip() for r in res]
for i in range(0, len(names), 3):
    if tuple(names[i:i + 3]) != SEQ:
        sys.exit(f"NCAA {i // 3 + 1} reads {names[i:i + 3]}, expected {list(SEQ)}")

n_ncaa = len(res) // 3
out, serial = [], 0
for j, r in enumerate(res):
    rn = names[j]
    heavy = [l for l in r if not is_hydrogen(l)]
    if len(heavy) != HEAVY[rn]:
        sys.exit(f"NCAA {j // 3 + 1} {rn}: {len(heavy)} heavy atoms, expected {HEAVY[rn]}")
    for l in heavy:
        serial += 1
        l = l.rstrip("\n").ljust(54)
        out.append(f"{l[:6]}{serial % 100000:5d}{l[11:22]}{(j + 1) % 10000:4d}{l[26:]}\n")
    if rn == "NME":
        out.append("TER\n")
out.append("END\n")
open(out_pdb, "w").writelines(out)

h = sum(1 for l in allat[n_prot:] if is_hydrogen(l))
print(f"protein: {n_prot} atoms skipped (taken from protein_noh.pdb instead)")
print(f"NCAAs: {n_ncaa} x {'-'.join(SEQ)}, {serial} heavy atoms "
      f"({serial // n_ncaa} each), {h} hydrogens dropped for tleap to rebuild")
print(f"wrote {out_pdb}: {n_ncaa} TER records, one after each NME")
