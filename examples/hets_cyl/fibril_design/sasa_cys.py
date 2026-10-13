# PyMOL:  run sasa_cys.py      (on the mutated m1 of step 2)
# Side-chain solvent-accessible area of the listed residues in every chain of m1:
# chain, residue, side-chain atom count, area in A^2.
from pymol import cmd
SITES = (239,)                     # residues to report
cmd.set("dot_solvent", 1)
cmd.set("dot_density", 3)
for ch in cmd.get_chains("m1"):
    for r in SITES:
        sel = f"m1 and chain {ch} and resi {r} and not name N+CA+C+O"
        print(ch, r, cmd.count_atoms(sel), round(cmd.get_area(sel), 2))
