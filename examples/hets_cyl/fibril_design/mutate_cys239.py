# PyMOL:  run mutate_cys239.py      (after the load / split_states / remove commands of §1)
# Mutates residue 239 to Cys in every chain of m1 with the Mutagenesis wizard and
# saves hets_3mer_mut.pdb. For each chain it takes the rotamer with the lowest strain
# (the wizard's bump score), which is the least-clashing one, and prints what it chose.
from pymol import cmd

OBJ, RESI, NEW = "m1", 239, "CYS"
OUT = "hets_3mer_mut.pdb"

cmd.wizard("mutagenesis")
cmd.do("refresh_wizard")
w = cmd.get_wizard()
w.set_mode(NEW)
for ch in cmd.get_chains(OBJ):
    w.do_select(f"/{OBJ}//{ch}/{RESI}/")
    scores = list(w.bump_scores)                 # one strain value per rotamer offered
    best = scores.index(min(scores)) + 1 if scores else 1
    cmd.frame(best)                              # rotamers are shown as states 1..N
    w.apply()
    print(f"chain {ch}: {len(scores)} rotamers, strain {[round(s, 1) for s in scores]}, took rotamer {best}")
cmd.set_wizard()

for ch in cmd.get_chains(OBJ):
    names = {a.resn for a in cmd.get_model(f"{OBJ} and chain {ch} and resi {RESI}").atom}
    if names != {NEW}:
        raise SystemExit(f"chain {ch} residue {RESI} is {sorted(names)}, not {NEW}: nothing saved")
cmd.save(OUT, OBJ)
print(f"wrote {OUT}: {NEW}{RESI} in chains {''.join(cmd.get_chains(OBJ))}")
