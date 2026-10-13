"""
check_cap_fold.py -- does the AlphaFold cap actually have the HET-s solenoid fold?

Run in PyMOL:  run check_cap_fold.py

HET-s(218-289) is disordered as a free monomer: the beta-solenoid exists only
inside the fibril, templated by its neighbours. AlphaFold predicts the cap out
of that context, so nothing should be built on it until the predicted fold has
been checked against the real one. Superimpose the cap's first-segment core on
a fibril chain and take the RMSD; a large value means the prediction is not the
solenoid and the cap cannot stack, whatever superimpose_cap.py reports later.

Runs on all five models of an AlphaFold job, so it also ranks them.
"""
from pymol import cmd
import glob

FIBRIL    = "hets11.pdb"
REF_CHAIN = "K"                  # chains B-K are copies of the interior monomer
MODELS    = sorted(glob.glob("*_model_[0-9].cif"))  # AlphaFold Server output
PASS_RMSD = 1.5                  # angstroms

# the four core beta-strands in cap numbering; cap residue n = HET-s residue n + 216,
# so these are HET-s 226-234, 237-245, 262-270, 273-282
CAPC = list(range(10, 19)) + list(range(21, 30)) + list(range(46, 55)) + list(range(57, 67))
CAP_SEL = "resi " + "+".join(str(c) for c in CAPC)

cmd.delete("all")
cmd.load(FIBRIL, "fib")
n_ref = cmd.count_atoms(f"fib and chain {REF_CHAIN} and resi "
                        + "+".join(str(c + 216) for c in CAPC) + " and name CA")
print(f"reference: chain {REF_CHAIN} of {FIBRIL}, {n_ref} core CA "
      f"(expected {len(CAPC)})")
if n_ref != len(CAPC):
    print("  !! core selection does not match the fibril numbering -- check CAPC")

if not MODELS:
    raise SystemExit("no *_model_N.cif found -- set MODELS to your AlphaFold output")

print(f"\n{'model':<26}{'residues':>9}{'core CA':>9}{'RMSD/A':>9}  notes")
print("-" * 70)
results = []
for i, f in enumerate(MODELS):
    tag = "cap%02d" % i
    cmd.load(f, tag)
    n_res = cmd.count_atoms(f"{tag} and name CA")
    n_core = cmd.count_atoms(f"{tag} and {CAP_SEL} and name CA")

    # the construct should carry Val at the HET-s 239 position of its first
    # segment (cap 23) and the conjugation Cys at its C terminus (cap 130)
    notes = []
    for resi, want in ((23, "VAL"), (130, "CYS")):
        got = [a.resn for a in cmd.get_model(f"{tag} and resi {resi} and name CA").atom]
        if got and got[0] != want:
            notes.append(f"resi {resi} is {got[0]}, expected {want}")
        elif not got:
            notes.append(f"resi {resi} missing")

    if n_core != len(CAPC):
        print(f"{f:<26}{n_res:>9}{n_core:>9}{'--':>9}  core incomplete, skipped")
        continue

    pairs = []
    for c in CAPC:
        pairs += [f"{tag} and resi {c} and name CA",
                  f"fib and chain {REF_CHAIN} and resi {c + 216} and name CA"]
    r = cmd.pair_fit(*pairs)
    results.append((r, f, tag))
    if r >= PASS_RMSD:
        notes.append("NOT the solenoid fold")
    print(f"{f:<26}{n_res:>9}{n_core:>9}{r:>9.2f}  {'; '.join(notes)}")

print()
if not results:
    raise SystemExit("no model could be fitted -- check the residue numbering")

results.sort()
best_r, best_f, best_tag = results[0]
if best_r < PASS_RMSD:
    print(f"GATE passed: {best_f} fits the fibril core to {best_r:.2f} A "
          f"(< {PASS_RMSD} A), so the predicted fold is the solenoid.")
    print("superimpose_cap.py reads all the models and picks the cap itself; the")
    print("RMSD only says the fold is right, not that the tether is clear of the fibril.")
    print(f"All models are loaded here superimposed on chain {REF_CHAIN} -- look at them before going on.")
else:
    print(f"GATE FAILED: the best model fits only to {best_r:.2f} A. AlphaFold has not")
    print("reproduced the solenoid. Build the cap's rung region by threading its")
    print("sequence onto a fibril chain's backbone instead of predicting it.")
