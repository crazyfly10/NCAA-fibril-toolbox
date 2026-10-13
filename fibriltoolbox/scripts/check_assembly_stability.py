#!/usr/bin/env python3
"""
check_assembly_stability.py -- did a protein assembly hold together and stop changing shape during a run?

    python3 check_assembly_stability.py prod.gro nojump_prod.xtc [every_ps]

prod.gro / nojump_prod.xtc: the system with every chain whole and the chains kept together (gmx trjconv -pbc
cluster, then -pbc nojump). Works for any assembly: it needs no chain names, residue numbers or sequence.
Chains are found from the structure itself (a peptide bond longer than 2 A, or a jump in the atom order, starts
a new one).

Three things are followed, frame by frame:
  RMSD        C-alpha RMSD of the whole assembly from its first frame, after fitting
  Rg          radius of gyration of its C-alpha atoms
  interfaces  for every pair of chains in contact at the start, the fraction of their C-alpha contacts
              (pairs within CONTACT_CA in the first frame) still within 1.2 x that distance

and judged over the last EQ_FRACTION of the run:
  the RMSD and Rg have stopped drifting (the drift across the window is inside the noise, or not
  distinguishable from zero), and every interface keeps at least Q_MIN of its starting contacts.

Writes assembly_stability.csv and assembly_stability_meta.json; plot_assembly_stability.py draws
assembly_stability.png from them. Needs MDAnalysis.
"""
import json
import sys
import numpy as np

EQ_FRACTION = 0.5      # judge the last half of the run
CONTACT_CA = 8.0       # A: C-alpha pairs of two chains closer than this in the first frame are a contact
Q_MIN = 0.7            # an interface should keep at least this fraction of its contacts
MIN_CONTACTS = 5       # chain pairs with fewer contacts than this at the start are not an interface
PEPTIDE = 2.0          # A: a C-N distance above this between consecutive residues starts a new chain


def stats(t, y):
    """mean, residual std, slope and its standard error; neighbouring frames are correlated, so the number of
    points is reduced by the lag-1 autocorrelation of the residuals"""
    a, b = np.polyfit(t, y, 1)
    res = y - (a * t + b)
    r1 = np.corrcoef(res[:-1], res[1:])[0, 1] if len(y) > 3 and res.std() > 0 else 0.0
    r1 = min(max(r1, 0.0), 0.99)
    neff = max(3.0, len(y) * (1 - r1) / (1 + r1))
    se = res.std() / np.sqrt(neff) / t.std() if t.std() > 0 else 0.0
    return y.mean(), res.std(), a, se


def settled(t, y):
    mean, sd, a, se = stats(t, y)
    drift = a * (t[-1] - t[0])
    if sd == 0:
        return True, "constant"
    ok = abs(drift) < sd or abs(a) < 3 * se
    return ok, f"drift {drift:+.3g} over the window = {abs(drift) / sd:.2f} sigma, {abs(a) / se if se else 0:.1f} x its error"


def find_chains(u):
    """lists of residue indices, one per chain, from the geometry of the first frame"""
    prot = u.select_atoms("protein")
    res = [r for r in prot.residues if "CA" in r.atoms.names]
    chains, cur = [], [res[0]]
    for a, b in zip(res[:-1], res[1:]):
        c, n = a.atoms.select_atoms("name C"), b.atoms.select_atoms("name N")
        linked = len(c) == 1 and len(n) == 1 and np.linalg.norm(c.positions[0] - n.positions[0]) < PEPTIDE
        if linked:
            cur.append(b)
        else:
            chains.append(cur)
            cur = [b]
    chains.append(cur)
    return chains


def label(k):
    s = ""
    k += 1
    while k:
        k, r = divmod(k - 1, 26)
        s = chr(65 + r) + s
    return s


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    import MDAnalysis as mda
    from MDAnalysis.analysis.rms import rmsd

    every = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
    u = mda.Universe(sys.argv[1], sys.argv[2])
    u.trajectory[0]
    chains = find_chains(u)
    cas = [u.atoms[[r.atoms.select_atoms("name CA")[0].index for r in ch]] for ch in chains]
    allca = u.atoms[np.concatenate([c.indices for c in cas])]
    names = [label(k) for k in range(len(chains))]
    print(f"{len(chains)} chains found: " + ", ".join(f"{n} ({len(c)} residues)" for n, c in zip(names, chains)))

    ref = allca.positions.copy()
    pairs = []
    for i in range(len(cas)):
        for j in range(i + 1, len(cas)):
            d = np.linalg.norm(cas[i].positions[:, None, :] - cas[j].positions[None, :, :], axis=2)
            ii, jj = np.where(d < CONTACT_CA)
            if len(ii) >= MIN_CONTACTS:
                pairs.append((i, j, ii, jj))
    if len(chains) > 1 and not pairs:
        print(f"no pair of chains has {MIN_CONTACTS} or more C-alpha contacts within {CONTACT_CA:g} A: nothing to follow between chains")
    pnames = [f"{names[i]}-{names[j]}" for i, j, _, _ in pairs]
    if pairs:
        print(f"{len(pairs)} interfaces: " + ", ".join(f"{n} ({len(p[2])} contacts)" for n, p in zip(pnames, pairs)))

    T, R, G, Q, last = [], [], [], [], -1e30
    for ts in u.trajectory:
        if ts.time - last < every - 1e-6:
            continue
        last = ts.time
        x = allca.positions
        T.append(ts.time)
        R.append(rmsd(x, ref, center=True, superposition=True))
        G.append(float(np.sqrt(((x - x.mean(0)) ** 2).sum(1).mean())))
        row = []
        for i, j, ii, jj in pairs:
            d = np.linalg.norm(cas[i].positions[ii] - cas[j].positions[jj], axis=1)
            row.append(float((d < 1.2 * CONTACT_CA).mean()))
        Q.append(row)
    T, R, G = np.array(T), np.array(R), np.array(G)
    Q = np.array(Q).reshape(len(T), len(pairs))
    if len(T) < 10:
        sys.exit(f"only {len(T)} frames read: too few to judge")
    cut = int(len(T) * (1 - EQ_FRACTION))
    t = T[cut:]
    print(f"{len(T)} frames, {T[0] / 1000:.4g}-{T[-1] / 1000:.4g} ns; judged over the last {EQ_FRACTION:.0%}, from {T[cut] / 1000:.4g} ns\n")

    verdicts = []
    ok, why = settled(t, R[cut:])
    verdicts.append(("C-alpha RMSD has settled", ok, f"mean {R[cut:].mean():.2f} A; {why}"))
    ok, why = settled(t, G[cut:])
    verdicts.append(("radius of gyration has settled", ok, f"mean {G[cut:].mean():.2f} A; {why}"))
    for k, n in enumerate(pnames):
        q = Q[cut:, k].mean()
        verdicts.append((f"interface {n} keeps its contacts", q >= Q_MIN, f"{q:.0%} of {len(pairs[k][2])} kept; limit {Q_MIN:.0%}"))
    print("verdict")
    for name, ok, detail in verdicts:
        print(f"  [{'PASS' if ok else 'FAIL'}]  {name}  ({detail})")
    bad = sum(1 for _, ok, _ in verdicts if not ok)
    print("\n" + ("the assembly held together and its shape has settled" if bad == 0
                  else f"{bad} line(s) did not pass: look at assembly_stability.png before using the relaxed structure"))

    cols = ["time_ns", "rmsd_A", "rg_A"] + [f"contacts_{n}" for n in pnames]
    np.savetxt("assembly_stability.csv", np.column_stack([T / 1000.0, R, G] + [Q[:, k] for k in range(len(pairs))]),
               delimiter=",", header=",".join(cols), comments="", fmt="%.4f")
    json.dump(dict(chains=names, sizes=[len(c) for c in chains], interfaces=pnames, contacts=[int(len(p[2])) for p in pairs],
                   EQ_FRACTION=EQ_FRACTION, Q_MIN=Q_MIN, CONTACT_CA=CONTACT_CA,
                   verdicts=[[n, bool(o), d] for n, o, d in verdicts]), open("assembly_stability_meta.json", "w"), indent=1)
    print(f"wrote assembly_stability.csv ({len(T)} frames, {len(cols)} columns) and assembly_stability_meta.json  "
          "(plot them with plot_assembly_stability.py)")


if __name__ == "__main__":
    main()
