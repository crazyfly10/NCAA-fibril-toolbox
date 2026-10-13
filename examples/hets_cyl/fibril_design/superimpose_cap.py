"""
superimpose_cap.py -- choose the AlphaFold cap model and its axial shift, and dock
it on the far end of the fibril as chain L.

Run in PyMOL, in the folder with hets11.pdb and the AlphaFold Server models:
    run superimpose_cap.py

1. Ghost: the screw transform between interior chains F and G, applied to chain K,
   gives a ghost monomer one step past the end of the fibril.
2. Dock: every model's 37 core CA (cap 10-18, 21-29, 46-54, 57-66 = HET-s 226-234,
   237-245, 262-270, 273-282) are superimposed on the ghost.
3. Screen: a model fails if its disordered tether segment (resi 75-101) threads
   through the fibril -- atoms below chain K's outer face and inside the fibril's
   cross-section that are still there after a 3.0 A shift, so no allowed shift
   can pull it back out.
4. Rank: the remaining models by the number of clashes at the dock.
5. Scan: the best model is stepped 0-3.0 A along the fibril axis. Past 3.0 A the
   rungs move out of hydrogen-bonding reach of chain K. At each shift, every
   heavy-atom contact < 2.5 A is counted, and a pair closer than HARD (1.5 A) is a
   clash; contacts of 1.5-2.5 A are tolerated, since minimisation and MD relax them.
   If clashes remain at 3.0 A (the N-terminal tail, resi 1-9, is the usual cause),
   resi 1-7 are deleted and the scan is repeated. If the model still has no usable
   shift, the next model in the ranking is tried.
6. Decide: the smallest shift with no clash while the core is still within 3.5 A
   of chain K. The cap is moved there, made chain L and saved with the fibril.

Everything it decides is printed: the dock screen, each scan, the decision, and
every contact < 2.5 A left at the chosen shift, so the result can be judged by eye.
Nothing is saved if no model passes.
"""
import glob
import numpy as np
from pymol import cmd

FIBRIL   = "hets11.pdb"
MODELS   = sorted(glob.glob("*_model_[0-9].cif"))    # AlphaFold Server output
OUT      = "hets11_capped.pdb"

HETS_CORE = list(range(226, 235)) + list(range(237, 246)) + list(range(262, 271)) + list(range(273, 283))
CAP_CORE  = [r - 216 for r in HETS_CORE]          # cap residue n = HET-s residue n + 216
TAIL      = (1, 9)        # N-terminal tail, outside the core
TAIL_DEL  = (1, 7)        # what is deleted if the tail still clashes at 3.0 A
DISORD    = (75, 101)     # disordered tether segment, checked for threading
SHIFTS    = np.round(np.arange(0.0, 3.01, 0.1), 2)
CONTACT_CUT = 2.5         # A: contacts listed
HARD      = 1.5           # A: closer than this is a clash
CORE_REACH = 3.5          # A: the core must stay this close to chain K
THREAD_DEPTH = 1.0        # A below chain K's outer face that counts as inside


def atoms(sel):
    m = cmd.get_model(sel)
    xyz = np.array([a.coord for a in m.atom], float)
    lab = [(a.chain, int(a.resi), a.resn, a.name) for a in m.atom]
    return xyz, lab


def ca(obj, chain, resis):
    sel = f"{obj} and name CA and resi " + "+".join(map(str, resis))
    if chain:
        sel += f" and chain {chain}"
    d = {int(a.resi): a.coord for a in cmd.get_model(sel).atom}
    return np.array([d[r] for r in resis], float)


def kabsch(P, Q):
    """R, t with R @ P_i + t ~= Q_i"""
    pc, qc = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    D = np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))])
    R = Vt.T @ D @ U.T
    return R, qc - R @ pc


def hull2d(p):
    """convex hull of 2-D points, counter-clockwise (monotone chain)"""
    p = sorted(map(tuple, np.round(p, 4)))
    cross = lambda o, a, b: (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lo, hi = [], []
    for q in p:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], q) <= 0:
            lo.pop()
        lo.append(q)
    for q in reversed(p):
        while len(hi) >= 2 and cross(hi[-2], hi[-1], q) <= 0:
            hi.pop()
        hi.append(q)
    return np.array(lo[:-1] + hi[:-1])


def in_hull(h, q):
    """True for each point of q inside the counter-clockwise hull h"""
    a, b = h, np.roll(h, -1, axis=0)
    c = (b[:, 0] - a[:, 0]) * (q[:, None, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (q[:, None, 0] - a[:, 0])
    return (c >= 0).all(axis=1)


# ------------------------------------------------------------------ fibril and ghost
cmd.delete("all")
cmd.load(FIBRIL, "fib")
cmd.remove("fib and hydro")
F, G, K = (ca("fib", c, HETS_CORE) for c in "FGK")
R, t = kabsch(F, G)                               # one screw step, F -> G
ghost = (R @ K.T).T + t                           # core CA of a 12th monomer
w, v = np.linalg.eig(R)
axis = np.real(v[:, np.argmin(abs(w - 1))])
axis /= np.linalg.norm(axis)
if axis @ (G.mean(0) - F.mean(0)) < 0:
    axis = -axis                                  # points out of the fibril, past K
rise = (G.mean(0) - F.mean(0)) @ axis

fxyz, flab = atoms("fib")
kmask = np.array([l[0] == "K" for l in flab])
o = K.mean(0)
k_top = ((fxyz[kmask] - o) @ axis).max()          # chain K's outer face
e1 = np.cross(axis, [1.0, 0, 0])
if np.linalg.norm(e1) < 0.1:
    e1 = np.cross(axis, [0, 1.0, 0])
e1 /= np.linalg.norm(e1)
e2 = np.cross(axis, e1)
footprint = hull2d(np.c_[(fxyz - o) @ e1, (fxyz - o) @ e2])   # fibril cross-section

print(f"fibril {FIBRIL}: {len(fxyz)} heavy atoms, rise {rise:.2f} A along the axis")
if not MODELS:
    raise SystemExit("no *_model_N.cif found -- put the AlphaFold Server models next to this script")
print(f"models: {', '.join(MODELS)}")


def contacts(x, sel):
    """(distance, cap atom, fibril atom) for every heavy-atom pair < CONTACT_CUT"""
    idx = np.where(sel)[0]
    out = []
    for a in range(0, len(idx), 200):                 # chunks keep memory small
        i = idx[a:a + 200]
        D = np.linalg.norm(x[i, None, :] - fxyz[None, :, :], axis=-1)
        for p, q in zip(*np.where(D < CONTACT_CUT)):
            out.append((D[p, q], i[p], q))
    return out


def inside(x, sel):
    """atoms of x[sel] below chain K's outer face and within the fibril's cross-section"""
    p = x[sel] - o
    below = (p @ axis) < k_top - THREAD_DEPTH
    return int((below & in_hull(footprint, np.c_[p @ e1, p @ e2])).sum())


def scan(x0, lab, keep):
    resi = np.array([l[1] for l in lab])
    core = keep & np.isin(resi, CAP_CORE)
    rows = []
    for s in SHIFTS:
        x = x0 + s * axis
        c = contacts(x, keep)
        clash = [q for q in c if q[0] < HARD]
        dk = np.linalg.norm(x[core][:, None, :] - fxyz[kmask][None, :, :], axis=-1).min(1)
        rows.append(dict(shift=s, n=len(c), clash=len(clash),
                         tail=sum(1 for q in clash if TAIL[0] <= resi[q[1]] <= TAIL[1]),
                         dmin=min([q[0] for q in c], default=float("nan")),
                         core_k=dk.min(), core_n=int((dk < CORE_REACH).sum()), contacts=c))
    return rows


def show(rows, title):
    print(f"\n  {title}")
    print(f"   shift | <{CONTACT_CUT} A  <{HARD} A (tail) | closest | core->K  core atoms <{CORE_REACH} A")
    print("  " + "-" * 72)
    for r in rows:
        print(f"   {r['shift']:5.1f} | {r['n']:6d}  {r['clash']:6d} ({r['tail']:4d}) | {r['dmin']:7.2f} |"
              f" {r['core_k']:7.2f}  {r['core_n']:8d}")


# ------------------------------------------------------------------ dock and screen
results = []
for f in MODELS:
    name = f.rsplit("_model_", 1)[1].replace(".cif", "")
    tag = "model_" + name
    cmd.load(f, tag)
    cmd.remove(f"{tag} and hydro")
    Rc, tc = kabsch(ca(tag, None, CAP_CORE), ghost)
    x, lab = atoms(tag)
    x = (Rc @ x.T).T + tc
    rmsd = np.sqrt((((Rc @ ca(tag, None, CAP_CORE).T).T + tc - ghost) ** 2).sum(1).mean())
    resi = np.array([l[1] for l in lab])
    dis = (resi >= DISORD[0]) & (resi <= DISORD[1])
    c0 = contacts(x, np.ones(len(x), bool))
    in0, in3 = inside(x, dis), inside(x + SHIFTS[-1] * axis, dis)
    fail = []
    if in3:
        fail.append(f"resi {DISORD[0]}-{DISORD[1]} threads through the fibril "
                    f"({in0} atoms inside at the dock, still {in3} at {SHIFTS[-1]:.1f} A)")
    results.append(dict(tag=tag, file=f, rmsd=rmsd, x=x, lab=lab, n0=len(c0),
                        clash0=sum(1 for q in c0 if q[0] < HARD), in0=in0, in3=in3, fail=fail))

print(f"\nDOCK SCREEN (each model at the ideal stacking position, shift 0)")
print(f" model     core RMSD | <{CONTACT_CUT} A  <{HARD} A | resi {DISORD[0]}-{DISORD[1]} inside (0 / {SHIFTS[-1]:.1f} A) | verdict")
print("-" * 82)
for r in results:
    print(f" {r['tag']:9s} {r['rmsd']:8.2f}  | {r['n0']:6d}  {r['clash0']:6d} |"
          f" {r['in0']:12d} / {r['in3']:<8d}  | {'FAIL' if r['fail'] else 'ok'}")
    for why in r["fail"]:
        print(f"            -> FAIL: {why}")

ranked = sorted([r for r in results if not r["fail"]], key=lambda r: (r["clash0"], r["n0"], r["rmsd"]))
if not ranked:
    raise SystemExit("\nno usable cap: every model threads through the fibril")
print("\nRANKING (fewest clashes at the dock): " + "  ".join(f"{r['tag']}({r['clash0']})" for r in ranked))

# ------------------------------------------------------------------ shift scan
chosen = None
for r in ranked:
    print(f"\n=== {r['tag']} ===")
    resi = np.array([l[1] for l in r["lab"]])
    keep = np.ones(len(resi), bool)
    rows = scan(r["x"], r["lab"], keep)
    show(rows, "shift scan, whole cap")
    ok = [q for q in rows if q["clash"] == 0 and q["core_k"] <= CORE_REACH]
    trimmed = False
    if not ok:
        last = rows[-1]
        print(f"\n  {last['clash']} clashes left at {SHIFTS[-1]:.1f} A, {last['tail']} of them in the "
              f"N-terminal tail: deleting resi {TAIL_DEL[0]}-{TAIL_DEL[1]} and scanning again")
        keep = ~((resi >= TAIL_DEL[0]) & (resi <= TAIL_DEL[1]))
        rows = scan(r["x"], r["lab"], keep)
        show(rows, f"shift scan, resi {TAIL_DEL[0]}-{TAIL_DEL[1]} deleted")
        ok = [q for q in rows if q["clash"] == 0 and q["core_k"] <= CORE_REACH]
        trimmed = True
    if ok:
        chosen = (r, ok[0], trimmed)
        break
    print(f"  -> no usable shift within {SHIFTS[-1]:.1f} A; trying the next model")

if chosen is None:
    raise SystemExit(f"\nno usable cap: no model clears {HARD} A within {SHIFTS[-1]:.1f} A -- nothing saved")

# ------------------------------------------------------------------ decision and save
r, best, trimmed = chosen
lab = r["lab"]
print("\nDECISION")
print(f"  cap model    {r['tag']} ({r['file']}): core RMSD {r['rmsd']:.2f} A, "
      f"{r['clash0']} clashes at the dock, the fewest of the models that passed")
print(f"  N-tail       {'resi %d-%d deleted' % TAIL_DEL if trimmed else 'kept'}")
print(f"  axial shift  {best['shift']:.1f} A: the smallest shift with no pair closer than {HARD} A "
      f"while the core stays within {CORE_REACH} A of chain K "
      f"(closest core-K {best['core_k']:.2f} A, {best['core_n']} core atoms within {CORE_REACH} A)")
print(f"  the cap core sits {rise + best['shift']:.2f} A past chain K (rise {rise:.2f} A)")
print(f"\n  contacts < {CONTACT_CUT} A at {best['shift']:.1f} A (all tolerated, none < {HARD} A):")
for d, i, j in sorted(best["contacts"]):
    a, b = lab[i], flab[j]
    print(f"    cap {a[2]}{a[1]:<4d}{a[3]:<4s} - {b[0]} {b[2]}{b[1]:<4d}{b[3]:<4s} {d:5.2f} A")

cmd.delete("all")
cmd.load(FIBRIL, "hets11")
cmd.load(r["file"], "cap")
cmd.remove("cap and hydro")
Rc, tc = kabsch(ca("cap", None, CAP_CORE), ghost)
M = np.eye(4)
M[:3, :3] = Rc
M[:3, 3] = tc + best["shift"] * axis
cmd.transform_selection("cap", list(M.flatten()), homogenous=1)
if trimmed:
    cmd.remove(f"cap and resi {TAIL_DEL[0]}-{TAIL_DEL[1]}")
cmd.alter("cap", "chain='L'; segi='L'")
cmd.create("capped", "hets11 or cap")
cmd.set("pdb_use_ter_records", 1)
cmd.save(OUT, "capped")
print(f"\nsaved {OUT}: hets11 + {r['tag']} as chain L, shifted {best['shift']:.1f} A")
