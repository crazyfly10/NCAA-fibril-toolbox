#!/usr/bin/env python3
"""
linked_eq_check.py -- from when is a run of an assembly with attached residues usable?

    python3 linked_eq_check.py solute.gro solute.xtc [every_ps] [mark_ns]

solute.gro / solute.xtc: the protein with its attached residues only, written with trjconv; water and ions
in the group are ignored. Every frame in the file is read unless every_ps is given (100 = one frame per
100 ps). mark_ns (optional): shade the figure from this time to the end.

Writes linked_eq.csv (every value, frame by frame) and linked_eq_meta.json; plot_linked_eq.py draws
linked_eq.png from them. The text report goes to the terminal.

What has to have stopped changing is the conformation of each attached residue. Three slow descriptors per
residue stand in for it:
  Rg        radius of gyration of its heavy atoms
  reach     distance from its ANCHOR atom to its TIP atom: how far it stretches
  on        its atoms beyond the anchor region within CONTACT of the protein: lying on the surface, or out
            in the water
and, alongside, whether the tip itself touches the protein.

Nothing here knows the chemistry. The residue is described by the settings below:
  RESNAME     its residue name
  CORE        the atoms it shares with the amino acid it derives from (they sit in the protein; not counted)
  ANCHOR      where the attached part starts: the last of those shared atoms (SG for a cysteine adduct)
  TIP         the atom at its far end, by name; None = the atom most bonds away from ANCHOR
  NEAR_BONDS  atoms within this many bonds of ANCHOR are too close to the protein to say anything, and are
              left out of "on"
  TIP_BONDS   atoms within this many bonds of TIP count as the tip
The toolbox fills these in from the config and from what ncaa-param recorded.

Where to start: the number of statistically independent samples left after a start time (Chodera, JCTC 12,
1799 (2016)) is low for a start still inside the initial relaxation, because the drift makes the data look
strongly correlated, and falls again for late starts, which throw good data away. For each descriptor the
script takes the earliest start that keeps at least 90% of the best number; the latest of these over all
descriptors and residues is the one to use. It then checks that the first and second half of the kept
stretch agree within their statistical uncertainty, and says how many independent samples each descriptor
has: a few dozen or more is comfortable, and under ten the run is too short to judge.

Needs MDAnalysis.
"""
import json
import sys
from collections import deque
import numpy as np

RESNAME = "CYL"
CORE = ["N", "CA", "C", "O", "CB", "SG"]
ANCHOR = "SG"
TIP = None             # atom name, or None for the atom most bonds from ANCHOR
NEAR_BONDS = 7
TIP_BONDS = 7
CONTACT = 4.5          # A
BOND = 1.95            # A: heavy-atom pairs closer than this are bonded
SKIP_RES = {"WAT", "SOL", "HOH", "Na+", "Cl-", "NA", "CL", "K+", "K"}      # not protein; RESNAME is added
ON_MIN = 5             # atoms on the protein to count the residue as lying on it (fewer for a small residue)
Z_MAX = 3.0            # halves may differ by this many standard errors
SETTLE_FRAC = 0.10     # contacts count as settled once their running mean stays within 10% of their final value
MIN_SAMPLES = 10       # fewer independent samples than this: too short to judge
# time windows in ns; None = a fraction of the run, which gives 5 / 20 / 20 / 5 ns for a 200 ns run
START_STEP_NS = None   # candidate start times                         (run / 40)
BLOCK_NS = None        # block averages in the report                  (run / 10)
TAIL_NS = None         # the stretch that defines the final value      (run / 10)
SMOOTH_NS = None       # running mean for the thick lines in the figure (run / 40)
SITE_NAMES = {}        # label -> name for the figure, e.g. {"CYL22": "chain A, residue 239"}
SITE_SHORT = {}        # label -> short name, e.g. {"CYL22": "chain A"}


def element(name):
    return name.lstrip("0123456789")[:1].upper()


def guess_bonds(x, box):
    L = np.asarray(box[:3], float)
    d = x[:, None, :] - x[None, :, :]
    d -= L * np.round(d / L)
    r = np.sqrt((d ** 2).sum(-1))
    i, j = np.where((r > 0.1) & (r < BOND))
    adj = [[] for _ in range(len(x))]
    for a, b in zip(i, j):
        adj[a].append(b)
    return adj


def bfs(adj, start):
    dist = np.full(len(adj), np.inf)
    parent = np.full(len(adj), -1)
    order = [start]
    dist[start] = 0
    q = deque([start])
    while q:
        a = q.popleft()
        for b in adj[a]:
            if dist[b] == np.inf:
                dist[b] = dist[a] + 1
                parent[b] = a
                order.append(b)
                q.append(b)
    return dist, parent, order


class Attached:
    """one attached residue: its anchor, its tip, which atoms are the distal part and which the tip group, and a
    bond tree for making it whole across the periodic boundary"""

    def __init__(self, names, x0, box):
        n = len(names)
        adj = guess_bonds(x0, box)
        a = [k for k in range(n) if names[k] == ANCHOR]
        if len(a) != 1:
            sys.exit(f"{RESNAME}: expected one atom named {ANCHOR} (ANCHOR), found {len(a)}; heavy atoms: {' '.join(names)}")
        self.anchor = a[0]
        d_a, self.parent, self.order = bfs(adj, self.anchor)
        if np.isinf(d_a).any():
            lost = [names[k] for k in np.where(np.isinf(d_a))[0]]
            sys.exit(f"{RESNAME}: bond guessing left atoms unconnected to {ANCHOR}: {' '.join(lost[:10])}")
        core = np.array([names[k] in set(CORE) for k in range(n)])
        if TIP is None:
            far = np.where(~core, d_a, -1)
            self.tip = int(np.flatnonzero(far == far.max())[-1])
        else:
            tp = [k for k in range(n) if names[k] == TIP]
            if len(tp) != 1:
                sys.exit(f"{RESNAME}: expected one atom named {TIP} (TIP), found {len(tp)}; heavy atoms: {' '.join(names)}")
            self.tip = tp[0]
        self.tip_name = names[self.tip]
        d_t = bfs(adj, self.tip)[0]
        near = NEAR_BONDS
        while True:
            part = np.array([0 if core[k] else 1 if d_a[k] <= near else 2 if d_t[k] <= TIP_BONDS else 3 for k in range(n)])
            if (part >= 2).any() or near == 0:
                break
            near = 0                       # a small residue: everything is near the anchor, so count all of it
        self.near_used = near
        self.sizes = np.bincount(part, minlength=4)
        self.distal = np.where(part >= 2)[0]
        self.tipgroup = np.where(part == 2)[0]
        if len(self.distal) == 0:
            sys.exit(f"{RESNAME}: no atoms beyond CORE ({' '.join(CORE)}); nothing to follow")

    def whole(self, x, box):
        L = np.asarray(box[:3], float)
        y = x.copy()
        for a in self.order[1:]:
            p = self.parent[a]
            d = y[a] - y[p]
            y[a] = y[p] + d - L * np.round(d / L)
        return y


def frame_values(link, x, prot, box, finder):
    """Rg, reach (anchor to tip), distal atoms on the protein, tip on the protein"""
    y = link.whole(x, box)
    rg = float(np.sqrt(((y - y.mean(0)) ** 2).sum(1).mean()))
    reach = float(np.linalg.norm(y[link.tip] - y[link.anchor]))
    pairs = finder(x[link.distal], prot, max_cutoff=CONTACT, box=box, return_distances=False)
    touching = np.unique(pairs[:, 0]) if len(pairs) else np.array([], int)
    on = len(touching)
    tip_on = bool(np.isin(link.distal[touching], link.tipgroup).any()) if on else False
    return rg, reach, on, tip_on


def stat_ineff(y, mintime=3):
    """statistical inefficiency g (in frames): the data hold len(y)/g independent samples"""
    y = np.asarray(y, float)
    n = len(y)
    dy = y - y.mean()
    var = dy.var()
    if n < 3 or var == 0:
        return 1.0
    g = 1.0
    for lag in range(1, n - 1):
        c = np.dot(dy[:n - lag], dy[lag:]) / ((n - lag) * var)
        if c <= 0 and lag > mintime:
            break
        g += 2.0 * c * (1.0 - lag / n)
    return max(1.0, g)


def detect_start(y, step, keep=0.9):
    """earliest start frame (multiple of step, within the first half) whose number
    of independent samples left is within `keep` of the best start's. The best
    start alone is noisy once the data are stationary; the earliest near-best one
    marks where the initial relaxation ends."""
    cand = []
    for i0 in range(0, len(y) // 2 + 1, step):
        g = stat_ineff(y[i0:])
        cand.append((i0, g, (len(y) - i0) / g))
    top = max(c[2] for c in cand)
    return next(c for c in cand if c[2] >= keep * top)


def halves(y):
    """difference of the two halves' means in standard errors, the means, and the
    number of independent samples in the whole stretch"""
    h = len(y) // 2
    a, b = y[:h], y[h:]
    ga, gb = stat_ineff(a), stat_ineff(b)
    se = np.sqrt(a.var() * ga / len(a) + b.var() * gb / len(b))
    z = abs(a.mean() - b.mean()) / se if se > 0 else 0.0
    return z, a.mean(), b.mean(), len(y) / stat_ineff(y)


def runs(s):
    d = np.diff(np.r_[0, np.asarray(s, int), 0])
    return list(zip(np.where(d == 1)[0], np.where(d == -1)[0]))


def switches(on, g):
    """number of lie-down / lift-off events, ignoring stretches shorter than g frames"""
    s = np.asarray(on, bool).copy()
    for a, b in runs(~s):
        if a > 0 and b < len(s) and b - a < g:
            s[a:b] = True
    for a, b in runs(s):
        if b - a < g:
            s[a:b] = False
    return max(0, len(np.flatnonzero(np.diff(s.astype(int)))))


def analyse(t, data, labels, on_min, out=print):
    """t: times (ps); data[k]: array (frames, 4) for residue k. Returns the start time in ns and whether the
    kept stretch can be trusted."""
    dt = float(np.median(np.diff(t)))
    names = ("Rg", "reach", "on")
    out(f"{len(t)} frames, {t[0]/1000:.0f}-{t[-1]/1000:.0f} ns, one every {dt:.0f} ps\n")

    for k, lab in enumerate(labels):
        d = data[k]
        on = d[:, 2] >= on_min
        out(f"{lab}: {d[:, 2].mean():.1f} distal atoms on the protein on average; lying on it "
            f"{on.mean():.0%} of the time, tip on it {d[:, 3].mean():.0%}; "
            f"{switches(on, max(1, int(round(1000 / dt))))} lie-down/lift-off events")
        out(f"   {'block/ns':>12}{'Rg/A':>8}{'reach/A':>9}{'on':>7}{'tip on':>9}")
        nb = int(np.ceil((t[-1] - t[0] + dt) / (BLOCK_NS * 1000)))
        for b in range(nb):
            t0 = t[0] + b * BLOCK_NS * 1000
            s = (t >= t0) & (t < t0 + BLOCK_NS * 1000)
            if s.sum() < 2:
                continue
            out(f"   {t0/1000:>5.4g}-{(t0 + BLOCK_NS*1000)/1000:<6.4g}{d[s, 0].mean():>8.2f}"
                f"{d[s, 1].mean():>9.2f}{d[s, 2].mean():>7.1f}{d[s, 3].mean():>9.0%}")
        out("")

    labs = [f"{lab} {n}" for lab in labels for n in names]
    series = [data[k][:, j] for k in range(len(labels)) for j in range(3)]
    step = max(1, int(round(START_STEP_NS * 1000 / dt)))
    out("earliest start keeping at least 90% of the most independent samples, per descriptor")
    out(f"   {'descriptor':<16}{'start/ns':>9}{'corr. time/ns':>15}{'independent samples':>21}")
    starts = []
    for lab, y in zip(labs, series):
        i0, g, neff = detect_start(y, step)
        starts.append(i0)
        out(f"   {lab:<16}{t[i0]/1000 - t[0]/1000:>9.4g}{g * dt / 1000:>15.2f}{neff:>21.0f}")
    i0 = max(starts)
    best = (t[i0] - t[0]) / 1000
    out(f"\nlatest of these: {best:.4g} ns -> keep {(t[-1] - t[i0]) / 1000:.4g} ns")

    out(f"\nfirst half vs second half of the kept stretch (difference in standard errors; under {Z_MAX:g} = agree)")
    out(f"   {'descriptor':<16}{'1st half':>10}{'2nd half':>10}{'difference':>12}{'independent samples':>21}")
    worst, fewest = 0.0, np.inf
    for lab, y in zip(labs, series):
        z, m1, m2, neff = halves(y[i0:])
        worst, fewest = max(worst, z), min(fewest, neff)
        out(f"   {lab:<16}{m1:>10.2f}{m2:>10.2f}{z:>12.1f}{neff:>21.0f}" + ("   DIFFER" if z >= Z_MAX else "")
            + ("   too few" if neff < MIN_SAMPLES else ""))
    out("")
    if fewest < MIN_SAMPLES:
        out(f"the run is too short to judge: a descriptor has only {fewest:.0f} independent samples in the kept stretch "
            f"(under {MIN_SAMPLES}). The start time above is not a finding; run for longer")
        return best, False
    if worst < Z_MAX:
        out(f"use {best:.4g} ns to the end for the bonded distributions")
        return best, True
    out(f"the halves still differ after {best:.4g} ns: the conformations have not settled, "
        "or switch too rarely for this run to average over")
    return best, False


def write_csv(t, data, labels, path="linked_eq.csv"):
    cols = ["time_ns"] + [f"{lab}_{q}" for lab in labels
                          for q in ("rg_A", "reach_A", "on_atoms", "tip_on")]
    table = [np.asarray(t) / 1000.0] + [data[k][:, j] for k in range(len(labels)) for j in range(4)]
    np.savetxt(path, np.column_stack(table), delimiter=",", header=",".join(cols),
               comments="", fmt="%.4f")
    print(f"wrote {path} ({len(t)} frames, {len(cols)} columns)")


def running_mean(y, w):
    if w <= 1:
        return np.asarray(y, float)
    k = np.ones(w) / w
    pad = np.r_[np.full(w // 2, y[0]), y, np.full(w - 1 - w // 2, y[-1])]
    return np.convolve(pad, k, mode="valid")


def bound_since(t, b, g, cut=0.5):
    """time from which the tip stays on the protein to the end (gaps shorter than g frames bridged), or None if
    it is off at the end"""
    s = np.asarray(b, float) >= cut
    for a, c in runs(~s):
        if a > 0 and c < len(s) and c - a < g:
            s[a:c] = True
    if not s[-1]:
        return None
    off = np.flatnonzero(~s)
    return float(t[0] if len(off) == 0 else t[off[-1] + 1])


def settled_from(t, d, w, g, cut=0.5):
    """When has one residue settled into its bound state? Two conditions, both needed: the tip has stuck for
    good (bound_since), AND the rest has stopped flattening onto the surface -- the running mean of the contact
    count stays within SETTLE_FRAC of its final value from then to the end. Returns (bind, contacts, settled) in
    ps, None where unmet."""
    bind = bound_since(t, d[:, 3], g, cut)
    t_ns = np.asarray(t, float) / 1000.0
    final = d[t_ns >= t_ns[-1] - TAIL_NS, 2].mean()
    sm = running_mean(d[:, 2], w)
    out = np.abs(sm - final) > SETTLE_FRAC * max(final, 1.0)
    contacts = None if out[-1] else float(t[0] if not out.any() else t[np.flatnonzero(out)[-1] + 1])
    # the final value is taken from the last TAIL_NS, so entering the band inside that stretch proves nothing
    if contacts is not None and contacts > t[-1] - TAIL_NS * 1000.0:
        contacts = None
    settled = None if bind is None or contacts is None else max(bind, contacts)
    return bind, contacts, settled


def main():
    global START_STEP_NS, BLOCK_NS, TAIL_NS, SMOOTH_NS
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import capped_distance

    every = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0     # 0 = every frame in the file
    u = mda.Universe(sys.argv[1], sys.argv[2])
    att = [r for r in u.residues if r.resname == RESNAME]
    if not att:
        sys.exit(f"no {RESNAME} residue found -- is this the system with the residue attached, and is RESNAME right?")
    skip = set(SKIP_RES) | {RESNAME}
    prot = u.atoms[[a.index for r in u.residues if r.resname not in skip for a in r.atoms
                    if not a.name.startswith("H")]]
    heavy = [r.atoms[np.array([not a.name.startswith("H") for a in r.atoms], bool)] for r in att]

    u.trajectory[0]
    links = [Attached(list(h.names), h.positions.copy(), u.dimensions) for h in heavy]
    labels = [f"{RESNAME}{r.resid}" for r in att]
    lk = links[0]
    n_distal = int(len(lk.distal))
    on_min = int(min(ON_MIN, max(1, n_distal // 2)))
    print(f"protein: {prot.n_atoms} heavy atoms; {len(att)} attached {RESNAME}: " + ", ".join(labels))
    print(f"each {RESNAME}: anchor {ANCHOR}, tip {lk.tip_name}" + (" (the atom most bonds from the anchor)" if TIP is None else "")
          + f"; {n_distal} distal atoms followed ({len(lk.tipgroup)} of them the tip group), "
          f"{int(lk.sizes[1])} near the anchor left out; lying on the protein = {on_min} or more distal atoms in contact")
    if lk.near_used != NEAR_BONDS:
        print(f"   every atom of {RESNAME} is within {NEAR_BONDS} bonds of {ANCHOR}, so none is left out as near the anchor")

    t, rows, last = [], [[] for _ in links], -1e30
    for ts in u.trajectory:
        if ts.time - last < every - 1e-6:
            continue
        last = ts.time
        t.append(ts.time)
        ppos = prot.positions
        for k, (h, l) in enumerate(zip(heavy, links)):
            rows[k].append(frame_values(l, h.positions, ppos, ts.dimensions, capped_distance))
    t = np.array(t)
    if len(t) < 20:
        sys.exit(f"only {len(t)} frames read: too few to analyse")
    data = [np.array(r, float) for r in rows]
    run_ns = (t[-1] - t[0]) / 1000.0
    START_STEP_NS = run_ns / 40 if START_STEP_NS is None else float(START_STEP_NS)
    BLOCK_NS = run_ns / 10 if BLOCK_NS is None else float(BLOCK_NS)
    TAIL_NS = run_ns / 10 if TAIL_NS is None else float(TAIL_NS)
    SMOOTH_NS = run_ns / 40 if SMOOTH_NS is None else float(SMOOTH_NS)
    print(f"time windows for this {run_ns:.4g} ns run: start candidates every {START_STEP_NS:.3g} ns, blocks of {BLOCK_NS:.3g} ns, "
          f"final value over the last {TAIL_NS:.3g} ns, smoothing {SMOOTH_NS:.3g} ns\n")
    best, ok = analyse(t, data, labels, on_min)
    print()
    dt = float(np.median(np.diff(t)))
    w = max(1, int(round(SMOOTH_NS * 1000 / dt)))
    g = max(1, int(round(1000 / dt)))
    print(f"adsorption: tip on the protein for good from / contact count within "
          f"{SETTLE_FRAC:.0%} of its last-{TAIL_NS:.3g}-ns value from")
    summary = []
    for lab, d in zip(labels, data):
        b, c, _ = settled_from(t, d, w, g)
        f = lambda x: f"{x / 1000:.4g} ns" if x is not None else "not by the end"
        print(f"   {lab:<8} tip from {f(b):<16} contacts level off from {f(c)}")
        summary.append(dict(label=lab, lying=float((d[:, 2] >= on_min).mean()), tip_on=float(d[:, 3].mean()),
                            tip_from_ns=None if b is None else b / 1000, contacts_from_ns=None if c is None else c / 1000))
    print()
    write_csv(t, data, labels)
    mark_ns = float(sys.argv[4]) if len(sys.argv) > 4 else None
    json.dump(dict(labels=labels, n_distal=n_distal, mark_ns=mark_ns, CONTACT=CONTACT, RESNAME=RESNAME, ANCHOR=ANCHOR,
                   tip=lk.tip_name, tip_atoms=int(len(lk.tipgroup)), on_min=on_min, start_ns=best, usable=bool(ok),
                   summary=summary, SITE_NAMES=SITE_NAMES, SITE_SHORT=SITE_SHORT, SMOOTH_NS=SMOOTH_NS, TAIL_NS=TAIL_NS),
              open("linked_eq_meta.json", "w"), indent=1)
    print("wrote linked_eq_meta.json  (plot linked_eq.csv with plot_linked_eq.py)")


if __name__ == "__main__":
    main()
