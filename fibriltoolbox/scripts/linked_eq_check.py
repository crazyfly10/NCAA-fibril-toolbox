#!/usr/bin/env python3
"""
linked_eq_check.py -- from when is the linked run usable for the bonded targets?

    python3 linked_eq_check.py solute.gro prot_cyl.xtc [every_ps] [mark_ns]

solute.gro / prot_cyl.xtc: fibril + attached linkers only (Prot_cyl), written with
trjconv; ions in the group are ignored. Every frame in the file is read unless
every_ps is given (e.g. 100 = one frame per 100 ps). Production saves a frame
every 10 ps, so write prot_cyl.xtc without -dt to get all of them.
mark_ns (optional): shade the figure from this time to the end.

Writes linked_eq.csv (every value, frame by frame); plot_linked_eq.py draws linked_eq.png from it (the
figure); the text report goes to the terminal.

The bonded distributions come from the attached CYL linkers, so what has to have
stopped changing is the linkers' conformation. Three slow descriptors per linker
stand in for it (the bond, angle and dihedral distributions of the coarse-grained
beads settle at least as fast as these):
  Rg        radius of gyration of the linker's heavy atoms
  reach     SG to biotin S distance: how far the linker stretches
  on        PEG + biotin heavy atoms within 4.5 A of the fibril: lying on the
            surface or out in the water

Where to start: the number of statistically independent samples left after a
start time (Chodera, JCTC 12, 1799 (2016)) is low for a start still inside the
initial relaxation, because the drift makes the data look strongly correlated,
and falls again for late starts, which throw good data away. For each descriptor
the script takes the earliest start that keeps at least 90% of the best number;
the latest of these over all descriptors and linkers is the one to use. It
then checks that the first and second half of the kept stretch agree within
their statistical uncertainty, and says how many independent samples each
descriptor has: a few dozen or more is comfortable.

It also reports how often each linker went on and off the fibril: a linker that
never left the surface gives bonded distributions for the lying-down state only.

Needs MDAnalysis.
"""
import sys
from collections import deque
import numpy as np

CONTACT = 4.5          # A
BOND = 1.95            # A: heavy-atom pairs closer than this are bonded
REACH = 7              # bonds from SG / biotin S that count as spacer / biotin
BACKBONE = {"N", "CA", "CB", "SG", "C", "O"}
SKIP_RES = {"CYL", "WAT", "SOL", "HOH", "Na+", "Cl-", "NA", "CL", "K+", "K"}
START_STEP_NS = 5.0    # candidate start times
Z_MAX = 3.0            # halves may differ by this many standard errors
BLOCK_NS = 20.0
ON_MIN = 5             # PEG+biotin atoms on the fibril to count the linker as lying on it
SETTLE_FRAC = 0.10     # contacts count as settled once their running mean stays within
                       # 10% of their final value (the mean of the last TAIL_NS)
TAIL_NS = 20.0
SMOOTH_NS = 5.0        # running-mean window for the thick lines in the figure

# names for the figure: residue numbers here run 1..N straight through the system
# (72 residues per HET-s chain, so chain A residue 239 is CYL22 and the cap's
# last residue, L130, is CYL915)
SITE_NAMES = {"CYL22": "chain A, residue 239", "CYL915": "cap, residue 130"}
SITE_SHORT = {"CYL22": "chain A", "CYL915": "cap"}


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


class Linker:
    """one attached CYL: which heavy atoms are SG, biotin S, PEG + biotin, and a
    bond tree for making it whole across the periodic boundary"""

    def __init__(self, names, x0, box):
        n = len(names)
        adj = guess_bonds(x0, box)
        sg = [k for k in range(n) if names[k] == "SG"]
        sb = [k for k in range(n) if element(names[k]) == "S" and names[k] != "SG"]
        if len(sg) != 1 or len(sb) != 1:
            sys.exit(f"CYL: expected one SG and one other sulfur, found {len(sg)} and {len(sb)}")
        self.sg, self.sb = sg[0], sb[0]
        d_sg, self.parent, self.order = bfs(adj, self.sg)
        if np.isinf(d_sg).any():
            lost = [names[k] for k in np.where(np.isinf(d_sg))[0]]
            sys.exit(f"CYL: bond guessing left atoms unconnected to SG: {' '.join(lost[:10])}")
        d_sb = bfs(adj, self.sb)[0]
        part = np.array([0 if names[k] in BACKBONE else 1 if d_sg[k] <= REACH
                         else 2 if d_sb[k] <= REACH else 3 for k in range(n)])
        self.sizes = np.bincount(part, minlength=4)
        self.distal = np.where(part >= 2)[0]          # biotin + PEG
        self.biotin = np.where(part == 2)[0]

    def whole(self, x, box):
        L = np.asarray(box[:3], float)
        y = x.copy()
        for a in self.order[1:]:
            p = self.parent[a]
            d = y[a] - y[p]
            y[a] = y[p] + d - L * np.round(d / L)
        return y


def frame_values(link, x, fib, box, finder):
    """Rg, reach (SG-biotin S), PEG+biotin atoms on the fibril, biotin on the fibril"""
    y = link.whole(x, box)
    rg = float(np.sqrt(((y - y.mean(0)) ** 2).sum(1).mean()))
    reach = float(np.linalg.norm(y[link.sb] - y[link.sg]))
    pairs = finder(x[link.distal], fib, max_cutoff=CONTACT, box=box, return_distances=False)
    touching = np.unique(pairs[:, 0]) if len(pairs) else np.array([], int)
    on = len(touching)
    biotin_on = bool(np.isin(link.distal[touching], link.biotin).any()) if on else False
    return rg, reach, on, biotin_on


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


def analyse(t, data, labels, out=print):
    """t: times (ps); data[k]: array (frames, 4) for linker k"""
    dt = float(np.median(np.diff(t)))
    names = ("Rg", "reach", "on")
    out(f"{len(t)} frames, {t[0]/1000:.0f}-{t[-1]/1000:.0f} ns, one every {dt:.0f} ps\n")

    for k, lab in enumerate(labels):
        d = data[k]
        on = d[:, 2] >= ON_MIN
        out(f"{lab}: {d[:, 2].mean():.1f} PEG+biotin atoms on the fibril on average; lying on it "
            f"{on.mean():.0%} of the time, biotin on it {d[:, 3].mean():.0%}; "
            f"{switches(on, max(1, int(round(1000 / dt))))} lie-down/lift-off events")
        out(f"   {'block/ns':>10}{'Rg/A':>8}{'reach/A':>9}{'on':>7}{'biotin on':>11}")
        nb = int(np.ceil((t[-1] - t[0] + dt) / (BLOCK_NS * 1000)))
        for b in range(nb):
            t0 = t[0] + b * BLOCK_NS * 1000
            s = (t >= t0) & (t < t0 + BLOCK_NS * 1000)
            if s.sum() < 2:
                continue
            out(f"   {t0/1000:>4.0f}-{(t0 + BLOCK_NS*1000)/1000:<5.0f}{d[s, 0].mean():>8.2f}"
                f"{d[s, 1].mean():>9.2f}{d[s, 2].mean():>7.1f}{d[s, 3].mean():>11.0%}")
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
        out(f"   {lab:<16}{t[i0]/1000 - t[0]/1000:>9.0f}{g * dt / 1000:>15.2f}{neff:>21.0f}")
    i0 = max(starts)
    best = (t[i0] - t[0]) / 1000
    out(f"\nlatest of these: {best:.0f} ns -> keep {(t[-1] - t[i0]) / 1000:.0f} ns")

    out(f"\nfirst half vs second half of the kept stretch (difference in standard errors; under {Z_MAX:g} = agree)")
    out(f"   {'descriptor':<16}{'1st half':>10}{'2nd half':>10}{'difference':>12}{'independent samples':>21}")
    worst = 0.0
    for lab, y in zip(labs, series):
        z, m1, m2, neff = halves(y[i0:])
        worst = max(worst, z)
        out(f"   {lab:<16}{m1:>10.2f}{m2:>10.2f}{z:>12.1f}{neff:>21.0f}" + ("   DIFFER" if z >= Z_MAX else ""))
    out("")
    if worst < Z_MAX:
        out(f"use {best:.0f} ns to the end for the bonded distributions")
    else:
        out(f"the halves still differ after {best:.0f} ns: the linker conformations have not settled, "
            "or switch too rarely for this run to average over")
    return best


def write_csv(t, data, labels, path="linked_eq.csv"):
    cols = ["time_ns"] + [f"{lab}_{q}" for lab in labels
                          for q in ("rg_A", "reach_A", "on_atoms", "biotin_on")]
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
    """time from which the biotin stays on the fibril to the end (gaps shorter
    than g frames bridged), or None if it is off at the end"""
    s = np.asarray(b, float) >= cut
    for a, c in runs(~s):
        if a > 0 and c < len(s) and c - a < g:
            s[a:c] = True
    if not s[-1]:
        return None
    off = np.flatnonzero(~s)
    return float(t[0] if len(off) == 0 else t[off[-1] + 1])


def settled_from(t, d, w, g, cut=0.5):
    """When has one linker settled into its bound state? Two conditions, both
    needed: the biotin has stuck for good (bound_since), AND the PEG has stopped
    flattening onto the surface -- the running mean of the contact count stays
    within SETTLE_FRAC of its final value from then to the end. Biotin binding
    alone is not enough: the PEG can keep spreading over the fibril for tens of
    ns after it, and that changes the linker's shape, which is what the bonded
    targets record. Returns (bind, contacts, settled) in ps, None where unmet."""
    bind = bound_since(t, d[:, 3], g, cut)
    t_ns = np.asarray(t, float) / 1000.0
    final = d[t_ns >= t_ns[-1] - TAIL_NS, 2].mean()
    sm = running_mean(d[:, 2], w)
    out = np.abs(sm - final) > SETTLE_FRAC * max(final, 1.0)
    contacts = None if out[-1] else float(t[0] if not out.any() else t[np.flatnonzero(out)[-1] + 1])
    # the final value is taken from the last TAIL_NS, so entering the band inside
    # that stretch proves nothing: the curve could still be climbing slowly
    if contacts is not None and contacts > t[-1] - TAIL_NS * 1000.0:
        contacts = None
    settled = None if bind is None or contacts is None else max(bind, contacts)
    return bind, contacts, settled


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import capped_distance

    every = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0     # 0 = every frame in the file
    u = mda.Universe(sys.argv[1], sys.argv[2])
    cyl = [r for r in u.residues if r.resname == "CYL"]
    if not cyl:
        sys.exit("no CYL residue found -- is this the linked system?")
    fib = u.atoms[[a.index for r in u.residues if r.resname not in SKIP_RES for a in r.atoms
                   if not a.name.startswith("H")]]
    heavy = [r.atoms[np.array([not a.name.startswith("H") for a in r.atoms], bool)] for r in cyl]

    u.trajectory[0]
    links = [Linker(list(h.names), h.positions.copy(), u.dimensions) for h in heavy]
    labels = [f"CYL{r.resid}" for r in cyl]
    print(f"fibril: {fib.n_atoms} heavy atoms; {len(cyl)} attached linkers: "
          + ", ".join(f"{lab} ({lk.sizes[2] + lk.sizes[3]} PEG+biotin atoms)" for lab, lk in zip(labels, links)))

    t, rows, last = [], [[] for _ in links], -1e30
    for ts in u.trajectory:
        if ts.time - last < every - 1e-6:
            continue
        last = ts.time
        t.append(ts.time)
        fpos = fib.positions
        for k, (h, lk) in enumerate(zip(heavy, links)):
            rows[k].append(frame_values(lk, h.positions, fpos, ts.dimensions, capped_distance))
    t = np.array(t)
    data = [np.array(r, float) for r in rows]
    analyse(t, data, labels)
    print()
    dt = float(np.median(np.diff(t)))
    w = max(1, int(round(SMOOTH_NS * 1000 / dt)))
    g = max(1, int(round(1000 / dt)))
    print(f"adsorption: biotin on the fibril for good from / contact count within "
          f"{SETTLE_FRAC:.0%} of its last-{TAIL_NS:g}-ns value from")
    for lab, d in zip(labels, data):
        b, c, _ = settled_from(t, d, w, g)
        f = lambda x: f"{x / 1000:.0f} ns" if x is not None else "not by the end"
        print(f"   {lab:<8} biotin from {f(b):<16} contacts level off from {f(c)}")
    print()
    write_csv(t, data, labels)
    mark_ns = float(sys.argv[4]) if len(sys.argv) > 4 else None
    import json
    json.dump(dict(labels=labels, n_distal=int(links[0].sizes[2] + links[0].sizes[3]), mark_ns=mark_ns, CONTACT=CONTACT,
                   SITE_NAMES=SITE_NAMES, SITE_SHORT=SITE_SHORT, SMOOTH_NS=SMOOTH_NS, TAIL_NS=TAIL_NS),
              open("linked_eq_meta.json", "w"), indent=1)
    print("wrote linked_eq_meta.json  (plot linked_eq.csv with plot_linked_eq.py)")


if __name__ == "__main__":
    main()
