#!/usr/bin/env python3
"""
free_ncaa_analysis.py -- is the free-NCAA run usable for RDF targets?
One pass over the trajectory, three questions.

    python3 free_ncaa_analysis.py solute.gro prod.xtc [--every PS] [--from NS] [--separate]

solute.gro / prod.xtc: the protein + NCAAs only (production writes only the solute group
to the .xtc). Every distance is measured across the periodic boundary where that
is shorter, so no no-jump trajectory is needed. A contact is any two heavy atoms
within 4.5 A.

  --every     ps between the frames analysed (default 100)
  --from      ns where parts 2 and 3 start (default: half the run, the same window
              part 1 judges)
  --separate  also write each graph as its own PNG (for a wiki or slides)

1. Has it settled?  Whole run; drift judged over the second half.
     NCAAs on the protein, mean NCAA radius of gyration, largest NCAA
     cluster -- each with the drift test of check_assembly_stability.py -- plus a
     check that the largest cluster stays below half the NCAAs. If it does
     not, the NCAAs are aggregating, and an NCAA-NCAA RDF from this run
     describes an aggregate, not a solution.
2. Do the NCAAs come and go?  From --from.
     per NCAA (NCAA1-NCAAn, in topology order): share of time on the protein, visits
     lasting at least MIN_VISIT_PS (shorter gaps bridged), arrivals and departures.
     A steady count on the protein does not by itself show exchange: the same few
     could stay on while others flicker. Then the largest cluster: its members
     in blocks of a twentieth of the run and how many sit on the protein; and the largest cluster
     among NCAAs off the protein (aggregation in solution).
3. Which part of the NCAA makes the contacts?  From --from.
     each NCAA split by its bonds into four parts: the backbone with its two caps, the
     atoms near the ANCHOR (where the attached part starts), the TIP group at its far
     end, and what lies between (check the sizes it prints). A part left with no atoms, as
     in a small residue, is left out of the tables and the figure. The settings RESNAME, CORE,
     ANCHOR, TIP, NEAR_BONDS, TIP_BONDS and PARTS below describe the residue; the toolbox
     fills them in from the config. Contacts per part, NCAA-NCAA
     and NCAA-protein, against the share contacts falling at random would give.

Writes (the CSVs open in Excel; none needs this script to replot):
  free_ncaa_timeseries.csv   one row per frame: adsorbed, mean Rg, cluster sizes
  free_ncaa_adsorption.csv   one row per frame, one column per NCAA: 1 = on the protein
  free_ncaa_per_ncaa.csv   one row per NCAA, over the window of parts 2-3
  free_ncaa_contacts.csv     part-by-part contact shares and ratios (part 3)
  free_ncaa_report.txt       everything printed
  free_ncaa_overview.png     all six graphs as panels a-f:
                                 a-c  NCAAs on the protein, largest cluster (and largest off
                                      the protein), mean Rg -- the judged half shaded
                                 d    each NCAA on/off the protein over time, most-bound first
                                 e    how long visits last; visits cut off by the window marked
                                 f    contact enrichment per NCAA part
  with --separate, also each graph on its own: free_ncaa_settling.png (a-c),
  free_ncaa_clusters.png (b), free_ncaa_exchange.png (d), free_ncaa_residence.png (e),
  free_ncaa_contacts.png (f)
Needs MDAnalysis. The figures are drawn from the CSV files by plot_free_ncaa.py (matplotlib).
"""
import argparse
import sys
import textwrap
from collections import deque
import numpy as np

CONTACT = 4.5          # A: heavy-atom contact
EQ_FRACTION = 0.5      # part 1 judges drift over this last fraction of the run
MIN_VISIT_PS = None    # part 2: a visit must last this long (ps); shorter gaps are bridged. None = a tenth of
                       #         the run, at most 1000 ps
MOST = 0.9             # part 2: "throughout" = more than 90% of the window
BLOCK_PS = None        # part 2: block length for the cluster table; None = a twentieth of the run (10 ns for 200 ns)
BOND = 1.95            # part 3: heavy atoms closer than this are bonded (C-S is 1.82 A)
RESNAME = "NCA"        # the free molecules are ACE-RESNAME-NME; the toolbox passes [ncaa] name
CORE = ["N", "CA", "CB", "SG", "C", "O"]     # atoms the residue shares with the amino acid it derives from
ANCHOR = "SG"          # where the attached part starts: the last of those shared atoms
TIP = None             # the atom at its far end, by name; None = the atom most bonds from ANCHOR
NEAR_BONDS = None      # part 3: atoms within this many bonds of ANCHOR are the "near anchor" part
TIP_BONDS = None       # part 3: atoms within this many bonds of TIP are the "tip" part. None = 7, or a third
                       #         of the bonds from ANCHOR to TIP when the residue is too short for 7
PARTS = ["backbone+caps", "near anchor", "tip", "middle"]     # names of the four parts, for the tables and figure
SMOOTH_NS = None       # running mean drawn over the raw traces; None = a fortieth of the run (5 ns for 200 ns)
RG_FLOOR = 0.02        # part 1: a drift in mean Rg smaller than this fraction of it counts as settled, however
                       #         small the noise (a near-rigid residue barely moves)
TINY = 1e-12

SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#1a1a19", "#52514e", "#8a8984", "#e8e7e3"
WINDOW = "#efeee9"             # shading for the window the drift test judges

LOG = []


def out(s=""):
    print(s)
    LOG.append(s)


# ----------------------------------------------------------------- helpers
def ns(ps):
    """a time in ps, written in ns with the digits a short run needs: 0.5, 12.5, 200"""
    return f"{ps / 1000:.3g}"


def settled(t, y, floor=0.0):
    """drift test of check_assembly_stability.py: a drift fails only if it is larger
    than the noise and than `floor`, and more than 3 x its own standard error"""
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    if len(y) < 4 or t.std() == 0:
        return True, 0.0, max(float(y.std()), TINY)
    a, b = np.polyfit(t, y, 1)
    res = y - (a * t + b)
    r1 = np.corrcoef(res[:-1], res[1:])[0, 1] if res.std() > 0 else 0.0
    r1 = min(max(r1 if np.isfinite(r1) else 0.0, 0.0), 0.99)
    neff = max(3.0, len(y) * (1 - r1) / (1 + r1))
    sigma = max(float(res.std()), TINY)
    se = max(sigma / np.sqrt(neff) / t.std(), TINY)
    drift = a * (t[-1] - t[0])
    return bool(abs(drift) < max(sigma, floor) or abs(a) < 3 * se), float(drift), sigma


def components(n, links):
    """connected-component label of each of n nodes"""
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for a, b in links:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    return np.array([find(i) for i in range(n)])


def mean_rg(x, starts, box):
    """mean radius of gyration of the NCAAs. Each is made whole by walking its
    heavy atoms in order and placing every atom at the minimum image of the one
    before it (consecutive heavy atoms of an Amber residue are 1-2 bonds apart)."""
    L = np.asarray(box[:3], float)
    d = np.diff(x, axis=0)
    d -= L * np.round(d / L)
    rg = []
    for a, b in zip(starts[:-1], starts[1:]):
        y = np.vstack([x[a], x[a] + np.cumsum(d[a:b - 1], axis=0)])
        rg.append(np.sqrt(((y - y.mean(0)) ** 2).sum(1).mean()))
    return float(np.mean(rg))


def element(name):
    return name.lstrip("0123456789")[:1].upper()


def guess_bonds(x, box):
    """bonded pairs among heavy atoms, by distance, minimum image"""
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
    dist[start] = 0
    q = deque([start])
    while q:
        a = q.popleft()
        for b in adj[a]:
            if dist[b] == np.inf:
                dist[b] = dist[a] + 1
                q.append(b)
    return dist


def classify(names, resnames, adj):
    """part index (0-3) of each heavy atom of one ACE-RESNAME-NME molecule, and the tip atom's name"""
    n = len(names)
    core = np.array([resnames[k] in ("ACE", "NME") or (resnames[k] == RESNAME and names[k] in set(CORE)) for k in range(n)])
    own = [k for k in range(n) if resnames[k] == RESNAME]
    a = [k for k in own if names[k] == ANCHOR]
    if len(a) != 1:
        sys.exit(f"expected one atom named {ANCHOR} (ANCHOR) in {RESNAME}, found {len(a)}")
    d_a = bfs(adj, a[0])
    if np.isinf(d_a).any():
        lost = [names[k] for k in np.where(np.isinf(d_a))[0]]
        sys.exit(f"bond guessing left {len(lost)} heavy atoms unconnected to {ANCHOR}: {' '.join(lost[:10])}")
    if TIP is None:
        far = np.where(~core, d_a, -1)
        tip = int(np.flatnonzero(far == far.max())[-1])
    else:
        tp = [k for k in own if names[k] == TIP]
        if len(tp) != 1:
            sys.exit(f"expected one atom named {TIP} (TIP) in {RESNAME}, found {len(tp)}")
        tip = tp[0]
    d_t = bfs(adj, tip)
    third = max(1, int(d_a[tip]) // 3)      # a short residue: near anchor, middle and tip a third each
    near = NEAR_BONDS if NEAR_BONDS is not None else min(7, third)
    tipb = TIP_BONDS if TIP_BONDS is not None else min(7, third)
    part = np.array([0 if core[k] else 1 if d_a[k] <= near else 2 if d_t[k] <= tipb else 3 for k in range(n)])
    return part, names[tip], near, tipb


def frame_state(prot, ncaa, ncaa_id, part_all, n_ncaa, box, finder):
    """everything the three parts need from one frame"""
    ads = np.zeros(n_ncaa, bool)
    F = np.zeros(4)
    pairs = finder(ncaa, prot, max_cutoff=CONTACT, box=box, return_distances=False)
    if len(pairs):
        ads[np.unique(ncaa_id[pairs[:, 0]])] = True
        np.add.at(F, part_all[pairs[:, 0]], 1)

    M = np.zeros((4, 4))
    links = set()
    pairs = finder(ncaa, ncaa, max_cutoff=CONTACT, box=box, return_distances=False)
    if len(pairs):
        a, b = pairs[:, 0], pairs[:, 1]
        keep = (a < b) & (ncaa_id[a] != ncaa_id[b])
        a, b = a[keep], b[keep]
        pa, pb = part_all[a], part_all[b]
        np.add.at(M, (np.minimum(pa, pb), np.maximum(pa, pb)), 1)
        la, lb = ncaa_id[a], ncaa_id[b]
        links = set(zip(np.minimum(la, lb).tolist(), np.maximum(la, lb).tolist()))

    root = components(n_ncaa, links)
    lab, cnt = np.unique(root, return_counts=True)
    big = root == lab[np.argmax(cnt)]
    off = [(i, j) for i, j in links if not ads[i] and not ads[j]]
    root2 = components(n_ncaa, off)
    free_max = int(np.unique(root2[~ads], return_counts=True)[1].max()) if (~ads).any() else 0
    return ads, big, free_max, M, F


def runs(s):
    """(start, end) of each stretch of True, end exclusive"""
    d = np.diff(np.r_[0, np.asarray(s, int), 0])
    return list(zip(np.where(d == 1)[0], np.where(d == -1)[0]))


def visits(s, g):
    """stretches on the protein lasting at least g frames, after bridging gaps
    shorter than g frames"""
    s = np.asarray(s, bool).copy()
    for a, b in runs(~s):
        if a > 0 and b < len(s) and b - a < g:
            s[a:b] = True
    for a, b in runs(s):
        if b - a < g:
            s[a:b] = False
    return runs(s)


def names_of(idx):
    return " ".join(f"NCAA{i + 1}" for i in idx) if len(idx) else "none"


def running_mean(y, w):
    if w <= 1:
        return y
    pad = np.r_[np.full(w // 2, y[0]), y, np.full(w - 1 - w // 2, y[-1])]
    return np.convolve(pad, np.ones(w) / w, mode="valid")


# ----------------------------------------------------------------- the three parts
def part1(T, ADS, BIG, RG, n_ncaa):
    cut = int(len(T) * (1 - EQ_FRACTION))
    t = T[cut:]
    series = [("NCAAs on the protein", ADS.sum(1)), ("mean NCAA Rg (A)", RG),
              ("largest cluster", BIG.sum(1))]
    out(f"\n1. has it settled?   drift judged over {ns(t[0])}-{ns(t[-1])} ns")
    out(f"   {'quantity':<24}{'mean':>9}{'drift':>10}{'sigma':>8}")
    verdicts = []
    for lab, y in series:
        floor = RG_FLOOR * abs(float(y[cut:].mean())) if "Rg" in lab else 0.0
        ok, drift, sigma = settled(t, y[cut:], floor)
        out(f"   {lab:<24}{y[cut:].mean():>9.3g}{drift:>+10.3g}{sigma:>8.3g}  "
            f"{'settled' if ok else 'NOT settled'}")
        verdicts.append((lab.split(' (')[0], ok, drift, float(y[cut:].mean())))
    big_mean = BIG.sum(1)[cut:].mean()
    agg = big_mean > n_ncaa / 2
    out(f"   largest cluster averages {big_mean:.1f} of {n_ncaa} NCAAs"
        + (" -- MORE THAN HALF: the NCAAs are aggregating" if agg else ""))
    return cut, verdicts, agg


def part2(T, ADS, BIG, FREE, dt):
    n_frames, n_ncaa = ADS.shape
    g = max(1, int(round(MIN_VISIT_PS / dt)))
    span = (T[-1] - T[0] + dt) / 1000.0
    out(f"\n2. do the NCAAs come and go?   {ns(T[0])}-{ns(T[-1])} ns, "
        f"visits must last {ns(MIN_VISIT_PS)} ns")

    frac = ADS.mean(0)
    rows, arrivals, departures, lengths, all_visits = [], 0, 0, [], []
    for k in range(n_ncaa):
        v = visits(ADS[:, k], g)
        all_visits += [(k, (b - a) * dt / 1000.0, a == 0 or b == n_frames) for a, b in v]
        arr = sum(1 for a, b in v if a > 0)
        dep = sum(1 for a, b in v if b < n_frames)
        arrivals += arr
        departures += dep
        lengths += [(b - a) * dt / 1000.0 for a, b in v if a > 0 and b < n_frames]
        longest = max(((b - a) * dt / 1000.0 for a, b in v), default=0.0)
        rows.append((k, frac[k], len(v), arr, dep, longest, BIG[:, k].mean()))

    visited = [r[0] for r in rows if r[2] > 0]
    stuck = [r[0] for r in rows if r[1] > MOST]
    never = [r[0] for r in rows if r[2] == 0]
    lab = lambda s: f"   {s:<32}"
    out(lab("on the protein on average") + f"{ADS.sum(1).mean():.2f} NCAAs")
    out(lab("visited the protein") + f"{len(visited)} of {n_ncaa}")
    out(lab(f"on for >{MOST:.0%} of the window") + f"{len(stuck)}  ({names_of(stuck)})")
    out(lab("never visited") + f"{len(never)}")
    out(lab("arrivals / departures") + f"{arrivals} / {departures} in {span:.3g} ns")
    if lengths:
        out(lab("completed visits") + f"{len(lengths)}, median {np.median(lengths):.3g} ns, "
            f"longest {max(lengths):.3g} ns")
    else:
        out(lab("completed visits") + "none (no NCAA both arrived and left inside the window)")
    out("\n   NCAAs on the protein for at least 10% of the window:")
    out(f"   {'NCAA':>7}{'on':>7}{'visits':>8}{'arrive':>8}{'leave':>7}{'longest/ns':>12}")
    for k, f, nv, arr, dep, lo, _ in sorted(rows, key=lambda r: -r[1]):
        if f >= 0.10:
            out(f"   {'NCAA' + str(k + 1):>7}{f:>7.0%}{nv:>8d}{arr:>8d}{dep:>7d}{lo:>12.3g}")

    size = BIG.sum(1)
    on = (BIG & ADS).sum(1)
    core = np.where(BIG.mean(0) > MOST)[0]
    out("\n   largest NCAA cluster")
    out(f"   mean size {size.mean():.2f}; of its members, on average {on.mean():.2f} "
        f"({(on / np.maximum(size, 1)).mean():.0%}) are on the protein")
    out(f"   in it for >{MOST:.0%} of the window: {len(core)}  ({names_of(core)})")
    if len(core):
        both = np.intersect1d(core, stuck)
        out(f"   of those, also on the protein for >{MOST:.0%}: {len(both)}  ({names_of(both)})")
    out(f"\n   {'block/ns':>10}{'size':>7}{'on protein':>11}{'same members':>14}{'off-protein max':>16}")
    prev = None
    for b in range(int(np.ceil((T[-1] - T[0] + dt) / BLOCK_PS))):
        t0 = T[0] + b * BLOCK_PS
        sel = (T >= t0) & (T < t0 + BLOCK_PS)
        if sel.sum() < 2:
            continue
        members = set(np.where(BIG[np.where(sel)[0][0]])[0])
        same = f"{len(members & prev)}/{len(members)}" if prev is not None else "-"
        prev = members
        out(f"   {t0/1000:>5.4g}-{(t0 + BLOCK_PS)/1000:<6.4g}{size[sel].mean():>7.2f}{on[sel].mean():>11.2f}"
            f"{same:>14}{FREE[sel].mean():>16.2f}")
    out("   (same members: how many of the cluster's members at the start of the block were")
    out("    also in it at the start of the previous block)")
    out(f"\n   largest cluster among NCAAs off the protein: mean {FREE.mean():.2f}, "
        f"largest seen {FREE.max()}")
    return rows, dict(arrivals=arrivals, departures=departures, stuck=stuck, visited=visited,
                      visits=all_visits)


def part3(M, F, part1_, n_frames):
    size = np.bincount(part1_, minlength=4)
    frac = size / size.sum()
    used = [a for a in range(4) if size[a]]              # a part with no atoms has no row
    tot, ft = M.sum(), F.sum()
    csv = []
    out(f"\n3. which part makes the contacts?   {n_frames} frames")
    out(f"   NCAA-NCAA: {tot / max(n_frames, 1):.0f} heavy-atom pairs per frame on average")
    w = max(22, max(len(PARTS[a] + PARTS[b]) for a, b in [(a, b) for a in used for b in used]) + 5)
    out(f"   {'parts':<{w}}{'observed':>10}{'expected':>10}{'ratio':>8}")
    E = np.zeros((4, 4))
    for a in range(4):
        for b in range(a, 4):
            E[a, b] = frac[a] ** 2 if a == b else 2 * frac[a] * frac[b]
    pairs = [(a, b) for a in used for b in used if b >= a]
    for a, b in sorted(pairs, key=lambda ab: -(M[ab] / tot / E[ab] if tot else 0)):
        obs = M[a, b] / tot if tot else 0.0
        out(f"   {PARTS[a] + ' - ' + PARTS[b]:<{w}}{obs:>10.1%}{E[a, b]:>10.1%}{obs / E[a, b]:>8.2f}")
        csv.append(("NCAA-NCAA pair", PARTS[a], PARTS[b], "", obs, E[a, b], obs / E[a, b]))

    inv = np.array([M[a, :].sum() + M[:, a].sum() for a in range(4)], float)
    inv /= 2 * tot if tot else 1
    wp = max(12, max(len(PARTS[a]) for a in used) + 2)
    out(f"\n   {'part':<{wp}}{'atoms':>7}{'share of atoms':>16}{'share of contacts':>19}{'ratio':>8}")
    for a in used:
        out(f"   {PARTS[a]:<{wp}}{size[a]:>7d}{frac[a]:>16.1%}{inv[a]:>19.1%}{inv[a] / frac[a]:>8.2f}")
        csv.append(("NCAA-NCAA part", PARTS[a], "", size[a], inv[a], frac[a], inv[a] / frac[a]))

    out(f"\n   NCAA-protein: {ft / max(n_frames, 1):.0f} heavy-atom pairs per frame on average")
    out(f"   {'part':<{wp}}{'share of atoms':>16}{'share of contacts':>19}{'ratio':>8}")
    for a in used:
        s = F[a] / ft if ft else 0.0
        out(f"   {PARTS[a]:<{wp}}{frac[a]:>16.1%}{s:>19.1%}{s / frac[a]:>8.2f}")
        csv.append(("NCAA-protein part", PARTS[a], "", size[a], s, frac[a], s / frac[a]))
    return csv


def _json(o):
    """numpy numbers, arrays and sets, for json.dump"""
    if isinstance(o, np.integer): return int(o)
    if isinstance(o, np.floating): return float(o)
    if isinstance(o, np.bool_): return bool(o)
    if isinstance(o, np.ndarray): return o.tolist()
    if isinstance(o, (set, frozenset)): return sorted(o)
    raise TypeError(type(o))


# ----------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gro")
    ap.add_argument("xtc")
    ap.add_argument("--every", type=float, default=100.0, help="ps between frames analysed (default 100)")
    ap.add_argument("--from", dest="start_ns", type=float, default=None,
                    help="ns where parts 2 and 3 start (default: half the run)")
    ap.add_argument("--separate", action="store_true",
                    help="also write each graph as its own PNG")
    args = ap.parse_args()

    import MDAnalysis as mda
    from MDAnalysis.lib.distances import capped_distance

    u = mda.Universe(args.gro, args.xtc)
    rn = [r.resname for r in u.residues]
    if "ACE" not in rn:
        sys.exit("no ACE residue found -- is this the protein + NCAA system?")
    first = rn.index("ACE")
    lres = u.residues[first:]
    if len(lres) % 3 or any(tuple(rn[first + i:first + i + 3]) != ("ACE", RESNAME, "NME")
                            for i in range(0, len(lres), 3)):
        sys.exit(f"residues after the protein are not a run of ACE-{RESNAME}-NME (is RESNAME right?)")
    n_ncaa = len(lres) // 3
    prot = u.residues[:first].atoms.select_atoms("not name H*")
    ncaa = lres.atoms.select_atoms("not name H*")
    counts = {len(lres[3 * k:3 * k + 3].atoms.select_atoms("not name H*")) for k in range(n_ncaa)}
    if len(counts) != 1:
        sys.exit("NCAAs do not all have the same number of heavy atoms")
    per = counts.pop()
    ncaa_id = np.repeat(np.arange(n_ncaa), per)
    starts = np.arange(n_ncaa + 1) * per

    # part 3's split of the NCAA, from the bonds in the first frame; the second NCAA as a check
    u.trajectory[0]
    split = []
    for k in (0, 1):
        g = ncaa[k * per:(k + 1) * per]
        split.append(classify(list(g.names), list(g.resnames), guess_bonds(g.positions, u.dimensions)))
    tip_name, near, tipb = split[0][1:]
    split = [x[0] for x in split]
    if not np.array_equal(split[0], split[1]):
        print("WARNING: NCAAs 1 and 2 split differently -- bond guessing may have failed; using NCAA 1")
    part1_ = split[0]
    part_all = np.tile(part1_, n_ncaa)
    size = np.bincount(part1_, minlength=4)

    n_traj = len(u.trajectory)
    dt_traj = float(u.trajectory.dt)
    stride = max(1, int(round(args.every / dt_traj)))
    dt = stride * dt_traj
    idx = list(range(0, n_traj, stride))
    out(f"protein: {first} residues, {prot.n_atoms} heavy atoms; NCAAs: {n_ncaa} x ACE-{RESNAME}-NME, "
        f"{per} heavy atoms each")
    out("NCAA parts: " + ", ".join(f"{PARTS[a]} {size[a]}" for a in range(4))
        + f"   (anchor {ANCHOR}, tip {tip_name}" + (", the atom most bonds from the anchor" if TIP is None else "")
        + f"; near anchor: within {near} bond{'s' * (near != 1)} of it, tip: within {tipb} bond{'s' * (tipb != 1)} of the tip)")
    out(f"trajectory: {n_traj} frames, one every {dt_traj:g} ps; analysing {len(idx)}, one every {dt:g} ps")

    T, ADS, BIG, FREE, RG, MM, FF = [], [], [], [], [], [], []
    for n, ts in enumerate(u.trajectory[::stride]):
        box = ts.dimensions
        a, b, f, M, F = frame_state(prot.positions, ncaa.positions, ncaa_id, part_all, n_ncaa, box,
                                    capped_distance)
        T.append(ts.time)
        ADS.append(a)
        BIG.append(b)
        FREE.append(f)
        MM.append(M)
        FF.append(F)
        RG.append(mean_rg(ncaa.positions, starts, box))
        if n % 100 == 0 or n == len(idx) - 1:
            print(f"\r  frame {n + 1}/{len(idx)}  ({ns(ts.time)} ns)", end="", file=sys.stderr)
    print("", file=sys.stderr)
    T, ADS, BIG, FREE, RG = map(np.asarray, (T, ADS, BIG, FREE, RG))
    MM, FF = np.asarray(MM), np.asarray(FF)
    if len(T) < 20:
        sys.exit(f"only {len(T)} frames analysed -- lower --every")
    global BLOCK_PS, SMOOTH_NS, MIN_VISIT_PS
    run_ps = float(T[-1] - T[0])
    MIN_VISIT_PS = min(1000.0, run_ps / 10) if MIN_VISIT_PS is None else float(MIN_VISIT_PS)
    BLOCK_PS = run_ps / 20 if BLOCK_PS is None else float(BLOCK_PS)
    SMOOTH_NS = run_ps / 40000 if SMOOTH_NS is None else float(SMOOTH_NS)
    out(f"time windows for this {run_ps / 1000:.4g} ns run: blocks of {BLOCK_PS / 1000:.3g} ns, smoothing {SMOOTH_NS:.3g} ns, "
        f"visits of at least {MIN_VISIT_PS / 1000:.3g} ns")

    cut, verdicts, agg = part1(T, ADS, BIG, RG, n_ncaa)
    start = args.start_ns * 1000.0 if args.start_ns is not None else T[cut]
    sel = T >= start - 1e-6
    if sel.sum() < 10:
        sys.exit(f"only {sel.sum()} frames after {ns(start)} ns -- lower --from")
    rows, ex = part2(T[sel], ADS[sel], BIG[sel], FREE[sel], dt)
    contacts = part3(MM[sel].sum(0), FF[sel].sum(0), part1_, int(sel.sum()))

    out("\nverdict")
    for lab, ok, _, _ in verdicts:
        out(f"  [{'PASS' if ok else 'FAIL'}]  {lab} has settled")
    out(f"  [{'FAIL' if agg else 'PASS'}]  NCAAs not aggregated into one cluster")
    last = T >= T[-1] - BLOCK_PS + 1e-6
    big_last = BIG.sum(1)[last].mean()
    if big_last > n_ncaa / 2:
        out(f"  [WARN]  final {BLOCK_PS / 1000:.3g} ns: the largest cluster averages {big_last:.1f} of "
            f"{n_ncaa} NCAAs, more than half")
    out(f"  [info]  from {ns(start)} ns: {len(ex['visited'])} of {n_ncaa} NCAAs visited the protein, "
        f"{ex['arrivals']} arrivals / {ex['departures']} departures, "
        f"{len(ex['stuck'])} on it for >{MOST:.0%} of the time")
    out("  take RDFs only from a stretch where all three part-1 traces are flat (panels a-c)")

    # ---- files
    np.savetxt("free_ncaa_timeseries.csv",
               np.column_stack([T / 1000.0, ADS.sum(1), RG, BIG.sum(1), (BIG & ADS).sum(1), FREE]),
               delimiter=",", comments="", fmt=["%.3f", "%d", "%.3f", "%d", "%d", "%d"],
               header="time_ns,ncaa_on_protein,mean_rg_A,largest_cluster,"
                      "largest_cluster_on_protein,largest_cluster_off_protein")
    np.savetxt("free_ncaa_adsorption.csv", np.column_stack([T / 1000.0, ADS.astype(int)]),
               delimiter=",", comments="", fmt=["%.3f"] + ["%d"] * n_ncaa,
               header=",".join(["time_ns"] + [f"NCAA{k + 1}" for k in range(n_ncaa)]))
    with open("free_ncaa_per_ncaa.csv", "w") as fh:
        fh.write("NCAA,frac_on_protein,visits,arrivals,departures,longest_visit_ns,"
                 "frac_in_largest_cluster,window_from_ns\n")
        for k, f, nv, arr, dep, lo, inbig in rows:
            fh.write(f"NCAA{k + 1},{f:.4f},{nv},{arr},{dep},{lo:.2f},{inbig:.4f},{start / 1000:.3f}\n")
    with open("free_ncaa_contacts.csv", "w") as fh:
        fh.write("kind,part_a,part_b,atoms,share_of_contacts,share_expected,ratio,window_from_ns,frames\n")
        for kind, a, b, n, obs, exp, r in contacts:
            fh.write(f"{kind},{a},{b},{n},{obs:.6f},{exp:.6f},{r:.6f},{start / 1000:.3f},{int(sel.sum())}\n")
    import json
    json.dump(dict(cut=cut, start=start, n_ncaa=n_ncaa, dt=dt, separate=bool(args.separate), verdicts=verdicts, ex=ex,
                   order=[r[0] for r in sorted(rows, key=lambda r: -r[1])],
                   SMOOTH_NS=SMOOTH_NS, MIN_VISIT_PS=MIN_VISIT_PS, MOST=MOST, BLOCK_PS=BLOCK_PS, PARTS=PARTS,
                   RESNAME=RESNAME, ANCHOR=ANCHOR, tip=tip_name),
              open("free_ncaa_meta.json", "w"), indent=1, default=_json)
    out("\nwrote free_ncaa_timeseries.csv, free_ncaa_adsorption.csv, free_ncaa_per_ncaa.csv,")
    out("      free_ncaa_contacts.csv, free_ncaa_report.txt")
    out("      free_ncaa_meta.json; draw the figures from these files with plot_free_ncaa.py")
    with open("free_ncaa_report.txt", "w") as fh:
        fh.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
