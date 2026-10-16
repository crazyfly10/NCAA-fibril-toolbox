#!/usr/bin/env python3
"""
Draw the free-NCAA figures from the files free_ncaa_analysis.py writes:

    free_ncaa_timeseries.csv   one row per frame: adsorbed, mean Rg, cluster sizes
    free_ncaa_adsorption.csv   one row per frame, one column per NCAA: 1 = on the protein
    free_ncaa_per_ncaa.csv   one row per NCAA
    free_ncaa_contacts.csv     part-by-part contact shares and ratios
    free_ncaa_meta.json        the window, the verdicts, the visits and the settings the analysis used

Usage:
    python3 plot_free_ncaa.py                 # free_ncaa_overview.png
    python3 plot_free_ncaa.py --separate      # also each graph as its own PNG

This script only draws. Change it freely: the numbers and the report come from free_ncaa_analysis.py.
"""

import csv
import json
import sys
import textwrap
import numpy as np

meta = json.load(open("free_ncaa_meta.json"))
SMOOTH_NS, MIN_VISIT_PS, MOST, BLOCK_PS, PARTS = (meta[k] for k in ("SMOOTH_NS", "MIN_VISIT_PS", "MOST", "BLOCK_PS", "PARTS"))
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#1a1a19", "#52514e", "#8a8984", "#e8e7e3"
WINDOW = "#efeee9"             # shading for the window the drift test judges


def ns(ps):
    """a time in ps, written in ns with the digits a short run needs: 0.5, 12.5, 200"""
    return f"{ps / 1000:.3g}"


def running_mean(y, w):
    if w <= 1:
        return y
    pad = np.r_[np.full(w // 2, y[0]), y, np.full(w - 1 - w // 2, y[-1])]
    return np.convolve(pad, np.ones(w) / w, mode="valid")


# ----------------------------------------------------------------- figures
SERIES1, SERIES2 = "#2a78d6", "#eb6834"      # categorical slots 1 and 2 (blue, orange)
NOTE_CONTACT = "On the protein: any NCAA heavy atom within 4.5 Å of a protein heavy atom. Cluster: NCAAs " \
               "joined by such contacts."


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                         "axes.linewidth": 0.8, "xtick.color": MUTED, "ytick.color": MUTED,
                         "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
                         "xtick.major.size": 3, "ytick.major.size": 3, "text.color": INK,
                         "legend.frameon": False, "legend.fontsize": 9.5})
    return plt


def _header(fig, title, stand, note=None, left=0.075, wrap=11.5):
    """title and standfirst at the top, a note at the bottom; returns the figure
    fraction the plots may reach up to"""
    W, H = fig.get_size_inches()
    y = 1 - 0.28 / H
    fig.text(left, y, title, fontsize=15, fontweight="bold", color=INK, va="top")
    y -= 0.38 / H
    for line in textwrap.wrap(stand, int(W * wrap)):
        fig.text(left, y, line, fontsize=10, color=INK2, va="top")
        y -= 0.19 / H
    if note:
        yn = 0.12 / H
        for line in reversed(textwrap.wrap(note, int(W * 13.5))):
            fig.text(left, yn, line, fontsize=8.5, color=MUTED, va="bottom")
            yn += 0.16 / H
    return y - 0.22 / H


def _axes(ax, grid="y"):
    ax.set_facecolor(SURFACE)
    if grid:
        ax.grid(True, axis=grid, color=GRID, lw=0.7)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def _title(ax, text):
    ax.set_title(text, loc="left", fontsize=11, fontweight="bold", color=INK, pad=6)


# ---- one function per panel; the separate figures and the overview both use them
def panel_trace(ax, D, y, name, unit, from_zero):
    t = D["t"]
    _axes(ax)
    ax.axvspan(t[D["cut"]], t[-1], color=WINDOW, lw=0, zorder=0)
    ax.plot(t, y, color=SERIES1, lw=0.6, alpha=0.3, zorder=2)
    ax.plot(t, running_mean(np.asarray(y, float), D["w"]), color=SERIES1, lw=2.0, zorder=3)
    if D["mark_start"]:
        ax.axvline(D["start"] / 1000.0, color=MUTED, lw=0.9, ls=(0, (4, 3)), zorder=4)
    ax.set_ylabel(unit, fontsize=9, labelpad=6)
    _title(ax, name)
    ax.ticklabel_format(axis="y", useOffset=False)
    if from_zero:
        ax.set_ylim(bottom=0)
    ax.set_xlim(t[0], t[-1])


def panel_clusters(ax, D, name="Largest NCAA cluster, and the largest off the protein"):
    t, n_ncaa = D["t"], D["n_ncaa"]
    _axes(ax)
    ax.axvspan(t[D["cut"]], t[-1], color=WINDOW, lw=0, zorder=0)
    ends = []
    for y, c, lab, end in ((D["size"], SERIES1, "largest cluster", "largest"),
                           (D["FREE"], SERIES2, "largest cluster off the protein", "off the protein")):
        ax.plot(t, y, color=c, lw=0.6, alpha=0.25, zorder=2)
        sm = running_mean(np.asarray(y, float), D["w"])
        ax.plot(t, sm, color=c, lw=2.0, zorder=3, label=lab)
        ends.append([sm[-1], end])
    top = max(n_ncaa / 2 + 2, D["size"].max() + 1)
    gap = 0.08 * top                                     # end labels closer than this are pushed apart
    if abs(ends[0][0] - ends[1][0]) < gap:
        mid = (ends[0][0] + ends[1][0]) / 2
        hi, lo = (0, 1) if ends[0][0] >= ends[1][0] else (1, 0)
        ends[hi][0], ends[lo][0] = mid + gap / 2, mid - gap / 2
    for yv, end in ends:
        ax.text(t[-1] + 0.01 * (t[-1] - t[0]), yv, end, color=INK2, fontsize=9, va="center", clip_on=False)
    if D["mark_start"]:
        ax.axvline(D["start"] / 1000.0, color=MUTED, lw=0.9, ls=(0, (4, 3)), zorder=4)
    ax.axhline(n_ncaa / 2, color=MUTED, lw=0.9, ls=(0, (4, 3)), zorder=1)
    ax.text(t[0] + 0.01 * (t[-1] - t[0]), n_ncaa / 2, "half the NCAAs", color=MUTED, fontsize=9,
            va="bottom")
    ax.set_xlim(t[0], t[-1])
    ax.set_ylim(0, top)
    ax.set_ylabel("NCAAs", fontsize=9, labelpad=6)
    _title(ax, name)
    ax.legend(loc="upper left", ncol=2, bbox_to_anchor=(0.18, 1.0), borderaxespad=0.2)


def panel_exchange(ax, D, label_size=None, only_bound=False):
    from matplotlib.colors import ListedColormap
    t, n_ncaa, order = D["t"], D["n_ncaa"], D["order"]
    ax.imshow(D["ADS"][:, order].T.astype(int), aspect="auto", interpolation="nearest",
              cmap=ListedColormap([SURFACE, SERIES1]), vmin=0, vmax=1,
              extent=[t[0], t[-1], n_ncaa - 0.5, -0.5])
    if D["start"] > D["T"][0] + 1e-6:
        ax.axvline(D["start"] / 1000.0, color=INK2, lw=0.9, ls=(0, (4, 3)))
    ax.set_yticks(range(n_ncaa))
    ax.set_yticklabels([f"NCAA{k + 1}" if not only_bound or D["frac"][k] >= 0.10 else "" for k in order],
                       fontsize=label_size or (7 if n_ncaa > 20 else 9))
    if only_bound:
        ax.set_ylabel("NCAAs, most-bound first\n(labelled if on >=10% of the window)", fontsize=8.5,
                      labelpad=6)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(t[0], t[-1])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    _title(ax, "Each NCAA on the protein (filled), most-bound first")


def panel_residence(ax, D):
    from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator, NullLocator
    comp, cens = D["comp"], D["cens"]
    if not len(comp) + len(cens):
        ax.set_axis_off()
        ax.text(0.5, 0.5, f"No visit lasted {ns(MIN_VISIT_PS)} ns or more\nin the window, so there is nothing to draw.",
                transform=ax.transAxes, ha="center", va="center", color=INK2, fontsize=10)
        return
    allv = np.r_[comp, cens]
    lo = MIN_VISIT_PS / 1000.0
    hi = max(allv.max() if len(allv) else 2 * lo, 2 * lo) * 1.05
    bins = np.logspace(np.log10(lo), np.log10(hi), 16)
    _axes(ax)
    ax.hist([comp, cens], bins=bins, stacked=True, color=[SERIES1, SERIES2], rwidth=0.88,
            label=["completed in the window", "cut off by the window edge (a lower bound)"], zorder=2)
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator([v for v in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500)
                                             if lo <= v <= hi]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("visit length  (ns)", fontsize=10, labelpad=6)
    ax.set_ylabel("visits", fontsize=9, labelpad=6)
    ax.legend(loc="lower left", ncol=1, bbox_to_anchor=(0, 1.0), borderaxespad=0.3)


def panel_contacts(ax, D):
    ll, lf = D["ll"], D["lf"]
    shown = [p for p in PARTS if p in lf]               # a part with no atoms has no row
    _axes(ax, grid="x")
    ys = np.arange(len(shown))[::-1]
    vals_all = [v for v in list(ll.values()) + list(lf.values()) if np.isfinite(v)]
    span = max([abs(v - 1) for v in vals_all] + [0.1])
    x0, x1 = max(-0.1, 1 - span - 0.15), 1 + span + 0.15     # a ratio is never below 0
    pad = 0.025 * (x1 - x0)                                  # gap between a dot and its value
    for d, (vals, c, lab) in zip((0.14, -0.14), ((lf, SERIES2, "NCAA–protein"), (ll, SERIES1, "NCAA–NCAA"))):
        x = np.array([vals[p] for p in shown])
        ax.hlines(ys + d, 1.0, x, color=c, lw=2.0, zorder=2)
        ax.scatter(x, ys + d, s=64, color=c, zorder=3, label=lab, edgecolors=SURFACE, linewidths=1.5)
        for xi, yi in zip(x, ys + d):
            if xi >= 1:
                at, ha = xi + pad, "left"
            elif xi - 5 * pad >= x0:
                at, ha = xi - pad, "right"
            else:                                            # no room left of the dot: past the line's end at 1
                at, ha = 1 + pad, "left"
            ax.text(at, yi, f"{xi:.2f}", color=INK2, fontsize=8.5, va="center", ha=ha, clip_on=False)
    ax.axvline(1.0, color=MUTED, lw=0.9, zorder=1)
    ax.set_xlim(x0, x1)
    ax.set_ylim(-0.6, len(shown) - 0.4)
    ax.set_yticks(ys)
    ax.set_yticklabels(shown, fontsize=10)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("contact share / atom share  (1 = random)", fontsize=10, labelpad=6)
    ax.legend(loc="lower left", ncol=2, bbox_to_anchor=(0, 1.0), borderaxespad=0.3)


# ---- captions, from the numbers
def _stand_settling(D):
    t, cut = D["t"], D["cut"]
    moved = [f"{lab} {'fell' if d < 0 else 'rose'} by {abs(d):.2g}{' Å' if 'Rg' in lab else ''}"
             for lab, ok, d, _ in D["verdicts"] if not ok]
    span = f"{t[cut]:.3g}-{t[-1]:.3g} ns"
    if moved:
        return (f"Over {span} (shaded), " + "; ".join(moved)
                + " -- more than their noise, so this stretch has not settled.")
    return f"Over {span} (shaded), all three stay level within their noise."


def _stand_clusters(D):
    t, cut, size, on = D["t"], D["cut"], D["size"], D["on"]
    last = D["T"] >= D["T"][-1] - BLOCK_PS + 1e-6
    off_share = 1 - (on[cut:] / np.maximum(size[cut:], 1)).mean()
    return (f"Over {t[cut]:.3g}-{t[-1]:.3g} ns the largest cluster averaged {size[cut:].mean():.1f} NCAAs, "
            f"{off_share:.0%} of its members off the protein; in the final {ns(BLOCK_PS)} ns it held "
            f"{size[last].mean():.0f} of {D['n_ncaa']}.")


def _stand_exchange(D):
    ex, comp = D["ex"], D["comp"]
    visited, stuck = len(ex["visited"]), len(ex["stuck"])
    rest = (" and the rest came and went" if visited > stuck else "")
    return (f"From {ns(D['start'])} ns, {visited} of {D['n_ncaa']} NCAAs stayed on the protein for "
            f"{ns(MIN_VISIT_PS)} ns or more at least once; {stuck} stayed on for more than {MOST:.0%} of the time{rest}"
            + (f" (median completed visit {np.median(comp):.3g} ns)." if len(comp) else "."))


def _stand_residence(D):
    comp, cens = D["comp"], D["cens"]
    return (f"{len(comp)} visits began and ended inside the window"
            + (f", median {np.median(comp):.3g} ns" if len(comp) else "")
            + f"; {len(cens)} were cut off by its start or end"
            + (f", the longest {cens.max():.3g} ns" if len(cens) else "")
            + ", so the slowest departures are not yet measured.")


def _stand_contacts(D):
    ll, lf = D["ll"], D["lf"]
    top_lf, top_ll = max(lf, key=lf.get), max(ll, key=ll.get)
    low = [p for p in PARTS if p in lf and ll[p] < 1 and lf[p] < 1]
    return (f"Most over-represented on the protein: {top_lf} ({lf[top_lf]:.2f}×); between NCAAs: {top_ll} "
            f"({ll[top_ll]:.2f}×)." + (f" Under-represented in both: {', '.join(low)}." if low else ""))


# ---- the figures
def fig_overview(plt, D, path):
    fig = plt.figure(figsize=(17, 14.5), facecolor=SURFACE)
    stand = " ".join([_stand_settling(D), _stand_clusters(D), _stand_exchange(D)])
    top = _header(fig, "Free-NCAA run: has it settled, and how do the NCAAs behave?", stand,
                  NOTE_CONTACT + " Thin: every frame analysed; thick: %.3g ns running mean. Panels e and f use "
                  "%s-%s ns." % (SMOOTH_NS, ns(D["start"]), ns(D["T"][-1])), left=0.05, wrap=9.5)
    left = fig.add_gridspec(4, 1, height_ratios=[1, 1.15, 1, 2.7], hspace=0.42, left=0.05, right=0.63,
                            top=top, bottom=0.075)
    right = fig.add_gridspec(2, 1, hspace=0.75, left=0.71, right=0.975, top=top - 0.05, bottom=0.075)
    a = fig.add_subplot(left[0])
    panel_trace(a, D, D["ADS"].sum(1), "a  NCAAs on the protein", "NCAAs", True)
    b = fig.add_subplot(left[1], sharex=a)
    panel_clusters(b, D, "b  Largest NCAA cluster, and the largest off the protein")
    c = fig.add_subplot(left[2], sharex=a)
    panel_trace(c, D, D["RG"], "c  Mean NCAA radius of gyration", "Rg  (Å)", False)
    d = fig.add_subplot(left[3], sharex=a)
    panel_exchange(d, D, label_size=7, only_bound=D["n_ncaa"] > 20)
    d.set_title("d  Each NCAA on the protein (filled), most-bound first", loc="left", fontsize=11,
                fontweight="bold", color=INK, pad=6)
    for ax in (a, b, c):
        ax.tick_params(labelbottom=False)
    d.set_xlabel("time  (ns)", fontsize=10, labelpad=6)
    e = fig.add_subplot(right[0])
    panel_residence(e, D)
    e.set_title("e  How long a visit to the protein lasts", loc="left", fontsize=11, fontweight="bold",
                color=INK, pad=34)
    f = fig.add_subplot(right[1])
    panel_contacts(f, D)
    f.set_title("f  Which part of the NCAA makes the contacts", loc="left", fontsize=11,
                fontweight="bold", color=INK, pad=22)
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def fig_settling(plt, D, path):
    fig = plt.figure(figsize=(11, 8.4), facecolor=SURFACE)
    top = _header(fig, "Has the free-NCAA run settled?", _stand_settling(D),
                  NOTE_CONTACT + " Thin: every frame analysed; thick: %g ns running mean." % SMOOTH_NS)
    gs = fig.add_gridspec(3, 1, hspace=0.45, left=0.075, right=0.985, top=top, bottom=0.125)
    panels = [(D["ADS"].sum(1), "NCAAs on the protein", "NCAAs", True),
              (D["size"], "Largest NCAA cluster", "NCAAs", True),
              (D["RG"], "Mean NCAA radius of gyration", "Rg  (Å)", False)]
    axes = []
    for i, (y, name, unit, z) in enumerate(panels):
        ax = fig.add_subplot(gs[i], sharex=axes[0] if axes else None)
        axes.append(ax)
        panel_trace(ax, D, y, name, unit, z)
        if i < 2:
            ax.tick_params(labelbottom=False)
    axes[-1].set_xlabel("time  (ns)", fontsize=10, labelpad=6)
    fig.savefig(path, dpi=170, facecolor=SURFACE)
    plt.close(fig)


def fig_single(plt, D, path, title, stand, note, size, draw, rect, xlabel=None):
    fig = plt.figure(figsize=size, facecolor=SURFACE)
    top = _header(fig, title, stand, note)
    ax = fig.add_axes([rect[0], rect[1], rect[2], top - rect[1] - rect[3]])
    draw(ax)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=10, labelpad=6)
    fig.savefig(path, dpi=170, facecolor=SURFACE)
    plt.close(fig)


def figures(T, ADS, size, on, FREE, RG, cut, start, n_ncaa, verdicts, rows, order, ex, contacts, dt, separate):
    plt = _plt()
    D = dict(T=T, t=T / 1000.0, ADS=ADS, FREE=FREE, RG=RG, cut=cut, start=start, n_ncaa=n_ncaa,
             verdicts=verdicts, ex=ex, size=size, on=on,
             w=max(1, int(round(SMOOTH_NS * 1000.0 / dt))), mark_start=abs(start - T[cut]) > 1e-6,
             order=order, frac={r[0]: r[1] for r in rows},
             comp=np.array([L for k, L, c in ex["visits"] if not c]),
             cens=np.array([L for k, L, c in ex["visits"] if c]),
             ll={r[1]: r[6] for r in contacts if r[0] == "NCAA-NCAA part"},
             lf={r[1]: r[6] for r in contacts if r[0] == "NCAA-protein part"})
    fig_overview(plt, D, "free_ncaa_overview.png")
    made = ["free_ncaa_overview.png"]
    if not separate:
        return made
    fig_settling(plt, D, "free_ncaa_settling.png")
    fig_single(plt, D, "free_ncaa_clusters.png", "Where do the NCAAs cluster?", _stand_clusters(D)
               + " Where the two lines meet, the largest cluster is floating in solution.",
               NOTE_CONTACT + " Thin: every frame; thick: %g ns running mean." % SMOOTH_NS, (11, 5.4),
               lambda ax: panel_clusters(ax, D, ""), (0.075, 0.17, 0.8, 0.0), "time  (ns)")
    fig_single(plt, D, "free_ncaa_exchange.png", "Do the NCAAs come and go, or stick?", _stand_exchange(D),
               "Filled: the NCAA has a heavy atom within 4.5 Å of the protein in that frame. Rows sorted by "
               "time on the protein in the window; labels are topology order.", (11, 8.6),
               lambda ax: panel_exchange(ax, D), (0.075, 0.12, 0.91, 0.0), "time  (ns)")
    if len(D["comp"]) + len(D["cens"]):
        fig_single(plt, D, "free_ncaa_residence.png", "How long does an NCAA stay on the protein?",
                   _stand_residence(D), "A visit lasts at least %g ns; gaps under %g ns are bridged. Bins are "
                   "equal widths on a log scale." % (MIN_VISIT_PS / 1000, MIN_VISIT_PS / 1000), (9.5, 5.4),
                   lambda ax: panel_residence(ax, D), (0.09, 0.18, 0.88, 0.07))
        made.append("free_ncaa_residence.png")
    fig_single(plt, D, "free_ncaa_contacts.png", "Which part of the NCAA makes the contacts?",
               _stand_contacts(D), "Ratio: a part's share of heavy-atom contacts (within 4.5 Å) divided by its "
               "share of the NCAA's heavy atoms.", (9.5, 5.0),
               lambda ax: panel_contacts(ax, D), (0.14, 0.19, 0.83, 0.07))
    made[1:1] = ["free_ncaa_settling.png", "free_ncaa_clusters.png", "free_ncaa_exchange.png"]
    made.append("free_ncaa_contacts.png")
    return made


def main():
    ts = np.loadtxt("free_ncaa_timeseries.csv", delimiter=",", skiprows=1)
    T, RG = np.round(ts[:, 0] * 1000.0, 3), ts[:, 2]         # ps, A
    size, on, FREE = ts[:, 3].astype(int), ts[:, 4].astype(int), ts[:, 5].astype(int)
    ADS = np.loadtxt("free_ncaa_adsorption.csv", delimiter=",", skiprows=1)[:, 1:].astype(bool)
    with open("free_ncaa_per_ncaa.csv") as fh:
        rows = [(int(r["NCAA"].removeprefix("NCAA")) - 1, float(r["frac_on_protein"]), int(r["visits"]), int(r["arrivals"]), int(r["departures"]),
                 float(r["longest_visit_ns"]), float(r["frac_in_largest_cluster"])) for r in csv.DictReader(fh)]
    with open("free_ncaa_contacts.csv") as fh:
        num = lambda v: float(v) if v not in ("", "nan") else float("nan")      # some rows leave a column empty
        contacts = [(r["kind"], r["part_a"], r["part_b"], r["atoms"], num(r["share_of_contacts"]),
                     num(r["share_expected"]), num(r["ratio"])) for r in csv.DictReader(fh)]
    separate = meta["separate"] or "--separate" in sys.argv[1:]
    made = figures(T, ADS, size, on, FREE, RG, meta["cut"], meta["start"], meta["n_ncaa"], [tuple(v) for v in meta["verdicts"]],
                   rows, meta["order"], meta["ex"], contacts, meta["dt"], separate)
    print("wrote " + ", ".join(made))


if __name__ == "__main__":
    main()
