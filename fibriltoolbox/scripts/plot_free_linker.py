#!/usr/bin/env python3
"""
Draw the free-linker figures from the files free_linker_analysis.py writes:

    free_linker_timeseries.csv   one row per frame: adsorbed, mean Rg, cluster sizes
    free_linker_adsorption.csv   one row per frame, one column per linker: 1 = on the fibril
    free_linker_per_linker.csv   one row per linker
    free_linker_contacts.csv     part-by-part contact shares and ratios
    free_linker_meta.json        the window, the verdicts, the visits and the settings the analysis used

Usage:
    python3 plot_free_linker.py                 # free_linker_overview.png
    python3 plot_free_linker.py --separate      # also each graph as its own PNG

This script only draws. Change it freely: the numbers and the report come from free_linker_analysis.py.
"""

import csv
import json
import sys
import textwrap
import numpy as np

meta = json.load(open("free_linker_meta.json"))
SMOOTH_NS, MIN_VISIT_PS, MOST, BLOCK_PS, PARTS = (meta[k] for k in ("SMOOTH_NS", "MIN_VISIT_PS", "MOST", "BLOCK_PS", "PARTS"))
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#1a1a19", "#52514e", "#8a8984", "#e8e7e3"
WINDOW = "#efeee9"             # shading for the window the drift test judges


def running_mean(y, w):
    if w <= 1:
        return y
    pad = np.r_[np.full(w // 2, y[0]), y, np.full(w - 1 - w // 2, y[-1])]
    return np.convolve(pad, np.ones(w) / w, mode="valid")


# ----------------------------------------------------------------- figures
SERIES1, SERIES2 = "#2a78d6", "#eb6834"      # categorical slots 1 and 2 (blue, orange)
NOTE_CONTACT = "On the fibril: any linker heavy atom within 4.5 Å of a fibril heavy atom. Cluster: linkers " \
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


def panel_clusters(ax, D, name="Largest linker cluster, and the largest off the fibril"):
    t, n_link = D["t"], D["n_link"]
    _axes(ax)
    ax.axvspan(t[D["cut"]], t[-1], color=WINDOW, lw=0, zorder=0)
    for y, c, lab, end in ((D["size"], SERIES1, "largest cluster", "largest"),
                           (D["FREE"], SERIES2, "largest cluster off the fibril", "off the fibril")):
        ax.plot(t, y, color=c, lw=0.6, alpha=0.25, zorder=2)
        sm = running_mean(np.asarray(y, float), D["w"])
        ax.plot(t, sm, color=c, lw=2.0, zorder=3, label=lab)
        ax.text(t[-1] + 0.01 * (t[-1] - t[0]), sm[-1], end, color=INK2, fontsize=9, va="center",
                clip_on=False)
    if D["mark_start"]:
        ax.axvline(D["start"] / 1000.0, color=MUTED, lw=0.9, ls=(0, (4, 3)), zorder=4)
    ax.axhline(n_link / 2, color=MUTED, lw=0.9, ls=(0, (4, 3)), zorder=1)
    ax.text(t[0] + 0.01 * (t[-1] - t[0]), n_link / 2, "half the linkers", color=MUTED, fontsize=9,
            va="bottom")
    ax.set_xlim(t[0], t[-1])
    ax.set_ylim(0, max(n_link / 2 + 2, D["size"].max() + 1))
    ax.set_ylabel("linkers", fontsize=9, labelpad=6)
    _title(ax, name)
    ax.legend(loc="upper left", ncol=2, bbox_to_anchor=(0.18, 1.0), borderaxespad=0.2)


def panel_exchange(ax, D, label_size=None, only_bound=False):
    from matplotlib.colors import ListedColormap
    t, n_link, order = D["t"], D["n_link"], D["order"]
    ax.imshow(D["ADS"][:, order].T.astype(int), aspect="auto", interpolation="nearest",
              cmap=ListedColormap([SURFACE, SERIES1]), vmin=0, vmax=1,
              extent=[t[0], t[-1], n_link - 0.5, -0.5])
    if D["start"] > D["T"][0] + 1e-6:
        ax.axvline(D["start"] / 1000.0, color=INK2, lw=0.9, ls=(0, (4, 3)))
    ax.set_yticks(range(n_link))
    ax.set_yticklabels([f"L{k + 1}" if not only_bound or D["frac"][k] >= 0.10 else "" for k in order],
                       fontsize=label_size or (7 if n_link > 20 else 9))
    if only_bound:
        ax.set_ylabel("linkers, most-bound first\n(labelled if on >=10% of the window)", fontsize=8.5,
                      labelpad=6)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(t[0], t[-1])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    _title(ax, "Each linker on the fibril (filled), most-bound first")


def panel_residence(ax, D):
    from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator, NullLocator
    comp, cens = D["comp"], D["cens"]
    allv = np.r_[comp, cens]
    lo = MIN_VISIT_PS / 1000.0
    hi = max(allv.max() if len(allv) else 2 * lo, 2 * lo) * 1.05
    bins = np.logspace(np.log10(lo), np.log10(hi), 16)
    _axes(ax)
    ax.hist([comp, cens], bins=bins, stacked=True, color=[SERIES1, SERIES2], rwidth=0.88,
            label=["completed in the window", "cut off by the window edge (a lower bound)"], zorder=2)
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator([v for v in (1, 2, 5, 10, 20, 50, 100, 200, 500) if lo <= v <= hi]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("visit length  (ns)", fontsize=10, labelpad=6)
    ax.set_ylabel("visits", fontsize=9, labelpad=6)
    ax.legend(loc="lower left", ncol=1, bbox_to_anchor=(0, 1.0), borderaxespad=0.3)


def panel_contacts(ax, D):
    ll, lf = D["ll"], D["lf"]
    _axes(ax, grid="x")
    ys = np.arange(len(PARTS))[::-1]
    for d, (vals, c, lab) in zip((0.14, -0.14), ((lf, SERIES2, "linker–fibril"), (ll, SERIES1, "linker–linker"))):
        x = np.array([vals[p] for p in PARTS])
        ax.hlines(ys + d, 1.0, x, color=c, lw=2.0, zorder=2)
        ax.scatter(x, ys + d, s=64, color=c, zorder=3, label=lab, edgecolors=SURFACE, linewidths=1.5)
        for xi, yi in zip(x, ys + d):
            ax.text(xi + (0.02 if xi >= 1 else -0.02), yi, f"{xi:.2f}", color=INK2, fontsize=8.5,
                    va="center", ha="left" if xi >= 1 else "right")
    ax.axvline(1.0, color=MUTED, lw=0.9, zorder=1)
    span = max(abs(v - 1) for v in list(ll.values()) + list(lf.values()))
    ax.set_xlim(1 - span - 0.15, 1 + span + 0.15)
    ax.set_ylim(-0.6, len(PARTS) - 0.4)
    ax.set_yticks(ys)
    ax.set_yticklabels(PARTS, fontsize=10)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("contact share / atom share  (1 = random)", fontsize=10, labelpad=6)
    ax.legend(loc="lower left", ncol=2, bbox_to_anchor=(0, 1.0), borderaxespad=0.3)


# ---- captions, from the numbers
def _stand_settling(D):
    t, cut = D["t"], D["cut"]
    moved = [f"{lab} {'fell' if d < 0 else 'rose'} by {abs(d):.1f}{' Å' if 'Rg' in lab else ''}"
             for lab, ok, d, _ in D["verdicts"] if not ok]
    if moved:
        return (f"Over {t[cut]:.0f}-{t[-1]:.0f} ns (shaded), " + "; ".join(moved)
                + " -- more than their noise, so this stretch has not settled.")
    return f"Over {t[cut]:.0f}-{t[-1]:.0f} ns (shaded), all three stay level within their noise."


def _stand_clusters(D):
    t, cut, size, on = D["t"], D["cut"], D["size"], D["on"]
    last = D["T"] >= D["T"][-1] - BLOCK_PS + 1e-6
    off_share = 1 - (on[cut:] / np.maximum(size[cut:], 1)).mean()
    return (f"Over {t[cut]:.0f}-{t[-1]:.0f} ns the largest cluster averaged {size[cut:].mean():.1f} linkers, "
            f"{off_share:.0%} of its members off the fibril; in the final {BLOCK_PS / 1000:.0f} ns it held "
            f"{size[last].mean():.0f} of {D['n_link']}.")


def _stand_exchange(D):
    ex, comp = D["ex"], D["comp"]
    return (f"From {D['start'] / 1000:.0f} ns, {len(ex['visited'])} of {D['n_link']} linkers visited the fibril; "
            f"{len(ex['stuck'])} stayed on for more than {MOST:.0%} of the time and the rest came and went"
            + (f" (median completed visit {np.median(comp):.1f} ns)." if len(comp) else "."))


def _stand_residence(D):
    comp, cens = D["comp"], D["cens"]
    return (f"{len(comp)} visits began and ended inside the window"
            + (f", median {np.median(comp):.1f} ns" if len(comp) else "")
            + f"; {len(cens)} were cut off by its start or end"
            + (f", the longest {cens.max():.0f} ns" if len(cens) else "")
            + ", so the slowest departures are not yet measured.")


def _stand_contacts(D):
    ll, lf = D["ll"], D["lf"]
    top_lf, top_ll = max(lf, key=lf.get), max(ll, key=ll.get)
    low = [p for p in PARTS if ll[p] < 1 and lf[p] < 1]
    return (f"Most over-represented on the fibril: {top_lf} ({lf[top_lf]:.2f}×); between linkers: {top_ll} "
            f"({ll[top_ll]:.2f}×)." + (f" Under-represented in both: {', '.join(low)}." if low else ""))


# ---- the figures
def fig_overview(plt, D, path):
    fig = plt.figure(figsize=(17, 14.5), facecolor=SURFACE)
    stand = " ".join([_stand_settling(D), _stand_clusters(D), _stand_exchange(D)])
    top = _header(fig, "Free-linker run: has it settled, and how do the linkers behave?", stand,
                  NOTE_CONTACT + " Thin: every frame analysed; thick: %g ns running mean. Panels e and f use "
                  "%.0f-%.0f ns." % (SMOOTH_NS, D["start"] / 1000, D["T"][-1] / 1000), left=0.05, wrap=9.5)
    left = fig.add_gridspec(4, 1, height_ratios=[1, 1.15, 1, 2.7], hspace=0.42, left=0.05, right=0.63,
                            top=top, bottom=0.075)
    right = fig.add_gridspec(2, 1, hspace=0.75, left=0.71, right=0.975, top=top - 0.05, bottom=0.075)
    a = fig.add_subplot(left[0])
    panel_trace(a, D, D["ADS"].sum(1), "a  Linkers on the fibril", "linkers", True)
    b = fig.add_subplot(left[1], sharex=a)
    panel_clusters(b, D, "b  Largest linker cluster, and the largest off the fibril")
    c = fig.add_subplot(left[2], sharex=a)
    panel_trace(c, D, D["RG"], "c  Mean linker radius of gyration", "Rg  (Å)", False)
    d = fig.add_subplot(left[3], sharex=a)
    panel_exchange(d, D, label_size=7, only_bound=D["n_link"] > 20)
    d.set_title("d  Each linker on the fibril (filled), most-bound first", loc="left", fontsize=11,
                fontweight="bold", color=INK, pad=6)
    for ax in (a, b, c):
        ax.tick_params(labelbottom=False)
    d.set_xlabel("time  (ns)", fontsize=10, labelpad=6)
    e = fig.add_subplot(right[0])
    if len(D["comp"]) + len(D["cens"]):
        panel_residence(e, D)
    e.set_title("e  How long a visit to the fibril lasts", loc="left", fontsize=11, fontweight="bold",
                color=INK, pad=34)
    f = fig.add_subplot(right[1])
    panel_contacts(f, D)
    f.set_title("f  Which part of the linker makes the contacts", loc="left", fontsize=11,
                fontweight="bold", color=INK, pad=22)
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def fig_settling(plt, D, path):
    fig = plt.figure(figsize=(11, 8.4), facecolor=SURFACE)
    top = _header(fig, "Has the free-linker run settled?", _stand_settling(D),
                  NOTE_CONTACT + " Thin: every frame analysed; thick: %g ns running mean." % SMOOTH_NS)
    gs = fig.add_gridspec(3, 1, hspace=0.45, left=0.075, right=0.985, top=top, bottom=0.125)
    panels = [(D["ADS"].sum(1), "Linkers on the fibril", "linkers", True),
              (D["size"], "Largest linker cluster", "linkers", True),
              (D["RG"], "Mean linker radius of gyration", "Rg  (Å)", False)]
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


def figures(T, ADS, size, on, FREE, RG, cut, start, n_link, verdicts, rows, order, ex, contacts, dt, separate):
    plt = _plt()
    D = dict(T=T, t=T / 1000.0, ADS=ADS, FREE=FREE, RG=RG, cut=cut, start=start, n_link=n_link,
             verdicts=verdicts, ex=ex, size=size, on=on,
             w=max(1, int(round(SMOOTH_NS * 1000.0 / dt))), mark_start=abs(start - T[cut]) > 1e-6,
             order=order, frac={r[0]: r[1] for r in rows},
             comp=np.array([L for k, L, c in ex["visits"] if not c]),
             cens=np.array([L for k, L, c in ex["visits"] if c]),
             ll={r[1]: r[6] for r in contacts if r[0] == "linker-linker part"},
             lf={r[1]: r[6] for r in contacts if r[0] == "linker-fibril part"})
    fig_overview(plt, D, "free_linker_overview.png")
    made = ["free_linker_overview.png"]
    if not separate:
        return made
    fig_settling(plt, D, "free_linker_settling.png")
    fig_single(plt, D, "free_linker_clusters.png", "Where do the linkers cluster?", _stand_clusters(D)
               + " Where the two lines meet, the largest cluster is floating in solution.",
               NOTE_CONTACT + " Thin: every frame; thick: %g ns running mean." % SMOOTH_NS, (11, 5.4),
               lambda ax: panel_clusters(ax, D, ""), (0.075, 0.17, 0.8, 0.0), "time  (ns)")
    fig_single(plt, D, "free_linker_exchange.png", "Do the linkers come and go, or stick?", _stand_exchange(D),
               "Filled: the linker has a heavy atom within 4.5 Å of the fibril in that frame. Rows sorted by "
               "time on the fibril in the window; labels are topology order.", (11, 8.6),
               lambda ax: panel_exchange(ax, D), (0.075, 0.12, 0.91, 0.0), "time  (ns)")
    if len(D["comp"]) + len(D["cens"]):
        fig_single(plt, D, "free_linker_residence.png", "How long does a linker stay on the fibril?",
                   _stand_residence(D), "A visit lasts at least %g ns; gaps under %g ns are bridged. Bins are "
                   "equal widths on a log scale." % (MIN_VISIT_PS / 1000, MIN_VISIT_PS / 1000), (9.5, 5.4),
                   lambda ax: panel_residence(ax, D), (0.09, 0.18, 0.88, 0.07))
        made.append("free_linker_residence.png")
    fig_single(plt, D, "free_linker_contacts.png", "Which part of the linker makes the contacts?",
               _stand_contacts(D), "Ratio: a part's share of heavy-atom contacts (within 4.5 Å) divided by its "
               "share of the linker's heavy atoms.", (9.5, 5.0),
               lambda ax: panel_contacts(ax, D), (0.14, 0.19, 0.83, 0.07))
    made[1:1] = ["free_linker_settling.png", "free_linker_clusters.png", "free_linker_exchange.png"]
    made.append("free_linker_contacts.png")
    return made


def main():
    ts = np.loadtxt("free_linker_timeseries.csv", delimiter=",", skiprows=1)
    T, RG = np.round(ts[:, 0] * 1000.0, 3), ts[:, 2]         # ps, A
    size, on, FREE = ts[:, 3].astype(int), ts[:, 4].astype(int), ts[:, 5].astype(int)
    ADS = np.loadtxt("free_linker_adsorption.csv", delimiter=",", skiprows=1)[:, 1:].astype(bool)
    with open("free_linker_per_linker.csv") as fh:
        rows = [(int(r["linker"][1:]) - 1, float(r["frac_on_fibril"]), int(r["visits"]), int(r["arrivals"]), int(r["departures"]),
                 float(r["longest_visit_ns"]), float(r["frac_in_largest_cluster"])) for r in csv.DictReader(fh)]
    with open("free_linker_contacts.csv") as fh:
        num = lambda v: float(v) if v not in ("", "nan") else float("nan")      # some rows leave a column empty
        contacts = [(r["kind"], r["part_a"], r["part_b"], r["atoms"], num(r["share_of_contacts"]),
                     num(r["share_expected"]), num(r["ratio"])) for r in csv.DictReader(fh)]
    separate = meta["separate"] or "--separate" in sys.argv[1:]
    made = figures(T, ADS, size, on, FREE, RG, meta["cut"], meta["start"], meta["n_link"], [tuple(v) for v in meta["verdicts"]],
                   rows, meta["order"], meta["ex"], contacts, meta["dt"], separate)
    print("wrote " + ", ".join(made))


if __name__ == "__main__":
    main()
