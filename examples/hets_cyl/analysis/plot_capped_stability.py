#!/usr/bin/env python3
"""
Draw stability.png from the files check_capped_stability_3.py writes:

    stability.csv         every trace, frame by frame (junction gaps, hydrogen bonds)
    stability_rmsd.csv    the core RMSD
    stability_meta.json   which junction is the cap, which the free end, the averaging window

Usage:
    python3 plot_capped_stability.py            # or: python3 plot_capped_stability.py stability.png

This script only draws. Change it freely: the numbers and the verdicts come from the check script.
"""

import json
import sys
import textwrap
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out = sys.argv[1] if len(sys.argv) > 1 else "stability.png"
meta = json.load(open("stability_meta.json"))
HB_OK, pair_names, core_chains = meta["HB_OK"], meta["pair_names"], meta["core_chains"]
CAP_PAIR, FREE_END_PAIR, EQ_FRACTION = meta["CAP_PAIR"], meta["FREE_END_PAIR"], meta["EQ_FRACTION"]

table = np.loadtxt("stability.csv", delimiter=",", skiprows=1)
n = len(pair_names)
times = np.round(table[:, 0] * 1000.0, 3)       # ps
spans = table[:, 2:2 + n]
hbonds = table[:, 2 + n:2 + 2 * n] if HB_OK else None
r = np.loadtxt("stability_rmsd.csv", delimiter=",", skiprows=1)
t_rmsd, rmsd = np.round(r[:, 0] * 1000.0, 3), r[:, 1]
cut = int(len(times) * (1 - EQ_FRACTION))
cut_r = int(len(rmsd) * (1 - EQ_FRACTION))

SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#1a1a19", "#52514e", "#8a8984", "#e8e7e3"
CAP_C, END_C, BAND_C, CORE_C = "#eb6834", "#2a78d6", "#9a9a94", "#1baf7a"
cap_nm, end_nm = pair_names[CAP_PAIR], pair_names[FREE_END_PAIR]

plt.rcParams.update({"font.size": 11, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "axes.linewidth": 0.8, "xtick.color": MUTED, "ytick.color": MUTED,
                     "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
                     "xtick.major.size": 3, "ytick.major.size": 3, "text.color": INK})

n_lower = 2 if HB_OK else 1
wide = n_lower == 2
FIG_W, FIG_H = 7.6 * n_lower, 9.0 if wide else 10.4
LEFT, RIGHT = 0.52 / FIG_W + 0.02, 0.985
fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor=SURFACE)

# --- header first, so the panels can start below however tall it turns out
avail_in = (RIGHT - LEFT) * FIG_W                 # text width, inches
line_h = 0.165 / FIG_H                            # one text line, figure fractions
title = "The capped fibril: is the engineered end holding like the fibril's own?"
sub = (f"Every junction between neighbouring chains, over the whole production run. "
       f"The strip beside each lower panel gives that junction's average over the "
       f"last {int(EQ_FRACTION * 100)}% of the run, in order along the rod, so the two "
       f"ends can be read against the middle.")
y = 0.985
for line in textwrap.wrap(title, width=max(20, int(avail_in / 0.1575))):
    fig.text(LEFT, y, line, fontsize=15, fontweight="bold", color=INK, va="top")
    y -= line_h * 1.45
y -= line_h * 0.3
for line in textwrap.wrap(sub, width=max(28, int(avail_in / 0.0920))):
    fig.text(LEFT, y, line, fontsize=10, color=INK2, va="top")
    y -= line_h
key = [(CORE_C, 2.0, "the fibril core, chains "
        + ", ".join(core_chains) + "   (panel A)", True),
       (CAP_C, 2.4, f"cap junction  {cap_nm}   (panels B and C)", True),
       (END_C, 2.4, f"the fibril's other free end  {end_nm}", True),
       (BAND_C, 1.2, f"the other {len(pair_names) - 2} junctions", False)]
handles = [plt.Line2D([], [], color=c, lw=lw, label=lab) for c, lw, lab, _ in key]
y -= line_h * 0.5
leg = fig.legend(handles=handles, loc="upper left", frameon=False, fontsize=10,
                 ncol=2 if wide else 1, handlelength=1.8, handletextpad=0.7,
                 columnspacing=3.0, borderaxespad=0, bbox_to_anchor=(LEFT, y))
for t, (c, _, _, bold) in zip(leg.get_texts(), key):
    t.set_color(c if bold else INK2)
    if bold:
        t.set_fontweight("bold")
y -= line_h * (4.6 if wide else 7.0)

gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.3], hspace=0.42,
                      left=LEFT, right=RIGHT, top=y, bottom=0.115 if wide else 0.135)
axA = fig.add_subplot(gs[0])
gsl = gs[1].subgridspec(1, n_lower, wspace=0.26)
lower = []
for k in range(n_lower):
    g = gsl[k].subgridspec(1, 2, width_ratios=[1, 0.34], wspace=0.07)
    ax = fig.add_subplot(g[0])
    lower.append((ax, fig.add_subplot(g[1], sharey=ax)))

for ax in [axA] + [a for pair in lower for a in pair]:
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, lw=0.7)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

# ---- A: interior RMSD
axA.set_title("A   Has the fibril core stopped changing shape?", loc="left", fontsize=12,
              fontweight="bold", color=INK, pad=10)
axA.plot(t_rmsd / 1000, rmsd, lw=1.6, color=CORE_C)
axA.axvline(t_rmsd[cut_r] / 1000, color=MUTED, lw=0.9, ls=(0, (4, 3)))
axA.text(t_rmsd[cut_r] / 1000 + 0.012 * t_rmsd[-1] / 1000, axA.get_ylim()[1],
         "averages taken over this part", fontsize=9.5, color=MUTED, va="top")
axA.set_xlabel("time  (ns)", fontsize=10, labelpad=6)
axA.set_ylabel("core Cα RMSD from the start  (Å)", fontsize=10, labelpad=6)
axA.set_xlim(t_rmsd[0] / 1000, t_rmsd[-1] / 1000)

# ---- B (and C): every junction, plus its average along the rod
panels = [(spans, "B   Gap between neighbouring chains",
           "distance between chain centres  (Å)")]
if HB_OK:
    panels.append((hbonds, "C   Hydrogen bonds at each junction",
                   "backbone hydrogen bonds"))

for (ax, st), (data, title, ylab) in zip(lower, panels):
    ax.set_title(title, loc="left", fontsize=12, fontweight="bold", color=INK, pad=10)
    for j in range(len(pair_names)):
        if j not in (CAP_PAIR, FREE_END_PAIR):
            ax.plot(times / 1000, data[:, j], lw=0.7, color=BAND_C, alpha=0.45)
    ax.plot(times / 1000, data[:, FREE_END_PAIR], lw=1.8, color=END_C, zorder=3)
    ax.plot(times / 1000, data[:, CAP_PAIR], lw=1.8, color=CAP_C, zorder=3)
    ax.axvline(times[cut] / 1000, color=MUTED, lw=0.9, ls=(0, (4, 3)))
    ax.set_xlabel("time  (ns)", fontsize=10, labelpad=6)
    ax.set_ylabel(ylab, fontsize=10, labelpad=6)
    ax.set_xlim(times[0] / 1000, times[-1] / 1000)

    m, sd = data[cut:].mean(0), data[cut:].std(0)
    xs = np.arange(len(pair_names))
    st.errorbar(xs, m, yerr=sd, fmt="none", ecolor=MUTED, elinewidth=0.9, capsize=0)
    st.scatter(xs, m, zorder=3, linewidths=0,
               s=[46 if j in (CAP_PAIR, FREE_END_PAIR) else 20 for j in xs],
               c=[CAP_C if j == CAP_PAIR else END_C if j == FREE_END_PAIR else BAND_C
                  for j in xs])
    st.set_xlim(-1.4, len(pair_names) + 0.4)
    st.set_xticks([0, len(pair_names) - 1])
    st.set_xticklabels([pair_names[0], pair_names[-1]], fontsize=9)
    st.tick_params(labelleft=False)
    st.set_xlabel("average per junction,\nalong the rod", fontsize=9.5, labelpad=6)

note = ("Error bars are the spread over the averaging window, not an uncertainty on "
        "the mean: each junction is one measurement in one simulation.")
yn = 0.016
for line in reversed(textwrap.wrap(note, width=max(30, int(avail_in / 0.083)))):
    fig.text(LEFT, yn, line, fontsize=9, color=MUTED, va="bottom")
    yn += line_h
fig.savefig(out, dpi=170, facecolor=SURFACE)
print(f"wrote {out}")
