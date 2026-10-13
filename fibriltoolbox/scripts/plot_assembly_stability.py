#!/usr/bin/env python3
"""
Draw assembly_stability.png from the files check_assembly_stability.py writes:

    assembly_stability.csv         one row per frame: C-alpha RMSD, radius of gyration, contacts kept per interface
    assembly_stability_meta.json   the chains and interfaces found, the window and the verdicts

Usage:
    python3 plot_assembly_stability.py          # or: python3 plot_assembly_stability.py out.png

This script only draws. Change it freely: the numbers and the verdicts come from check_assembly_stability.py.
"""
import json
import sys
import textwrap
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out = sys.argv[1] if len(sys.argv) > 1 else "assembly_stability.png"
meta = json.load(open("assembly_stability_meta.json"))
d = np.loadtxt("assembly_stability.csv", delimiter=",", skiprows=1)
t, rmsd, rg, Q = d[:, 0], d[:, 1], d[:, 2], d[:, 3:]
cut = int(len(t) * (1 - meta["EQ_FRACTION"]))
names = meta["interfaces"]

SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#1a1a19", "#52514e", "#8a8984", "#e8e7e3"
BLUE, ORANGE, BAND = "#2a78d6", "#eb6834", "#9a9a94"
plt.rcParams.update({"font.size": 11, "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "axes.linewidth": 0.8,
                     "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
                     "xtick.major.size": 3, "ytick.major.size": 3, "text.color": INK})
n_panels = 3 if len(names) else 2
fig = plt.figure(figsize=(9.6, 3.4 * n_panels + 1.6), facecolor=SURFACE)
failed = [v[0] for v in meta["verdicts"] if not v[1]]
title = "Did the assembly hold together during the run?"
stand = (f"{len(meta['chains'])} chains" + (f", {len(names)} interfaces between them" if names else "") + ". "
         + ("Every check passes over the shaded part of the run." if not failed
            else "Over the shaded part of the run these checks did not pass: " + "; ".join(failed) + "."))
LEFT, RIGHT = 0.12, 0.975
line_h = 0.165 / fig.get_size_inches()[1]
y = 0.985
fig.text(LEFT, y, title, fontsize=15, fontweight="bold", color=INK, va="top")
y -= line_h * 1.9
for line in textwrap.wrap(stand, width=100):
    fig.text(LEFT, y, line, fontsize=10, color=INK2, va="top")
    y -= line_h
gs = fig.add_gridspec(n_panels, 1, hspace=0.42, left=LEFT, right=RIGHT, top=y - line_h * 3.2, bottom=0.09)
axes = [fig.add_subplot(gs[k]) for k in range(n_panels)]
for ax in axes:
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, lw=0.7)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.axvspan(t[cut], t[-1], color="#efeee9", lw=0, zorder=0)
    ax.set_xlim(t[0], t[-1])
axes[0].plot(t, rmsd, lw=1.6, color=BLUE)
axes[0].set_title("A   Has the assembly stopped changing shape?", loc="left", fontsize=12, fontweight="bold", color=INK, pad=10)
axes[0].set_ylabel("Cα RMSD from the start  (Å)", fontsize=10, labelpad=6)
axes[1].plot(t, rg, lw=1.6, color=BLUE)
axes[1].set_title("B   Has its size stopped changing?", loc="left", fontsize=12, fontweight="bold", color=INK, pad=10)
axes[1].set_ylabel("radius of gyration  (Å)", fontsize=10, labelpad=6)
if len(names):
    ax = axes[2]
    worst = int(np.argmin(Q[cut:].mean(0)))
    for k in range(len(names)):
        if k != worst:
            ax.plot(t, Q[:, k], lw=0.8, color=BAND, alpha=0.6)
    ax.plot(t, Q[:, worst], lw=1.8, color=ORANGE, label=f"{names[worst]}, the interface that kept the fewest")
    ax.axhline(meta["Q_MIN"], color=MUTED, lw=0.9, ls=(0, (4, 3)))
    ax.text(t[0] + 0.01 * (t[-1] - t[0]), meta["Q_MIN"], f"limit {meta['Q_MIN']:.0%}", fontsize=9, color=MUTED, va="bottom")
    ax.set_ylim(min(0.0, Q.min()), 1.03)
    ax.set_title("C   Are the chains still in contact with their neighbours?", loc="left", fontsize=12, fontweight="bold", color=INK, pad=10)
    ax.set_ylabel("starting Cα contacts kept", fontsize=10, labelpad=6)
    leg = ax.legend(loc="lower left", frameon=False, fontsize=9.5)
    for txt in leg.get_texts():
        txt.set_color(ORANGE)
axes[-1].set_xlabel("time  (ns)", fontsize=10, labelpad=6)
fig.text(LEFT, 0.012, f"Shaded: the last {meta['EQ_FRACTION']:.0%} of the run, which the checks judge. A contact is a pair of Cα atoms of "
         f"two chains within {meta['CONTACT_CA']:g} Å in the first frame.", fontsize=9, color=MUTED, va="bottom", wrap=True)
fig.savefig(out, dpi=170, facecolor=SURFACE)
print(f"wrote {out}")
