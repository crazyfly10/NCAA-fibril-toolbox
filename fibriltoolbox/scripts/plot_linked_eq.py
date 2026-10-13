#!/usr/bin/env python3
"""
Draw linked_eq.png from the files linked_eq_check.py writes:

    linked_eq.csv         one row per frame; for each linker its Rg, reach, atoms on the fibril, biotin on
    linked_eq_meta.json   the linker labels and the settings the analysis used

Usage:
    python3 plot_linked_eq.py                   # or: python3 plot_linked_eq.py linked_eq.png

This script only draws. Change it freely: the numbers and the report come from linked_eq_check.py.
"""

import json
import sys
import numpy as np

meta = json.load(open("linked_eq_meta.json"))
CONTACT, SMOOTH_NS, TAIL_NS = meta["CONTACT"], meta["SMOOTH_NS"], meta["TAIL_NS"]
SITE_NAMES, SITE_SHORT = meta["SITE_NAMES"], meta["SITE_SHORT"]


def runs(s):
    d = np.diff(np.r_[0, np.asarray(s, int), 0])
    return list(zip(np.where(d == 1)[0], np.where(d == -1)[0]))


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


def plot(t, data, labels, n_distal, mark_ns=None, path="linked_eq.png",
         note=None, bound_cut=0.5):
    """t in ps; data[k] = (frames, 4): Rg, reach, on, biotin on (0/1, or a
    fraction when the input is already averaged)"""
    import textwrap
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#1a1a19", "#52514e", "#8a8984", "#e8e7e3"
    COLORS = ["#2a78d6", "#eb6834"]              # categorical slots 1 and 2
    t_ns = np.asarray(t, float) / 1000.0
    dt = float(np.median(np.diff(t))) if len(t) > 1 else 100.0
    w = max(1, int(round(SMOOTH_NS * 1000 / dt)))
    g = max(1, int(round(1000 / dt)))
    names = [SITE_NAMES.get(lab, lab) for lab in labels]
    short = [SITE_SHORT.get(lab, lab) for lab in labels]

    plt.rcParams.update({"font.size": 11, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                         "axes.linewidth": 0.8, "xtick.color": MUTED, "ytick.color": MUTED,
                         "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
                         "xtick.major.size": 3, "ytick.major.size": 3, "text.color": INK})
    FIG_W, FIG_H = 9.6, 11.2
    LEFT, RIGHT = 0.13, 0.975
    fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor=SURFACE)
    avail_in = (RIGHT - LEFT) * FIG_W
    line_h = 0.165 / FIG_H

    # ---- header: title, a standfirst written from the data, the key
    since = [bound_since(t, data[k][:, 3], g, bound_cut) for k in range(len(labels))]
    tail = t_ns >= t_ns[-1] - TAIL_NS
    parts = []
    for nm, b in zip(short, since):
        parts.append(f"the {nm} linker's from {b / 1000:.0f} ns" if b is not None
                     else f"the {nm} linker's not by the end")
    stand = ("Both linkers start out in the water and end up lying on the fibril. Each biotin "
             "sticks and stays: " + " and ".join(parts) + ". The PEG keeps spreading over the "
             f"surface after that: over the last {TAIL_NS:g} ns, "
             + " and ".join(f"{data[k][tail, 2].mean():.0f}" for k in range(len(labels)))
             + f" of each linker's {n_distal} PEG and biotin atoms touch the fibril.")
    title = "The two linkers lie down on the fibril during the production run"
    y = 0.985
    for line in textwrap.wrap(title, width=max(20, int(avail_in / 0.118))):
        fig.text(LEFT, y, line, fontsize=15, fontweight="bold", color=INK, va="top")
        y -= line_h * 1.45
    y -= line_h * 0.3
    for line in textwrap.wrap(stand, width=max(28, int(avail_in / 0.078))):
        fig.text(LEFT, y, line, fontsize=10, color=INK2, va="top")
        y -= line_h
    y -= line_h * 0.6
    handles = [plt.Line2D([], [], color=c, lw=2.4, label=f"{nm}   ({lab})")
               for c, nm, lab in zip(COLORS, names, labels)]
    leg = fig.legend(handles=handles, loc="upper left", frameon=False, fontsize=10,
                     ncol=2, handlelength=1.8, handletextpad=0.7, columnspacing=3.0,
                     borderaxespad=0, bbox_to_anchor=(LEFT, y))
    for txt, c in zip(leg.get_texts(), COLORS):
        txt.set_color(c)
        txt.set_fontweight("bold")
    y -= line_h * 3.4

    # ---- panels
    outer = fig.add_gridspec(3, 1, height_ratios=[1.45, 0.85, 0.85], hspace=0.34,
                             left=LEFT, right=RIGHT, top=y, bottom=0.115)
    top = outer[0].subgridspec(2, 1, height_ratios=[1.0, 0.16], hspace=0.10)
    axA = fig.add_subplot(top[0])
    axT = fig.add_subplot(top[1], sharex=axA)
    axB = fig.add_subplot(outer[1], sharex=axA)
    axC = fig.add_subplot(outer[2], sharex=axA)

    for ax in (axA, axB, axC):
        ax.set_facecolor(SURFACE)
        ax.grid(True, color=GRID, lw=0.7)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    panels = [(axA, 2, "A   How much of each linker lies on the fibril",
               f"PEG + biotin atoms touching\nthe fibril  (of {n_distal})"),
              (axB, 1, "B   How far each linker reaches",
               "SG to biotin S  (Å)"),
              (axC, 0, "C   How compact each linker is",
               "radius of gyration  (Å)")]
    for ax, j, ttl, ylab in panels:
        ax.set_title(ttl, loc="left", fontsize=12, fontweight="bold", color=INK,
                     pad=10)
        for k, c in enumerate(COLORS[:len(labels)]):
            yk = data[k][:, j]
            if w > 1:
                ax.plot(t_ns, yk, lw=0.5, color=c, alpha=0.14)
            ax.plot(t_ns, running_mean(yk, w), lw=2.0, color=c, solid_capstyle="round")
        ax.set_ylabel(ylab, fontsize=10, labelpad=6)
    axA.set_ylim(0, n_distal)
    for ax in (axA, axB):
        plt.setp(ax.get_xticklabels(), visible=False)
    axC.set_xlabel("time  (ns)", fontsize=10, labelpad=6)
    axC.set_xlim(t_ns[0], t_ns[-1])

    # ---- biotin strip under A: share of the time the biotin touches the fibril
    for k, c in enumerate(COLORS[:len(labels)]):
        frac = np.clip(running_mean(data[k][:, 3], w), 0, 1)
        cmap = LinearSegmentedColormap.from_list("b", [SURFACE, c])
        axT.imshow(frac[None, :], aspect="auto", cmap=cmap, vmin=0, vmax=1,
                   extent=[t_ns[0], t_ns[-1], len(labels) - 1 - k + 0.12,
                           len(labels) - k - 0.12], origin="lower",
                   interpolation="nearest")
    axT.set_ylim(0, len(labels))
    axT.set_yticks([len(labels) - 1 - k + 0.5 for k in range(len(labels))])
    axT.set_yticklabels([f"biotin, {s}" for s in short], fontsize=9)
    axT.tick_params(axis="y", length=0)
    axT.tick_params(axis="x", length=0)
    plt.setp(axT.get_xticklabels(), visible=False)
    for side in ("top", "right", "left", "bottom"):
        axT.spines[side].set_visible(False)
    axT.set_facecolor(SURFACE)

    if mark_ns is not None:
        for ax in (axA, axB, axC):
            ax.axvspan(mark_ns, t_ns[-1], color="#efeee9", lw=0, zorder=0)
        for ax in (axA, axT, axB, axC):
            ax.axvline(mark_ns, color=MUTED, lw=0.9, ls=(0, (4, 3)), zorder=1)
        axA.text(mark_ns + 0.01 * (t_ns[-1] - t_ns[0]), n_distal * 0.97,
                 f"from {mark_ns:.0f} ns", fontsize=9.5, color=INK2, ha="left", va="top")

    if note is None:
        note = (f"A contact is any heavy atom within {CONTACT:g} Å of a fibril heavy atom. "
                f"Thin lines: every frame ({dt:.0f} ps apart); thick lines: {SMOOTH_NS:g} ns "
                f"running mean. Biotin strips: share of the time the biotin touches the fibril, "
                f"{SMOOTH_NS:g} ns window (darker = more).")
    yn = 0.014
    for line in reversed(textwrap.wrap(note, width=max(30, int(avail_in / 0.072)))):
        fig.text(LEFT, yn, line, fontsize=9, color=MUTED, va="bottom")
        yn += line_h
    fig.savefig(path, dpi=170, facecolor=SURFACE)
    print(f"wrote {path}")
    print(f"wrote {path}")


def main():
    labels = meta["labels"]
    table = np.loadtxt("linked_eq.csv", delimiter=",", skiprows=1)
    t = np.round(table[:, 0] * 1000.0, 3)       # ps
    data = [table[:, 1 + 4 * k:5 + 4 * k] for k in range(len(labels))]
    plot(t, data, labels, meta["n_distal"], meta["mark_ns"], path=sys.argv[1] if len(sys.argv) > 1 else "linked_eq.png")


if __name__ == "__main__":
    main()
