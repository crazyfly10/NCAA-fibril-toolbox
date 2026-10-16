#!/usr/bin/env python3
"""
Draw nvt.png from nvt.csv, the traces written by check_nvt.py.

Usage:
    python3 plot_nvt.py                       # or: python3 plot_nvt.py nvt.csv nvt.png

The reference values and the averaging window are read from nvt_meta.json, which check_nvt.py
writes beside the CSV, so the figure marks what the check judged against.
This script only draws. Change it freely: the numbers and the verdict come from check_nvt.py.
"""

import csv
import json
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

src = sys.argv[1] if len(sys.argv) > 1 else "nvt.csv"
out = sys.argv[2] if len(sys.argv) > 2 else "nvt.png"
meta_file = os.path.splitext(src)[0] + "_meta.json"
meta = json.load(open(meta_file)) if os.path.exists(meta_file) else {}
REF_T = meta.get("REF_T", 298.0)
EQ_FRACTION = meta.get("EQ_FRACTION", 0.5)

with open(src) as fh:
    rows = list(csv.reader(fh))
labels = rows[0][1:]
data = np.array([[float(v) for v in r] for r in rows[1:]])
time = data[:, 0]
series = {lab: data[:, i + 1] for i, lab in enumerate(labels)}
cut = int(len(time) * (1 - EQ_FRACTION))

# ------------------------------------------------------------ plotting
n = len(labels)
fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 4.2), squeeze=False)
axes = axes[0]


for ax, lab in zip(axes, labels):
    y = series[lab]
    ax.plot(time, y, lw=0.7, alpha=0.5, color="#1f77b4")

    # running mean to make drift visible through the noise
    w = max(1, len(y) // 50)
    if w > 1:
        ker = np.ones(w) / w
        ax.plot(time[w - 1:], np.convolve(y, ker, mode="valid"),
                lw=1.8, color="#d62728", label=f"running mean ({w} pts)")

    ax.axvline(time[cut], ls=":", color="grey")
    ax.axhline(y[cut:].mean(), ls="--", color="k", lw=1,
               label=f"mean(last {int(EQ_FRACTION*100)}%) = {y[cut:].mean():.3g}")

    if "Temperature" in lab or lab.startswith("T-"):
        ax.axhline(REF_T, ls="-", color="green", lw=1, label=f"ref_t = {REF_T}")

    ax.set_xlabel("time (ps)")
    ax.set_ylabel(lab)
    ax.set_title(lab)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig(out, dpi=150)
print(f"wrote {out}")
