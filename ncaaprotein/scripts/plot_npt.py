#!/usr/bin/env python3
"""
Draw npt.png from npt.csv, the traces written by check_npt.py.

Usage:
    python3 plot_npt.py                       # or: python3 plot_npt.py npt.csv npt.png

The reference values and the averaging window are read from npt_meta.json, which check_npt.py
writes beside the CSV, so the figure marks what the check judged against.
This script only draws. Change it freely: the numbers and the verdict come from check_npt.py.
"""

import csv
import json
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

src = sys.argv[1] if len(sys.argv) > 1 else "npt.csv"
out = sys.argv[2] if len(sys.argv) > 2 else "npt.png"
meta_file = os.path.splitext(src)[0] + "_meta.json"
meta = json.load(open(meta_file)) if os.path.exists(meta_file) else {}
REF_T = meta.get("REF_T", 298.0)
EQ_FRACTION = meta.get("EQ_FRACTION", 0.5)
REF_P = meta.get("REF_P", 1.0)
RHO_LO, RHO_HI = meta.get("RHO_LO", 995.0), meta.get("RHO_HI", 1020.0)

with open(src) as fh:
    rows = list(csv.reader(fh))
labels = rows[0][1:]
data = np.array([[float(v) for v in r] for r in rows[1:]])
time = data[:, 0]
series = {lab: data[:, i + 1] for i, lab in enumerate(labels)}
cut = int(len(time) * (1 - EQ_FRACTION))

# ------------------------------------------------------------ plotting
n = len(labels)
ncol = min(n, 2)
nrow = int(np.ceil(n / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(6.5 * ncol, 4.2 * nrow),
                         squeeze=False)
flat = axes.ravel()

for ax, lab in zip(flat, labels):
    y = series[lab]
    ax.plot(time, y, lw=0.7, alpha=0.45, color="#1f77b4")

    w = max(1, len(y) // 50)
    if w > 1:
        ax.plot(time[w - 1:], np.convolve(y, np.ones(w) / w, mode="valid"),
                lw=1.8, color="#d62728", label=f"running mean ({w} pts)")

    m = y[cut:].mean()
    ax.axvline(time[cut], ls=":", color="grey")
    ax.axhline(m, ls="--", color="k", lw=1,
               label=f"mean(last {int(EQ_FRACTION*100)}%) = {m:.5g}")

    if "Temperature" in lab or lab.startswith("T-"):
        ax.axhline(REF_T, color="green", lw=1, label=f"ref_t = {REF_T}")
    if lab.startswith("Pressure"):
        ax.axhline(REF_P, color="green", lw=1, label=f"ref_p = {REF_P}")
    if "Density" in lab:
        ax.axhspan(RHO_LO, RHO_HI, color="green", alpha=0.12,
                   label=f"expected {RHO_LO:g}-{RHO_HI:g}")

    ax.set_xlabel("time (ps)")
    ax.set_ylabel(lab)
    ax.set_title(lab)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

for ax in flat[n:]:
    ax.axis("off")

plt.tight_layout()
plt.savefig(out, dpi=150)
print(f"wrote {out}")
