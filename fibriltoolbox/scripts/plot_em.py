#!/usr/bin/env python3
"""
Draw em.png from em.csv, the energy curve written by check_em.py.

Usage:
    python3 plot_em.py                        # or: python3 plot_em.py em.csv em.png

This script only draws. Change it freely: the numbers and the verdict come from check_em.py.
"""

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

src = sys.argv[1] if len(sys.argv) > 1 else "em.csv"
out = sys.argv[2] if len(sys.argv) > 2 else "em.png"

data = np.loadtxt(src, delimiter=",", skiprows=1)
step, pot = data[:, 0], data[:, 1]

# ---------------------------------------------------------------- plot
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))

ax[0].plot(step, pot, lw=1.0, color="#1f77b4")
ax[0].set_xlabel("step")
ax[0].set_ylabel("potential energy (kJ/mol)")
ax[0].set_title("full minimisation")

# log-scale view of the distance still to fall, reveals early behaviour
floor = pot.min()
above = pot - floor + 1.0
ax[1].semilogy(step, above, lw=1.0, color="#d62728")
ax[1].set_xlabel("step")
ax[1].set_ylabel("E - E_min + 1  (kJ/mol)")
ax[1].set_title("approach to minimum (log)")

# last 20% - the panel that actually shows convergence
n = max(10, len(step) // 5)
ax[2].plot(step[-n:], pot[-n:], lw=1.0, color="#2ca02c")
ax[2].set_xlabel("step")
ax[2].set_ylabel("potential energy (kJ/mol)")
ax[2].set_title("last 20% (zoom)")

for a in ax:
    a.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(out, dpi=150)
print(f"wrote {out}")
