#!/usr/bin/env python3
"""
Judge whether a GROMACS NVT equilibration equilibrated, and write its traces to nvt.csv.
The figure is drawn from nvt.csv by plot_nvt.py, a separate script.

Extract the data first (the T- terms are the per-group temperatures; leave them
out for a single tc-grps = System run):
    printf "Temperature\nT-Solute\nT-Water_and_ions\nPotential\nPressure\n0\n" \
        | gmx_mpi energy -f nvt.edr -o nvt.xvg

Then:
    python3 check_nvt.py            # or: python3 check_nvt.py nvt.xvg

Reads the legend from the .xvg header, so column order does not matter and
extra terms are handled automatically.
"""

import sys
import numpy as np

REF_T = 298.0          # <-- must match ref_t in nvt.mdp
LARGE_STEP_MDP = ""    # set by the toolbox when this run used a time step above 2 fs with repartitioned hydrogen
                       # masses: the mdp file to change if the solute is off ref_t while the water is not
EQ_FRACTION = 0.5      # analyse the last half as the "equilibrated" window

xvg = sys.argv[1] if len(sys.argv) > 1 else "nvt.xvg"


def stats(t, y):
    """mean, residual std, slope (per ps) and its standard error, effective N.
    Neighbouring frames are correlated, so N is reduced by the lag-1
    autocorrelation of the residuals before any uncertainty is computed."""
    a, b = np.polyfit(t, y, 1)
    res = y - (a * t + b)
    r1 = np.corrcoef(res[:-1], res[1:])[0, 1] if len(y) > 3 else 0.0
    r1 = min(max(r1, 0.0), 0.99)
    neff = max(3.0, len(y) * (1 - r1) / (1 + r1))
    se = res.std() / np.sqrt(neff) / t.std()
    return y.mean(), res.std(), a, se, neff


def settled(t, y):
    """drift across the window is inside the noise, or not distinguishable from zero"""
    mean, sd, a, se, neff = stats(t, y)
    drift = a * (t[-1] - t[0])
    ok = abs(drift) < sd or abs(a) < 3 * se
    return ok, f"drift {drift:+.3g} over window = {abs(drift) / sd:.2f} sigma, {abs(a) / se:.1f} x its error"


# ------------------------------------------------------------ read .xvg
labels, rows = [], []
with open(xvg) as fh:
    for line in fh:
        s = line.strip()
        if s.startswith("@") and " legend " in s and s.startswith("@ s"):
            labels.append(s.split('"')[1])
        elif s and not s.startswith(("#", "@")):
            rows.append([float(v) for v in s.split()])

data = np.asarray(rows)
time = data[:, 0]
series = {lab: data[:, i + 1] for i, lab in enumerate(labels)}

if not series:
    sys.exit(f"No data columns found in {xvg} - check the file.")

print(f"read {len(time)} frames, {time[0]:.1f} to {time[-1]:.1f} ps")
print("terms:", ", ".join(series) + "\n")

cut = int(len(time) * (1 - EQ_FRACTION))

# ------------------------------------------------------------ data for plot_nvt.py
import csv, json
with open("nvt.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["time_ps"] + labels)
    for r in data:
        w.writerow([repr(float(v)) for v in r[:len(labels) + 1]])
json.dump(dict(REF_T=REF_T, EQ_FRACTION=EQ_FRACTION), open("nvt_meta.json", "w"), indent=1)
print("wrote nvt.csv and nvt_meta.json  (plot them with plot_nvt.py)\n")

# ------------------------------------------------------------- verdict
print(f"{'term':<18}{'mean':>14}{'std':>12}{'drift/ns':>12}{'+/-':>10}")
print("-" * 66)

verdict = []
info = []
for lab in labels:
    y = series[lab][cut:]
    t = time[cut:]
    mean, sd, a, se, neff = stats(t, y)
    print(f"{lab:<18}{mean:>14.6g}{sd:>12.4g}{a * 1000:>12.4g}{se * 1000:>10.3g}")

    if "Temperature" in lab or lab.startswith("T-"):
        off, sem = abs(mean - REF_T), sd / np.sqrt(neff)
        verdict.append((f"{lab} mean at ref_t", off < max(1.0, 3 * sem),
                        f"off by {off:.2f} K, standard error {sem:.2f} K"))
        verdict.append((f"{lab} not drifting",) + settled(t, y))
    if "Potential" in lab:
        info.append(("potential energy",) + settled(t, y))

# a temperature line that failed: is it the solute or the water?
t_fail = {n for n, ok, _ in verdict if n.startswith("T-") and n.endswith("mean at ref_t") and not ok}
t_all = {n for n, _, _ in verdict if n.startswith("T-") and n.endswith("mean at ref_t")}
water = {n for n in t_all if "Water" in n or "SOL" in n}
for k, (n, ok, detail) in enumerate(verdict):
    if n not in t_fail:
        continue
    if n in water:
        detail += ("; the water is off ref_t, which the time step does not explain: look at the run (bad contacts, "
                   "the thermostat settings) before going on")
    elif water and not (water & t_fail) and LARGE_STEP_MDP:
        detail += (f"; the solute is off ref_t while the water is not. That is the known effect of a large time step "
                   f"with repartitioned hydrogen masses: run this step at 2 fs (dt = 0.002 and mass-repartition-factor = 1 "
                   f"in your copy of {LARGE_STEP_MDP}; README, \"Your own mdp files\")")
    verdict[k] = (n, ok, detail)

print("\nverdict")
print("-" * 66)
for name, ok, detail in verdict:
    print(f"  [{'PASS' if ok else 'FAIL'}]  {name}  ({detail})")
for name, ok, detail in info:
    print(f"  [info]  {name} {'settled' if ok else 'still drifting'}  ({detail})")

print("""
Notes
  A drift only fails if it is both larger than the noise (1 sigma across the
  window) and clearly real (more than 3 x its own standard error). A per-ns
  slope from a 50 ps window is mostly noise, so it is not used as a threshold.
  Pressure fluctuates enormously (hundreds of bar) in NVT and is not a
  useful equilibration criterion here; judge it in the NPT stage instead.
  The potential energy is shown for information only. At fixed volume it keeps
  drifting until the box reaches its density, which a tleap-built box only
  does in NPT; check_npt.py judges it there.
""")
