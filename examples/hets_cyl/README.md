# Example: HET-s fibril with Cys–maleimide–PEG12–biotin

An 11-chain HET-s(218–289) amyloid fibril with an engineered "poison" cap docked on one end, carrying the
non-canonical residue CYL (cysteine–maleimide–PEG12–biotin) at two sites: Cys239 of chain A, and residue 130 of the
cap (chain L). The example holds the config for the whole pipeline and the results of the study's own runs, one
folder per toolbox stage, laid out and named as the toolbox writes them.

For a step-by-step guide to producing all of it, the trajectories included, see [WALKTHROUGH.md](WALKTHROUGH.md).

## What is here

| Folder | Toolbox steps | Holds |
|---|---|---|
| `config.toml` | all | every setting of the runs below |
| `protein_relax/` | 1–4 | the capped fibril on its own, 150 ns: starting structure, equilibration checks, stability and interface analyses, the relaxed structure |
| `protein_relax/fibril_design/` | before step 1 | how the capped fibril was designed, by hand in PyMOL and AlphaFold: [its own README](protein_relax/fibril_design/README.md) |
| `ncaa/` | 5–11 | the CYL molecule (`12EO_linker.mol2`) and its parameters: `CYL.prepin`, `frcmod.CYL`, `frcmod.ff19SB_CYL` |
| `ncaa_relax/` | 12 | capped CYL alone in water, 1 ns: the frame grafted onto the fibril (`cyl.pdb`), the frames scattered in Run 2, its reach |
| `linked/` | 13–21 | Run 1, CYL attached at the two sites, 200 ns: the built structure, checks, image distance, how each CYL behaves |
| `free/` | 22–29 | Run 2, the fibril among 40 free capped CYL, 200 ns: the built structure, checks, settling, exchange and contacts |

Each stage folder has the minimisation, NVT and NPT checks (`potential.xvg`, `nvt.xvg`, `npt.xvg`, their `.csv` and
`.png`), the GROMACS logs and energy files, the structures that start and end each run, and the analysis results.
The trajectories are not included: they are several GB per run.

## Running it

From this folder:

```bash
python -m ncaaprotein protein-prepare config.toml
nohup python -m ncaaprotein protein-run config.toml > protein_run.out 2>&1 &
python -m ncaaprotein ncaa-param ncaa-relax config.toml
nohup python -m ncaaprotein linked-build linked-run linked-analyse config.toml > run1.out 2>&1 &
nohup python -m ncaaprotein free-build free-run free-analyse config.toml > run2.out 2>&1 &
```

Your run goes into `work/`, next to the study's folders, so each file can be compared with its counterpart here.
Keep `workdir` a separate folder: set to `"."`, the toolbox would write into these folders and over the study's
results. Set `pymol` under `[tools]` to your PyMOL if it is not on your `PATH`.

The config starts Run 1 and Run 2 from the fibrils the study used (`[linked] protein`, `[free] protein`). Leave those
two lines out to start both from the end of your own `protein-run` instead.

## How the study's runs were set up

- **Capped fibril:** 150 ns, nothing restrained in production, 12 Å of water around the fibril.
- **Run 1:** built from the 9.5 ns frame of an earlier capped-fibril run (`linked/hets11_capped_9.5ns.pdb`), with
  Met217 kept in every chain, and every atom of chains F and G restrained at 200 kJ/mol/nm².
- **Run 2:** built from the 101.2 ns frame of the capped-fibril run (`free/hets11_capped_101.2ns.pdb`), with Met217
  removed, and the backbone of chains F and G restrained at 200 kJ/mol/nm².
- **Both runs:** 12 Å buffer, 0.15 M NaCl, 298 K, 4 fs steps with repartitioned hydrogen masses.

The config removes Met217 in both builds (`[remove] when = "builds"`), so Run 1 rebuilt from it lacks that residue in
every chain, unlike the study's. Met217 comes from the expression construct of PDB 2KJ3 and was kept in the
capped-fibril run; its S–CH₃ trips LINCS warnings at 4 fs, which is why the config sets
`mdrun_env = {GMX_MAXCONSTRWARN = "-1"}`.

## Results

| What | Result | File |
|---|---|---|
| Capped fibril, interfaces | every chain junction keeps 75–92% of its contacts over 150 ns (F–G lowest); the cap junction (K–L) is about 1.2 Å wider than the fibril's own and holds about half the backbone hydrogen bonds (5.9 against 11.5 at the free end A–B), but stays below the width at which it would count as coming apart | `protein_relax/stability.png`, `assembly_stability.png` |
| Capped fibril, hydrogen bonds per interface (75–150 ns) | A–B 21.7, E–F 18.1, F–G 19.6, H–I 22.0, K–L 11.8 | `protein_relax/hb_*.xvg` |
| Image distance, capped fibril | **FAIL**: 1.04 nm at 65.9 ns, under 1.2 nm in 9 of 1,501 frames, all between 58 and 66 ns; clean after 66 ns | `protein_relax/mindist_pi.xvg` |
| Image distance, Run 1 | **FAIL**: 0.91 nm at 20.7 ns, under 1.2 nm in 7 of 2,001 frames | `linked/linked.amb2gmx/mindist_pi.xvg` |
| CYL in water | reach from Cα: median 3.38 nm, 95% 4.13 nm, largest 4.38 nm; fully stretched 7.07 nm | `ncaa_relax/ace_cyl_nme.amb2gmx/reach.json` |
| Run 1, the attached CYL | both lie on the fibril about 78% of the time, with the biotin tip on it 75% (chain A) and 63% (cap) of the time | `linked/linked.amb2gmx/linked_eq.png`, `linked_eq_report.txt` |
| Run 2, the free CYL | 9.0 on the fibril on average; 32 of 40 visited it, with 84 arrivals and 80 departures; 2 stayed on it for more than 90% of the run | `free/free.amb2gmx/free_ncaa_overview.png`, `free_ncaa_report.txt` |

Two things to know when reading these:

- **The 12 Å buffer was too small.** The capped fibril and Run 1 both came within 1.0–1.2 nm of their own periodic
  image, inside the reach of the 1.0 nm cutoff. The toolbox's default automatic box (`box = "auto"`) avoids this, at
  the cost of a larger system.
- **The solute runs about 4 K warm at 4 fs** (301.4–302.3 K against 298 K, with the water at 297.6–297.9 K), the
  known effect of a large time step with repartitioned hydrogen masses. In the free run's NPT the check fails on this,
  and the config accepts it (`[free] accept = ["npt"]`).

## Files in a different form

- `prod.log.gz`: the production logs, gzipped (`gunzip prod.log.gz`).
- `prod_0-100ns.edr` and `prod_100-200ns.edr` (Run 1 and Run 2): the production energy files, split at 100 ns.
  Join them with `gmx eneconv -f prod_0-100ns.edr prod_100-200ns.edr -o prod.edr`. The temperatures, pressure and
  density are also in `prod_energy.xvg`.
- The capped fibril's production energies are in `protein_relax/check_energy.xvg`.

## Scripts specific to this example

- `protein_relax/check_capped_stability.py` and `plot_capped_stability.py`: the stability check written for a capped
  HET-s rod, run in step 3 through `[relax.stability]`. It compares the cap junction with the fibril's own junctions.
- `protein_relax/image_contacts.py`: which part of the fibril comes close to its own periodic image.
- `protein_relax/fibril_design/`: the five PyMOL scripts of the design.
