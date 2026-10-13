# NCAA-fibril-toolbox

This toolbox is designed to allow users to run full atom trajectories of custom NCAA attached to custom fibrils. The aim
is to generate trajectories that allow MS-IBI to sample bonded and non-bonded interactions for a coarse-grained model.

It was written by Team Imperial College, iGEM 2026.

## Description

`fibriltoolbox` builds and runs all-atom simulations of a protein fibril that carries a non-canonical amino acid
(NCAA): a small molecule attached to a residue, such as a linker on a cysteine. It produces the three kinds of
trajectory a coarse-grained (MS-IBI) model is fitted to:

1. the fibril on its own,
2. the fibril with the NCAA attached at chosen residues,
3. the fibril surrounded by free copies of the capped NCAA.

Everything runs in the Amber ff19SB force field with GAFF2 for the small molecule and OPC water, in GROMACS.

The toolbox is the pipeline of our protocol, "ff19SB parameterisation pipeline: Cys–maleimide–PEG12–biotin" (35 steps,
worked out on a HET-s fibril with a Cys–maleimide–PEG12–biotin residue), with its steps 7 to 35 grouped into ten stages
that run from one config file. What was fixed in the protocol is a setting here: the starting structure, the residue and the
residue it derives from, the attachment sites, the restrained chains, the box.

What it does that the hand-run protocol does not:

- **No atom numbers typed in.** The backbone, the two caps and the atoms the residue shares with its parent amino acid
  are found in the molecule's bond graph and put on ff19SB names, types and charges.
- **Any number of attachment sites**, on any chains.
- **Optional stages.** Bring an already relaxed structure and start further down.
- **Box size from the structure**: every stage that makes a box works its size out by default, so that the fibril
  and whatever is attached to it cannot touch their own periodic images. Give a number instead to set it yourself.
- **The bookkeeping**: ion counts for 0.15 M, acpype's known repairs, restraint and index files, checks after every
  MD stage (with the expected density worked out from what is in the box), a log of every command, and a report
  per stage.

The protocol's own scripts are in `fibriltoolbox/scripts/` unchanged and still run by hand. The toolbox runs them
through a wrapper that replaces the settings at the top of a script with values from the config.


## Installation

### Requirements

| Software | Version | Used for |
|---|---|---|
| GROMACS | 2026.3 or newer (`gmx_mpi`) | every MD step; ships ff19SB as `amber19sb.ff` with OPC |
| AmberTools | conda-forge | antechamber, prepgen, parmchk2, tleap, ParmEd |
| acpype | 2026.9.3 or newer | Amber → GROMACS with 4-site OPC and the CMAP terms |
| PyMOL | 2 or 3, callable as `pymol` | building the fibril, docking a cap, grafting the residue (run headless) |
| Python | 3.11 or newer | with NumPy, NetworkX, ParmEd, MDAnalysis, matplotlib |

A GPU is assumed for the MD stages (`mdrun_flags` in the config).

### Install

```bash
conda activate <your environment with GROMACS, AmberTools and acpype>
git clone <this repository>
cd NCAA-fibril-toolbox
pip install -e .          # or: export PYTHONPATH=$PWD
python -m fibriltoolbox -h
```

## Usage

```bash
python -m fibriltoolbox <stage> [<stage> ...] config.toml
```

| Stage | Protocol steps | What it does |
|---|---|---|
| `ncaa-param` | 10–16 | capped NCAA `.mol2` → `XXX.prepin`, `frcmod.XXX`, `frcmod.ff19SB_XXX`, checked against the parent residue |
| `ncaa-relax` | 17 | the capped NCAA in water: the frame to graft, the frames to scatter, its reach |
| `fibril-prepare` | 7 | your protein-only structure → topology, box, water and ions *(optional)* |
| `fibril-run` | 7–9 | minimisation, NVT, NPT, production, checks → relaxed `fibril.pdb` *(optional)* |
| `linked-build` | 18–21, 23 | relaxed fibril + the NCAA at chosen sites → solvated GROMACS system |
| `linked-run` | 22, 24–26 | minimisation, NVT, NPT, production, periodic-image check |
| `linked-analyse` | 27 | linker conformations and adsorption |
| `free-build` | 28–32 | relaxed fibril + N free capped NCAAs → solvated GROMACS system |
| `free-run` | 33–34 | minimisation, NVT, NPT, production |
| `free-analyse` | 35 | settling, exchange, contacts |

The toolbox starts from a protein-only structure of the whole assembly. Designing that structure (protocol steps
1 to 6 for our fibril: template, mutation, n-mer, cap) is specific to each project and is not part of the toolbox.

### The worked example

[examples/hets_cyl/](examples/hets_cyl/) holds the protocol's own system:

- [fibril_design/](examples/hets_cyl/fibril_design/): how the capped HET-s fibril was designed by hand (protocol
  steps 1 to 6), with the PyMOL scripts used and the resulting `hets11_capped.pdb`;
- [config.toml](examples/hets_cyl/config.toml): the toolbox config that takes it from there, every setting commented.

```bash
cd examples/hets_cyl
python -m fibriltoolbox ncaa-param ncaa-relax config.toml
python -m fibriltoolbox fibril-prepare config.toml
nohup python -m fibriltoolbox fibril-run config.toml > fibril_run.out 2>&1 &
# when that has finished:
nohup python -m fibriltoolbox linked-build linked-run linked-analyse config.toml > run1.out 2>&1 &
nohup python -m fibriltoolbox free-build free-run free-analyse config.toml > run2.out 2>&1 &
```

### Your own system

Copy the example config and change:

- `[ncaa]`: `mol2` (your molecule drawn with an acetyl cap on the N and an N-methyl cap on the C), `name`
  (three characters), `net_charge`, `parent` (the amino acid it derives from, as named in ff19SB: `CYX`, `LYS` …).
- `[attach] sites`: chain and residue number of each attachment. After placing the residue, `linked-build` runs a
  clash test: a perfect backbone fit says nothing about whether the rest of the residue has room. If any of its
  heavy atoms lies within 1.5 Å of the protein at the best rotamer (`clash_within`: 1.5, 2.0 or 2.5), the stage stops
  and says the site is buried: choose an exposed residue or redesign the protein.
- the starting structure, in one of two ways:
  - `[relax] structure = "construct.pdb"`: your protein-only model, which `fibril-prepare` solvates and `fibril-run` relaxes;
  - `[linked] fibril = "relaxed.pdb"` and `[free] fibril = "relaxed.pdb"`: an already relaxed one, skipping those two stages.
- `restrain` in `[relax]`, `[linked]` and `[free]`: the interior chains held in production so the rod cannot turn.
- the box, in every stage that makes one: worked out from the structure by default (`"auto"`), or your own number.
  - `[relax] box_d`: nm between the protein and the box edge. Automatic: the protein's extent plus
    `image_margin_nm` (3.0), with any segment listed in `flexible` allowed to stretch out by its `reach_nm`.
  - `[linked] box`: buffer in Å. Automatic: each attached residue may point its full reach (`reach`: `"median"`,
    `"p95"`, `"max"` or nm; measured in `ncaa-relax`) in any direction, plus `image_margin_nm` (2.0).
  - `[free] box`: buffer in Å. Automatic: the fibril stays one residue's reach plus `image_margin_nm` (2.0) from
    its own image, so that no free molecule can touch both.

  The trade-off, which each of these stages also prints with an atom estimate for your system:

  - **The automatic box probably removes periodic-image artefacts, but costs computer time.** The run time per
    nanosecond grows roughly with the number of atoms. On the worked example: about 1.4 times the atoms of a 12 Å
    margin for the protein-only fibril (204,000 against 145,000), about 3 times for the fibril with two attached
    residues (440,000 against 145,000; 587,526 in another build of it), and no more for the free residues (518,000
    against 564,000).
  - **A small fixed margin, such as 12 Å, leads to periodic-image problems when a flexible region sits near the box
    edge.** An attached residue that reaches 4 nm, or a disordered tail, crosses 12 Å of water easily and meets the
    image of the structure on the other side. With 12 Å the worked example came within 0.22 nm of its image.

  After the build, each stage reports the box as built (x, y, z in nm), the water between your construct and the box
  edge on each side, and the distance from the construct to its own periodic image, whichever way the box was set.

  A number is the right choice when you know how far your structure can reach and what distance to its image you
  need. Either way, the periodic-image check after production says whether the box held.

### Checks that hold a run

The run stages (`ncaa-relax`, `fibril-run`, `linked-run`, `free-run`) apply the protocol's criteria after
minimisation, NVT and NPT (`check_em.py`, `check_nvt.py`, `check_npt.py`). If every line passes, the stage goes on to
the next run by itself. If one fails, it stops before the next run, prints the failed lines, and tells you where the
plot is. Look at it, then either:

- add `accept = ["npt"]` (or `"em"`, `"nvt"`) to that stage's table in the config and run the same command again:
  it carries on from where it stopped; or
- fix the cause and repeat the step.

`stop_before = "production"` makes a stage always stop before production, so the equilibration can be looked at
first. Minimisation does not hold `ncaa-relax` and `free-run`: the protocol expects it to end unconverged on those
systems and judges by the next step.

### Analysis and plotting are separate scripts

Every analysis writes its numbers to CSV files, and a second script draws the figure from those files. To change
how a figure looks, edit only the plotting script, or write your own against the CSV: the numbers, verdicts and
reports do not depend on it. A plotting script can be run again by hand in the folder that holds the CSV files.

| Analysis (numbers, verdicts, CSV) | Writes | Plotting (figure only) |
|---|---|---|
| `check_em.py` | `em.csv` | `plot_em.py` → `em.png` |
| `check_nvt.py` | `nvt.csv`, `nvt_meta.json` | `plot_nvt.py` → `nvt.png` |
| `check_npt.py` | `npt.csv`, `npt_meta.json` | `plot_npt.py` → `npt.png` |
| `check_capped_stability_3.py` | `stability.csv`, `stability_rmsd.csv`, `stability_meta.json` | `plot_capped_stability.py` → `stability.png` |
| `linked_eq_check.py` | `linked_eq.csv`, `linked_eq_meta.json` | `plot_linked_eq.py` → `linked_eq.png` |
| `free_linker_analysis.py` | four `free_linker_*.csv`, `free_linker_meta.json` | `plot_free_linker.py` → `free_linker_overview.png` |

The small `*_meta.json` beside a CSV holds what is not a table: the reference values, the averaging window, labels
and the verdicts the figure quotes. In the protocol document each of these pairs is still one script; here the
figure code was moved out unchanged. A figure that fails to draw is reported and does not stop a stage. In the
config, `plot_script` under `[relax.stability]`, `[linked.analysis]` or `[free.analysis]` names your own plotting
script.

### What a stage leaves behind

Each stage writes into its own folder under `workdir`:

- its files, named as in the protocol (`CYL.prepin`, `full_fib_noh.pdb`, `linked.amb2gmx/prod.xtc` …);
- `REPORT.md`: the checks and their results, step by step;
- `commands.log`: every command run, in order, with the settings handed to each script;
- `logs/`: the full output of each command;
- `steps_done.json`: finished steps. Run the stage again and it carries on after the last finished step; a
  production run that stopped is continued from its checkpoint.

A stage never deletes a file. Each `mdrun` starts in its own session with its process ID in `<name>.pid`, so it
outlives the terminal, and it can be stopped on its own with `kill $(cat prod.pid)`.

### Limits of this version

- The toolbox does not build the protein assembly: bring a protein-only PDB of it (see `examples/hets_cyl/fibril_design/` for how ours was made).
- `linked-analyse`, `free-analyse` and the stability check in `fibril-run` are the study's scripts. Their atom
  groups (cysteine and caps, spacer, PEG, biotin) describe our linker; another linker needs them adapted.
- The NCAA must attach through a side chain and derive from a canonical residue. Proline-like and
  Cα-disubstituted residues are not handled.

## Data and large files
This repository holds source code and the small input files of the worked example. The files of the study's own
runs (trajectories, topologies) are not here. [TODO: put them on Zenodo and link the record here.]

## Contributing
Pull requests are welcome. Every change is reviewed by the maintainer, Harvey Fung, who decides whether it is
merged.

To check a change, run the worked example's first stages, which need no GPU and take about a minute, and compare
`work/ncaa/CYL.prepin` and `work/ncaa/frcmod.CYL` with those from before your change:

```bash
cd examples/hets_cyl && python -m fibriltoolbox ncaa-param config.toml
```

## Authors and acknowledgment

**Authors**

- Harvey Fung (Team Imperial College, iGEM 2026)
- Claude Code (Anthropic's AI coding assistant), which drafted the code in `fibriltoolbox/` from the team's protocol

**Please cite**

The toolbox drives other people's programs. These are the references those programs ask to be cited, as they
print them when the pipeline runs.

GROMACS, printed in every `mdrun` log:

1. M. J. Abraham, T. Murtola, R. Schulz, S. Páll, J. C. Smith, B. Hess, E. Lindahl. GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers. *SoftwareX* 1–2, 19–25 (2015). doi:10.1016/j.softx.2015.06.001
2. S. Páll, M. J. Abraham, C. Kutzner, B. Hess, E. Lindahl. Tackling exascale software challenges in molecular dynamics simulations with GROMACS. In S. Markidis, E. Laure (eds.), *Solving Software Challenges for Exascale*, Lecture Notes in Computer Science, 3–27 (2015). doi:10.1007/978-3-319-15976-8_1
3. S. Pronk, S. Páll, R. Schulz, P. Larsson, P. Bjelkmar, R. Apostolov, M. R. Shirts, J. C. Smith, P. M. Kasson, D. van der Spoel, B. Hess, E. Lindahl. GROMACS 4.5: a high-throughput and highly parallel open source molecular simulation toolkit. *Bioinformatics* 29, 845–854 (2013). doi:10.1093/bioinformatics/btt055
4. B. Hess, C. Kutzner, D. van der Spoel, E. Lindahl. GROMACS 4: Algorithms for highly efficient, load-balanced, and scalable molecular simulation. *J. Chem. Theory Comput.* 4, 435–447 (2008). doi:10.1021/ct700301q
5. D. van der Spoel, E. Lindahl, B. Hess, G. Groenhof, A. E. Mark, H. J. C. Berendsen. GROMACS: Fast, flexible and free. *J. Comput. Chem.* 26, 1701–1718 (2005). doi:10.1002/jcc.20291
6. E. Lindahl, B. Hess, D. van der Spoel. GROMACS 3.0: A package for molecular simulation and trajectory analysis. *J. Mol. Model.* 7, 306–317 (2001). doi:10.1007/s008940100045
7. H. J. C. Berendsen, D. van der Spoel, R. van Drunen. GROMACS: A message-passing parallel molecular dynamics implementation. *Comput. Phys. Commun.* 91, 43–56 (1995). doi:10.1016/0010-4655(95)00042-E

GROMACS, for the algorithms the runs use:

8. U. Essmann, L. Perera, M. L. Berkowitz, T. Darden, H. Lee, L. G. Pedersen. A smooth particle mesh Ewald method. *J. Chem. Phys.* 103, 8577–8593 (1995). doi:10.1063/1.470117
9. S. Miyamoto, P. A. Kollman. SETTLE: An analytical version of the SHAKE and RATTLE algorithm for rigid water models. *J. Comput. Chem.* 13, 952–962 (1992). doi:10.1002/jcc.540130805
10. B. Hess, H. Bekker, H. J. C. Berendsen, J. G. E. M. Fraaije. LINCS: A linear constraint solver for molecular simulations. *J. Comput. Chem.* 18, 1463–1472 (1997). doi:10.1002/(SICI)1096-987X(199709)18:12<1463::AID-JCC4>3.0.CO;2-H
11. G. Bussi, D. Donadio, M. Parrinello. Canonical sampling through velocity rescaling. *J. Chem. Phys.* 126, 014101 (2007). doi:10.1063/1.2408420
12. M. Bernetti, G. Bussi. Pressure control using stochastic cell rescaling. *J. Chem. Phys.* 153, 114107 (2020). doi:10.1063/5.0020514

sqm (AmberTools), which computes the AM1 charges in `ncaa-param`:

13. R. C. Walker, M. F. Crowley, D. A. Case. The implementation of a fast and accurate QM/MM potential method in Amber. *J. Comput. Chem.* 29, 1019–1031 (2008). doi:10.1002/jcc.20857

GROMACS tools the toolbox also runs (solvation and insertion, a net-charge notice, `gmx hbond`), given as printed:

14. A. Bondi. van der Waals Volumes and Radii. *J. Phys. Chem.* (1964). doi:10.1021/j100785a001
15. J. S. Hub, B. L. de Groot, H. Grubmüller, G. Groenhof. Quantifying Artifacts in Ewald Simulations of Inhomogeneous Systems with a Net Charge. *J. Chem. Theory Comput.* (2014). doi:10.1021/ct400626b
16. S. Gorelov, A. Titov, O. Tolicheva, A. Konevega, A. Shvetsov. Determination of Hydrogen Bonds in GROMACS: A New Implementation to Overcome Memory Limitation. *J. Chem. Inf. Model.* (2024). doi:10.1021/acs.jcim.3c02087

The structure the worked example starts from, PDB entry 2KJ3:

17. H. Van Melckebeke, C. Wasmer, A. Lange, E. Ab, A. Loquet, A. Böckmann, B. H. Meier. Atomic-resolution three-dimensional structure of HET-s(218–289) amyloid fibrils by solid-state NMR spectroscopy. *J. Am. Chem. Soc.* 132, 13765 (2010).

References 1–13 are those listed in the team's protocol, "ff19SB parameterisation pipeline:
Cys–maleimide–PEG12–biotin". The toolbox also uses AmberTools (antechamber, parmchk2, tleap), acpype, the ff19SB,
GAFF2 and OPC models, PyMOL and MDAnalysis; these did not print a reference in our runs, so cite them as their own
documentation asks.

## License
This repository is licensed under the [Apache License 2.0](LICENSE) — a permissive
open-source license recommended for software (Creative Commons licenses are *not*
intended for source code). You are free to use, modify, and distribute this software,
provided you keep the license and attribution notices. If you prefer different terms
for your own tool, you may replace this license, but it must remain an
[OSI-approved open-source license](https://opensource.org/licenses).
