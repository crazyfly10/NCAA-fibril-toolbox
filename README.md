# NCAA-protein-toolbox

This toolbox is designed to allow users to run full atom trajectories of custom NCAAs attached to custom proteins. The aim
is to generate trajectories that allow MS-IBI to sample bonded and non-bonded interactions for a coarse-grained model.

It was written by Team Imperial College, iGEM 2026.

## Description

`ncaaprotein` builds and runs all-atom simulations of a protein, or an assembly of proteins, that carries a
non-canonical amino acid (NCAA): a small molecule attached to a residue, such as a synthetic side chain on a
cysteine. It produces the three kinds of trajectory a coarse-grained (MS-IBI) model is fitted to:

1. the protein on its own,
2. the protein with the NCAA attached at chosen residues,
3. the protein surrounded by free copies of the capped NCAA.

Everything runs in the Amber ff19SB force field with GAFF2 for the small molecule and OPC water, in GROMACS.

The toolbox groups the work into ten stages that run from one config file. Everything about your system is a
setting there: the starting structure, the residue and the residue it derives from, the attachment sites, the
restrained chains, the box.

What it takes care of:

- **No atom numbers typed in.** The backbone, the two caps and the atoms the residue shares with its parent amino acid
  are found in the molecule's bond graph and put on ff19SB names, types and charges.
- **Any number of attachment sites**, on any chains.
- **Optional stages.** Bring an already relaxed structure and start further down.
- **Box size from the structure**: every stage that makes a box works its size out by default, so that the protein
  and whatever is attached to it cannot touch their own periodic images. Give a number instead to set it yourself.
- **The bookkeeping**: ion counts for 0.15 M, acpype's known repairs, restraint and index files, checks after every
  MD stage (with the expected density worked out from what is in the box), a log of every command, and a report
  per stage.

The analysis and helper scripts are in `ncaaprotein/scripts/` and can also be run by hand. The toolbox runs
them with their settings taken from the config.


## Installation

### Requirements

| Software | Version | Used for |
|---|---|---|
| GROMACS | 2026.3 or newer (`gmx_mpi`) | every MD step; ships ff19SB as `amber19sb.ff` with OPC |
| AmberTools | conda-forge | antechamber, prepgen, parmchk2, tleap, ParmEd |
| acpype | 2026.9.3 or newer | Amber → GROMACS with 4-site OPC and the CMAP terms |
| PyMOL | 2 or 3, callable as `pymol` | placing the residue on the protein in `linked-build` (run headless) |
| Python | 3.11 or newer | with NumPy, NetworkX, ParmEd, MDAnalysis, matplotlib |

A GPU is assumed for the MD stages (`mdrun_flags` in the config).

### Install

```bash
conda activate <your environment with GROMACS, AmberTools and acpype>
git clone <URL of this repository>
cd <the cloned folder>
pip install -e .          # or: export PYTHONPATH=$PWD
python -m ncaaprotein -h
```

## Usage

```bash
python -m ncaaprotein <stage> [<stage> ...] config.toml
```

| Stage | Steps | What it does |
|---|---|---|
| `protein-prepare` | 1 | your protein-only structure → topology, box, water and ions *(optional)* |
| `protein-run` | 2–4 | minimisation, NVT, NPT, production, checks → relaxed `protein.pdb` *(optional)* |
| `ncaa-param` | 5–11 | capped NCAA `.mol2` → `XXX.prepin`, `frcmod.XXX`, `frcmod.ff19SB_XXX`, checked against the parent residue |
| `ncaa-relax` | 12 | the capped NCAA in water: the frame to graft, the frames to scatter, its reach |
| `linked-build` | 13–17 | relaxed protein + the NCAA at chosen sites → solvated GROMACS system |
| `linked-run` | 18–20 | minimisation, NVT, NPT, production, periodic-image check |
| `linked-analyse` | 21 | NCAA conformations and adsorption |
| `free-build` | 22–26 | relaxed protein + N free capped NCAAs → solvated GROMACS system |
| `free-run` | 27–28 | minimisation, NVT, NPT, production |
| `free-analyse` | 29 | settling, exchange, contacts |

The toolbox starts from a protein-only structure of the whole assembly. Designing that structure is specific
to each project and is not part of the toolbox.

### Quick start

```bash
mkdir my_run && cp config_template.toml my_run/config.toml     # then put your structure and molecule in my_run/
cd my_run
# fill in the lines marked <...> in config.toml, then:
python -m ncaaprotein protein-prepare config.toml
nohup python -m ncaaprotein protein-run config.toml > protein_run.out 2>&1 &
python -m ncaaprotein ncaa-param ncaa-relax config.toml
# when both have finished:
nohup python -m ncaaprotein linked-build linked-run linked-analyse config.toml > run1.out 2>&1 &
nohup python -m ncaaprotein free-build free-run free-analyse config.toml > run2.out 2>&1 &
```

You bring two files: a protein-only PDB of the whole assembly, and a `.mol2` of the molecule drawn capped. The
first two lines of commands do not depend on each other and can run in either order. Several stages given in one
command run one after the other.

[config_template.toml](config_template.toml) lists every setting with its default.

### Your own system

In your copy of the config, set:

- `[ncaa]`: `mol2` (your molecule drawn with an acetyl cap on the N and an N-methyl cap on the C), `name`
  (three characters), `net_charge`, `parent` (the amino acid it derives from, as named in ff19SB: `CYX`, `LYS` …).
- `[attach] sites`: chain and residue number of each attachment. **You choose these yourself.** The toolbox does not
  look for attachment sites: open your structure in PyMOL (or another viewer), decide which residues the molecule
  should replace, and give their chain letters and residue numbers as they appear in that file. Pick residues whose
  side chains point into the solvent; measuring their solvent-accessible area first saves a failed build. The toolbox replaces whatever residue is at the position you give and does
  not check its type, and a glycine cannot be a site because the fit needs a CB. After placing the residue, `linked-build` runs a
  clash test: a perfect backbone fit says nothing about whether the rest of the residue has room. If any of its
  heavy atoms lies within 1.5 Å of the protein at the best rotamer (`clash_within`: 1.5, 2.0 or 2.5), the stage stops
  and says the site is buried: choose an exposed residue or redesign the protein.
- the starting structure, in one of two ways:
  - `[relax] structure = "construct.pdb"`: your protein-only model, which `protein-prepare` solvates and `protein-run` relaxes;
  - `[linked] protein = "relaxed.pdb"` and `[free] protein = "relaxed.pdb"`: an already relaxed one, skipping those two stages.
- **an N-terminal methionine, if you designed the protein from a deposited structure.** Structures from the PDB
  often start with the methionine of the expression construct. Whether it is still there in your protein depends
  on how it is made, so if you are not sure, you can keep it. Keeping it has one consequence here: production runs
  use 4 fs steps with repartitioned hydrogen masses, and the S–CH₃ of a methionine at a flexible chain end can turn
  far enough in one step to raise LINCS warnings; `gmx mdrun` stops after 1,000 of them. If your structure keeps
  such a methionine, do one of two things:
  - set `mdrun_env = {GMX_MAXCONSTRWARN = "-1"}` under `[tools]`. This lifts the limit, and it also silences every
    other constraint warning, including one that would have mattered, so read the LINCS count each stage reports;
  - or run at 2 fs without mass repartitioning: in your own copy of the mdp files (see "Your own mdp files"
    below), set `dt = 0.002` and `mass-repartition-factor = 1` in `protein_prod.mdp`, `linked_prod.mdp`,
    `free_prod.mdp` and `free_npt.mdp`. Runs then take about twice as long; the toolbox works the number of steps
    out from `dt`.

  If the methionine is not in your construct, remove it with `[remove]` rather than by editing the structure:

  ```toml
  [remove]
  resn = "MET"
  resi = 217          # residue number as in your structure; taken out of every chain that has it
  # when = "all"      # "all" (default): from protein-prepare on, so no stage has it. "builds": only in linked-build
                      # and free-build, e.g. for a relaxed structure given with [linked]/[free] protein that still has it
  ```
- `restrain` in `[relax]`, `[linked]` and `[free]`: what is held in production so the assembly cannot turn in its
  box. Leave it out to hold nothing. It works the way GROMACS position restraints do: an index group of atoms, each
  held to its starting position with force constant `fc` (kJ/mol/nm²). The simple form names chains and atoms:

  ```toml
  restrain = {chains = ["F", "G"], atoms = "N CA C O", fc = 200}
  ```

  `chains` takes whole chains (`"G"`) or stretches of one (`"G:226-282"`, residue numbers as in your structure), any
  number of them. `atoms` is a list of atom names, or `"heavy"` (every atom but hydrogens) or `"all"`. For anything
  else, give the group as a `gmx select` expression instead (see the GROMACS documentation of `gmx select`, or of
  `gmx make_ndx` for the same groups built interactively):

  ```toml
  restrain = {select = "resnr 226 to 282 and name CA", name = "core", fc = 1000}
  ```

  In `[relax]` the expression is applied to the protein as pdb2gmx built it (residue numbers as in your structure,
  chain by chain); in `[linked]` and `[free]` to the system tleap built, where residues are numbered 1 to N straight
  through the protein and its NCAAs. `name` labels the restraint files and the mdp `define`.
- `salt` in `[relax]`, `[linked]` and `[free]` (default 0.15 M). In `protein-prepare` this goes to
  `gmx genion -conc`. The NCAA stages build their systems in tleap instead, because the NCAA is parameterised with
  AmberTools, and tleap's `addIonsRand` takes ion counts, not a concentration. The toolbox therefore builds each
  system once to count the waters, and adds round(concentration × waters / 55.5) Na⁺/Cl⁻ pairs on top of the
  ions that neutralise it. Set the concentration as you would for `genion`; the counts used are in the stage report.
- the box, in every stage that makes one: worked out from the structure by default (`"auto"`), or your own number.
  - `[relax] box_d`: nm between the protein and the box edge. Automatic: the protein's extent plus
    `image_margin_nm` (3.0), with any segment listed in `flexible` allowed to stretch out by its `reach_nm`.
  - `[linked] box`: buffer in Å. Automatic: each attached residue may point its full reach (`reach`: `"median"`,
    `"p95"`, `"max"` or nm; measured in `ncaa-relax`) in any direction, plus `image_margin_nm` (2.0).
  - `[free] box`: buffer in Å. Automatic: the protein stays one residue's reach plus `image_margin_nm` (2.0) from
    its own image, so that no free molecule can touch both.

  The trade-off, which each of these stages also prints with an atom estimate for your system:

  - **The automatic box probably removes periodic-image artefacts, but costs computer time.** The run time per
    nanosecond grows roughly with the number of atoms, and a box sized for a long attached residue can hold
    several times the atoms of one with a small fixed margin.
  - **A small fixed margin, such as 12 Å, leads to periodic-image problems when a flexible region sits near the box
    edge.** An attached residue that reaches several nanometres, or a disordered tail, crosses 12 Å of water easily
    and meets the image of the structure on the other side.

  After the build, each stage reports the box as built (x, y, z in nm), the water between your construct and the box
  edge on each side, and the distance from the construct to its own periodic image, whichever way the box was set.

  The reach used by the automatic box is measured in `ncaa-relax`, which is why that run is 10 ns by default. Even
  so it is not a fixed number: for a long flexible residue it can differ by a few tenths of a nanometre between two
  runs, so the box differs by a few ångström per side from one run to the next. For a box that is identical every
  time, `reach = "extended"` uses the residue's fully stretched length, worked out from its bonds with no
  simulation. That gives a much larger box.
  If you set your own box in `[linked]` and `[free]`, the reach is not used at all, so shorten the residue's run as
  well (`[ncaa.relax] production_ns`): it then only has to be long enough for the frames it exports.

  A number is the right choice when you know how far your structure can reach and what distance to its image you
  need. Either way, the periodic-image check after production says whether the box held.

### Your own mdp files

The config holds the settings most runs change: run lengths, temperature, restraints, box and salt. Everything
else GROMACS is told (cutoffs, thermostat and barostat, time step, constraints, how often output is written) is in
thirteen mdp files in `ncaaprotein/mdp/`, one per simulation: `protein_em`, `protein_nvt`, `protein_npt`,
`protein_prod`, `ncaa_em`, `ncaa_eq`, `ncaa_prod`, `linked_em`, `linked_nvt`, `linked_npt`, `linked_prod`,
`free_npt` and `free_prod` (`free-run` minimises with `linked_em` and runs NVT with `linked_nvt`). To change any
of them, copy them out, edit the copies and name the folder in your config:

```bash
python -m ncaaprotein copy-mdp my_mdp      # copies the thirteen files; never overwrites one already there
```

```toml
[project]
mdp_dir = "my_mdp"                           # relative to the config file
```

A stage then reads your file where your folder has one of that name and the packaged file where it does not, so
you can delete the ones you have not changed. Each stage prints which of your files it used. Keep the file names.

The toolbox still writes these keys itself, whatever your file says, because they come from the config or from
the system being built:

| Key | Set from |
|---|---|
| `nsteps` | `production_ns` always, and `nvt_ps` / `npt_ps` when the config gives them, divided by the `dt` in your file |
| `define` | which chains are restrained (`restrain`) |
| `tc-grps`, `compressed-x-grps` | the solute group of the linked and free systems |
| `ref_t`, `gen_temp` | `temperature`, only when the config gives one |

The checks after each simulation read the energy file, so an edit that writes energies less often gives them
fewer points to judge. `gmx grompp` refuses an inconsistent file, and the stage then stops and prints what GROMACS
objected to.

### Using more cores and GPUs

The simulations already use the whole machine: `gmx mdrun` takes the GPU and every CPU thread by default. Everything
else the toolbox runs uses one CPU core at a time. With more cores or GPUs, two things shorten a run:

- **The image-distance check in steps 3 and 20** (`gmx mindist -pi`) is the slowest of those single-core steps; on a
  long run it takes longer than the rest of its stage. Split it across cores under `[tools]`:

  ```toml
  [tools]
  cpu_jobs = 6        # the trajectory is cut into 6 stretches of time, measured side by side; the result is the same
  ```

- **Stages that do not depend on each other can run at the same time:** `protein-run` beside `ncaa-param ncaa-relax`
  (as in the quick start), and Run 1 beside Run 2 once both stages before them have finished. With two GPUs, give
  each run its own:

  ```bash
  CUDA_VISIBLE_DEVICES=0 nohup python -m ncaaprotein linked-build linked-run linked-analyse config.toml > run1.out 2>&1 &
  CUDA_VISIBLE_DEVICES=1 nohup python -m ncaaprotein free-build free-run free-analyse config.toml > run2.out 2>&1 &
  ```

  On one GPU the two runs share it, and each goes at roughly half speed.

### Checks that hold a run

The run stages (`ncaa-relax`, `protein-run`, `linked-run`, `free-run`) apply fixed criteria after
minimisation, NVT and NPT (`check_em.py`, `check_nvt.py`, `check_npt.py`). If every line passes, the stage goes on to
the next run by itself.

If an NVT or NPT check fails, the stage first gives the run more time: it continues that run for its own
length again and checks afresh, up to `extend_max` times (default 2; 0 switches this off). The checks look at the
last half of the run, so each time they judge new data. This matters because the checks are statistical: a run that is
equilibrated can still fail one line by a small margin purely by chance. A chance failure passes
on the next window; a system that is really not equilibrated keeps failing.

If a check still fails, or the minimisation ends "NOT converged", the stage stops before the next run, prints the
failed lines, and tells you where the plot is. Look at it, then either:

- add `accept = ["npt"]` (or `"em"`, `"nvt"`) to that stage's table in the config and run the same command again:
  it carries on from where it stopped; or
- fix the cause and repeat the step.

The temperature of each thermostat group must sit within 1 K of `ref_t` (or 3 standard errors, if more). When a
temperature line fails, the check says whether it is the solute or the water that is off. With a time step above
2 fs and repartitioned hydrogen masses (with the packaged mdp files: `free-run`'s NPT), the solute can read a few K
warm while the water stays at `ref_t`. This is a known effect of large time steps (Eastwood et al., *J. Chem. Theory
Comput.* 2010, 6, 2045; Jung, Kobayashi & Sugita, *J. Chem. Theory Comput.* 2019, 15, 84), and the check then says
to run that step at 2 fs (the stage also warns before such a run starts): in your own copy of the mdp file it names, set `dt = 0.002` and
`mass-repartition-factor = 1` (see "Your own mdp files" below). A warm solute and water at `ref_t` in a 2 fs run, or
water off `ref_t` in any run, is not this effect: look at the run. Production runs at 4 fs show the same offset, about
4 K for the solute; if that matters for what you measure, run production at 2 fs too.

`stop_before = "production"` makes a stage always stop before production, so the equilibration can be looked at
first. Minimisation is judged by its verdict line: "converged" and "usable" go on, "NOT converged" stops. In
`ncaa-relax` that is after a second pass, which is tried once if the first does not converge. `free-run` is the
exception: minimisation is expected to end unconverged on a system that size, and the NVT that follows is the judge.

### Analysis and plotting are separate scripts

Every analysis writes its numbers to CSV files, and a second script draws the figure from those files. To change
how a figure looks, edit only the plotting script, or write your own against the CSV: the numbers, verdicts and
reports do not depend on it. A plotting script can be run again by hand in the folder that holds the CSV files.

| Analysis (numbers, verdicts, CSV) | Writes | Plotting (figure only) |
|---|---|---|
| `check_em.py` | `em.csv` | `plot_em.py` → `em.png` |
| `check_nvt.py` | `nvt.csv`, `nvt_meta.json` | `plot_nvt.py` → `nvt.png` |
| `check_npt.py` | `npt.csv`, `npt_meta.json` | `plot_npt.py` → `npt.png` |
| `check_assembly_stability.py` | `assembly_stability.csv`, `assembly_stability_meta.json` | `plot_assembly_stability.py` → `assembly_stability.png` |
| `linked_eq_check.py` | `linked_eq.csv`, `linked_eq_meta.json` | `plot_linked_eq.py` → `linked_eq.png` |
| `free_ncaa_analysis.py` | four `free_ncaa_*.csv`, `free_ncaa_meta.json` | `plot_free_ncaa.py` → `free_ncaa_overview.png` |

The small `*_meta.json` beside a CSV holds what is not a table: the reference values, the averaging window, labels
and the verdicts the figure quotes. A figure that fails to draw is reported and does not stop a stage. In the
config, `plot_script` under `[relax.stability]`, `[linked.analysis]` or `[free.analysis]` names your own plotting
script.

### What a stage leaves behind

Each stage writes into its own folder under `workdir`:

- its files (`<name>.prepin`, `protein_ncaa_noh.pdb`, `linked.amb2gmx/prod.xtc` …);
- `REPORT.md`: the checks and their results, step by step;
- `commands.log`: every command run, in order, with the settings handed to each script;
- `logs/`: the full output of each command;
- `steps_done.json`: finished steps. Run the stage again and it carries on after the last finished step; a
  production run that stopped is continued from its checkpoint.

A stage never deletes a file. Each `mdrun` starts in its own session with its process ID in `<name>.pid`, so it
outlives the terminal, and it can be stopped on its own with `kill $(cat prod.pid)`.

### Scope and limits

- The toolbox does not build the protein assembly: bring a protein-only PDB of it.
- `linked-analyse` and `free-analyse` are general: they follow any attached or free residue from its anchor (the
  last atom it shares with its parent amino acid) to its tip, which you can name with `[linked.analysis] tip`;
  otherwise the atom most bonds from the anchor is taken. Their time windows scale with the length of the run, and
  both say so when a run is too short to judge.
- `protein-run` always runs a general stability check that needs no settings: it finds the chains from the structure
  and reports whether the Cα RMSD and the radius of gyration have settled and whether every pair of neighbouring
  chains kept its contacts. A check written for your own system can be run beside it (`[relax.stability]`), and hydrogen bonds between chosen chains are counted only if you ask
  (`[relax] hbond`).
- The NCAA must attach through a side chain and derive from a canonical residue. Proline-like and
  Cα-disubstituted residues are not handled.

## Data and large files
This repository holds source code and small input files. Trajectories are not stored here.

## Contributing
Pull requests are welcome. Every change is reviewed by the maintainer, Harvey Fung, who decides whether it is
merged.

To check a change to the parameter stages, run `ncaa-param` on a molecule of your own before and after it and
compare the residue's `.prepin` and `frcmod` files in `work/ncaa/`. It needs no GPU and takes about a minute.

## Authors and acknowledgment

**Authors**

- Harvey Fung (Team Imperial College, iGEM 2026)
- Claude Code (Anthropic's AI coding assistant), which drafted the code in `ncaaprotein/`

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

The toolbox also uses AmberTools (antechamber, parmchk2, tleap), acpype, the ff19SB,
GAFF2 and OPC models, PyMOL and MDAnalysis; these print no reference of their own, so cite them as their
documentation asks.

## License
This repository is licensed under the [Apache License 2.0](LICENSE). You are free to use, modify and distribute
this software, provided you keep the license and attribution notices.
