# Designing the capped HET-s fibril

These six steps produce `hets11_capped.pdb`, the protein-only structure the toolbox starts from (`[relax] structure`
in `../../config.toml`). Designing the assembly is specific to each project, so it is done by hand in PyMOL and
AlphaFold and is not part of the toolbox. This folder shows how it was done for the worked example.

| Step | What | Script | In | Out |
|---|---|---|---|---|
| 1 | Template structure: 2KJ3, model 1 | PyMOL commands below | `2KJ3.pdb` | object `m1` |
| 2 | Mutate residue 239 to Cys in all three chains | `mutate_cys239.py` | `m1` | `hets_3mer_mut.pdb` |
| 3 | Solvent accessibility of Cys239 | `sasa_cys.py` | mutated `m1` | side-chain SASA per chain |
| 4 | Build the 11-mer | `build_nmer.py` | `hets_3mer_mut.pdb` | `hets11.pdb` |
| 5 | Fold the cap construct, check the fold | AlphaFold Server, `check_cap_fold.py` | cap sequence, `hets11.pdb` | five `*_model_N.cif`, fold RMSD |
| 6 | Dock the cap on the fibril end | `superimpose_cap.py` | `hets11.pdb`, the five models | `hets11_capped.pdb` |

Run each in PyMOL with `run`, from a
folder that holds the script and its input files. Never paste a script into PyMOL: pasted loops run twice.

## 1 · Template structure

```
wget https://files.rcsb.org/download/2KJ3.pdb
```

In PyMOL:

```
load 2KJ3.pdb
split_states 2KJ3, 1, 1
delete 2KJ3
set_name 2KJ3_0001, m1
```

Model 1 is kept.

**The N-terminal methionine is kept.** 2KJ3 was deposited as Met217 + HET-s(218–289) + His tag. We were not sure
whether that methionine is removed when the protein is made in *E. coli*, so we kept Met217 in every chain; the His
tag is removed in step 4. Keeping it has a cost in the simulations. With 4 fs steps and repartitioned hydrogen
masses, the S–CH₃ at this floppy chain end turns far enough in one step to raise LINCS warnings, and `gmx mdrun`
stops after 1,000 of them. The example therefore runs with `GMX_MAXCONSTRWARN=-1` (set as `mdrun_env` in
`../../config.toml`), which lifts that limit. If your construct has no such methionine, delete it here instead
(`remove m1 and resi 217`) and leave `mdrun_env` out.

## 2 · Mutate residue 239 to Cys

```
run mutate_cys239.py
```

It ends with `wrote hets_3mer_mut.pdb: CYS239 in chains ABC`. On 2KJ3 model 1 rotamer 1 had the lowest strain in each
chain (27.6, 29.3 and 28.7).

## 3 · Solvent accessibility of Cys239

```
run sasa_cys.py
```

Expected: chain A exposed, chains B and C buried (A 21.5, B 0.0, C 0.0 Å²). Only the chain at the free end of the
fibril has an exposed Cys239, which is why the linker goes on A239.

## 4 · Build the 11-mer

```
run build_nmer.py
```

Settings at the top of the script: `SRC = "hets_3mer_mut.pdb"`, `OUT = "hets11.pdb"`, `N_CHAINS = 11`,
`STRIP = "resi 290-295"` (the His tag). It saves only if no interface clashes; expect 11 chains and 803 Cα (73 residues per chain, Met217 to Asn289).

## 5 · Fold the cap construct

Fold this sequence with AlphaFold Server and keep the five `.cif` models next to the scripts:

```
MKIDAIVGRNSAKDIRTEERARVQLGNVVTAAALHGGIRISDQTTNSVETVVGKGESRVLIGNEYGGKGFWDNGGSGGSGGSGGSGGSKIDAPVGRPSSKDIRTEESSRVQLGNVVTSSSLHHHHHHGSC
```

Then check that the predicted fold is the fibril's:

```
run check_cap_fold.py
```

Below about 1.5 Å the cap can stack; here the models fitted to 1.06 Å. The cap models are not in this
repository.

## 6 · Dock the cap

```
run superimpose_cap.py
```

Expected: model 3 fails the screen, model 0 is chosen, residues 1 to 7 are deleted and the cap is shifted 1.1 Å. The
result is `hets11_capped.pdb`: chains A to K and the cap as chain L, 12,933 atoms.

## The result

The result, `hets11_capped.pdb`, was copied into `protein_relax/`, where step 1 reads it. From here the toolbox
takes over:

```
cd ../..
python -m ncaaprotein protein-prepare protein-run config.toml
```

## References

The template structure is PDB entry **2KJ3**; the PDB asks that an entry be cited by its ID and its publication:

1. H. Van Melckebeke, C. Wasmer, A. Lange, E. AB, A. Loquet, A. Böckmann, B. H. Meier. Atomic-Resolution Three-Dimensional Structure of HET-s(218–289) Amyloid Fibrils by Solid-State NMR Spectroscopy. *J. Am. Chem. Soc.* 132, 13765–13775 (2010). doi:10.1021/ja104213j

The cap was folded with AlphaFold Server, whose developers ask that work using its outputs cite:

2. J. Abramson, J. Adler, J. Dunger, et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* 630, 493–500 (2024). doi:10.1038/s41586-024-07487-w

The design steps run in PyMOL: see the toolbox's [README](../../../../README.md#citing-what-the-toolbox-uses) for how to cite it.
