# Designing the capped HET-s fibril

These six steps produce `hets11_capped.pdb`, the protein-only structure the toolbox starts from (`[relax] structure`
in `../config.toml`). Designing the assembly is specific to each project, so it is done by hand in PyMOL and
AlphaFold and is not part of the toolbox. This folder shows how it was done for the worked example.

| Step | What | Script | In | Out |
|---|---|---|---|---|
| 1 | Template structure: 2KJ3, model 1, Met217 deleted | PyMOL commands below | `2KJ3.pdb` | object `m1` |
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
remove m1 and resi 217
```

Model 1 is kept and the cloning methionine (217) is removed.

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
`STRIP = "resi 290-295"` (the His tag). It saves only if no interface clashes; expect 11 chains and 792 Cα.

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
result is `hets11_capped.pdb`: chains A to K and the cap as chain L.

## The result

`hets11_capped.pdb` in this folder is the result. From here the toolbox takes over:

```
cd ..
python -m fibriltoolbox fibril-prepare fibril-run config.toml
```
