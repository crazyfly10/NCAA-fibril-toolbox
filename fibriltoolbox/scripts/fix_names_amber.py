"""
§19: GROMACS -> AMBER names, strip hydrogens, write a PDB tleap reads correctly.
Run in PyMOL:  run fix_names_amber.py
"""
from pymol import cmd

IN  = "full_fib.pdb"    # from place_linker.py
OUT = "full_fib_noh.pdb"

cmd.delete("all")
cmd.load(IN, "f")
cmd.remove("f and (solvent or inorganic)")    # water and ions, if any: protein + CYL only
cmd.remove("f and hydro")      # all H, CYL's too: tleap rebuilds them from the
                               # templates (amino19.lib, CYL.prepin)
# GROMACS terminal residue names are 4 characters (NMET, NLYS, CASN, NGLY, CCYS) and tleap
# reads only 3 -- NMET would become NME, the N-methyl cap. Use the plain names:
# tleap makes the first/last residue of each chain N-/C-terminal by itself.
for old, new in {"NMET": "MET", "NME": "MET", "NLYS": "LYS", "NLY": "LYS",
                 "CASN": "ASN", "CAS": "ASN", "NGLY": "GLY", "NGL": "GLY",
                 "CCYS": "CYS", "CCY": "CYS"}.items():
    cmd.alter(f"f and resn {old}", f"resn='{new}'")
cmd.alter("f and resn ILE and name CD", "name='CD1'")
cmd.alter("f and name OC1", "name='O'")
cmd.alter("f and name OC2", "name='OXT'")
cmd.sort("f")

n_cyl = cmd.count_atoms("f and resn CYL")
print(f"CYL heavy atoms: {n_cyl}  (expected 72 per CYL, 144 for two)")

# PyMOL's own TER records also go where a standard residue is followed by a
# non-standard one (CYL), and tleap reads every TER as a chain break -- CYL would
# lose its peptide bond to the residue before it. So: no PyMOL TERs, no CONECT;
# write a TER only where the chain ID changes.
cmd.set("pdb_use_ter_records", 0)
cmd.save(OUT, "f")

out, prev, chains = [], None, []
for line in open(OUT):
    if line.startswith(("ATOM", "HETATM")):
        ch = line[21]
        if prev is not None and ch != prev:
            out.append("TER\n")
        if ch != prev:
            chains.append(ch)
        prev = ch
        out.append(line)
    elif line.startswith(("TER", "END", "CONECT", "MASTER")):
        continue
    else:
        out.append(line)
out += ["TER\n", "END\n"]
open(OUT, "w").writelines(out)
print(f"wrote {OUT}: chains {''.join(chains)} ({len(chains)} chains, {len(chains)} TER records)")
