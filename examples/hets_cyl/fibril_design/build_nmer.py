"""
build_nmer.py — tessellate a mutated 2KJ3 trimer into an N-mer fibril.

Run from the PyMOL command line:
    run build_nmer.py

Do NOT paste this file into the prompt. Pasted multi-line blocks nest the
interpreter and silently re-execute loops — one such paste produced 66
monomers instead of 12.

WHY IT IS BUILT THIS WAY
  1. The screw transform is EXTRACTED from the structure, not assumed.
  2. It is fitted on the ORDERED CORE only. About 30 of 72 residues are
     disordered (the N-terminal 218-225, the 247-261 loop) and a whole-chain
     superposition is dominated by them — fitting all atoms gave 4.50 and
     6.27 A for two pairs that should agree.
  3. The transform comes from an INTERIOR pair (B->C). Chain A is terminal
     and has relaxed differently.
  4. An INTERIOR monomer is replicated. A terminal monomer's side chains have
     relaxed into space a neighbour would occupy; replicating chain C gave
     0.96 A clashes at every new interface.
  5. Transforms are COMPOSED (Rk, tk), not applied by chaining copies, so
     numerical error does not accumulate.
"""

from pymol import cmd
import numpy as np
import os

# ------------------------------------------------------------------ config
SRC      = "hets_3mer_mut.pdb"     # mutated 2KJ3 trimer
OUT      = "hets11.pdb"
N_CHAINS = 11                      # 11 for the capped construct, 12 for uncapped
LETTERS  = "ABCDEFGHIJKLMNOPQRSTUVWX"

# four beta-strands of the solenoid core, in 2KJ3 author numbering
CORE = "resi 226-234+237-245+262-270+273-282 and name CA"

# set to a range to delete before replicating, e.g. "resi 290-295" for a
# His tag. Replication CANNOT represent a disordered region: twelve identical
# copies of one guessed conformation all point the same way and collide.
STRIP = "resi 290-295"


def kabsch(P, Q):
    """rotation + translation mapping point set P onto Q"""
    cP, cQ = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((P - cP).T @ (Q - cQ))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, cQ - R @ cP


def interfaces(obj, label=""):
    """minimum heavy-atom distance across every consecutive chain pair"""
    ch = cmd.get_chains(obj)
    print(f"\n  {label}")
    bad = 0
    for a, b in zip(ch, ch[1:]):
        P = cmd.get_coords(f"{obj} and chain {a} and not hydro")
        Q = cmd.get_coords(f"{obj} and chain {b} and not hydro")
        D = np.linalg.norm(P[:, None, :] - Q[None, :, :], axis=-1)
        n25, n35 = int((D < 2.5).sum()), int((D < 3.5).sum())
        flag = "  <-- CLASH" if D.min() < 2.0 else ""
        if D.min() < 2.0:
            bad += 1
        print(f"    {a}-{b}: min {D.min():.2f}  n<2.5 {n25:4d}  n<3.5 {n35:4d}{flag}")
    return bad


# -------------------------------------------------------------------- build
cmd.reinitialize()

assert os.path.exists(SRC), f"not found: {SRC}"
cmd.load(SRC, "src")

if cmd.count_states("src") > 1:
    cmd.split_states("src", 1, 1)
    cmd.delete("src")
    cmd.set_name("src_0001", "src")

if STRIP:
    n_before = cmd.count_atoms("src")
    cmd.remove("src and " + STRIP)
    n_after = cmd.count_atoms("src")
    print(f"stripped {STRIP}: removed {n_before - n_after} atoms "
          f"({n_before} -> {n_after} left)")

ch = cmd.get_chains("src")
print(f"\nsource: {len(ch)} chains {ch}, "
      f"{cmd.count_atoms('src and name CA')} CA, "
      f"{cmd.count_atoms('src and name CA') // max(1, len(ch))} per chain")
assert len(ch) >= 3, "need at least 3 chains to take an interior transform"

# verify the mutation landed where intended
for c in ch:
    cys = [int(a.resi) for a in
           cmd.get_model(f"src and chain {c} and resn CYS and name CA").atom]
    print(f"  chain {c}: CYS at {sorted(cys)}")

# transform from the interior pair B -> C
B = cmd.get_coords("src and chain B and " + CORE)
C = cmd.get_coords("src and chain C and " + CORE)
assert B.shape == C.shape, f"core selection mismatch: {B.shape} vs {C.shape}"
assert len(B) > 20, f"core selection only {len(B)} atoms — check CORE numbering"
R, t = kabsch(B, C)

rise = np.linalg.norm(C.mean(0) - B.mean(0))
ang  = np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1)))
print(f"\ncore atoms per chain: {len(B)}")
print(f"rise  {rise:.2f} A   twist {ang:.2f} deg")
assert 8.5 < rise < 10.5, "rise outside the expected 9-10 A — check CORE"

# A stays as the real terminus; B is the repeat unit
cmd.create("m00", "src and chain A")
cmd.create("m01", "src and chain B")
parts = ["m00", "m01"]

Rk, tk = np.eye(3), np.zeros(3)
for i in range(2, N_CHAINS):
    Rk, tk = R @ Rk, R @ tk + t
    nm = "m%02d" % i
    cmd.create(nm, "src and chain B")
    cmd.load_coords(cmd.get_coords(nm) @ Rk.T + tk, nm)
    cmd.alter(nm, "chain='%s'" % LETTERS[i])
    cmd.alter(nm, "segi='%s'" % LETTERS[i])
    cmd.sort(nm)
    parts.append(nm)

cmd.create("nmer", " or ".join(parts))
for p in parts:
    cmd.delete(p)

# ------------------------------------------------------------------ verify
out_ch = cmd.get_chains("nmer")
print(f"\nbuilt: {len(out_ch)} chains {out_ch}, "
      f"{cmd.count_atoms('nmer and name CA')} CA")

print("\n  axial spacing of core centroids")
prev = None
for c in out_ch:
    com = cmd.get_coords(f"nmer and chain {c} and {CORE}").mean(0)
    if prev is not None:
        print(f"    {c}: {np.linalg.norm(com - prev):.2f} A")
    prev = com

bad = interfaces("nmer", "inter-chain contacts")

print(f"""
  GATE
    chains            {len(out_ch)}  (expected {N_CHAINS})
    spacing           every step within ~0.1 A of {rise:.2f}
    interfaces        min ~2.6 A, n<2.5 == 0, n<3.5 ~120
    interfaces failed {bad}

  Contacts at 2.0-2.5 A may be legitimate structure, not clashes. Identify
  each one before removing anything — in this system they were a Leu/His
  packing contact, an Asp-Arg SALT BRIDGE and an Asn-Asn LADDER: interactions
  worth keeping.
""")

if bad:
    print("  !! not saving — resolve the clashes first")
else:
    cmd.set("pdb_use_ter_records", 1)   # or downstream tools read one chain
    cmd.save(OUT, "nmer")
    print(f"  saved {OUT}")
