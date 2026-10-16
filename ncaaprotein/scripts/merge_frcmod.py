#!/usr/bin/env python3
"""
merge_frcmod.py  --  build ONE clean frcmod for a non-standard residue in ff19SB
whose attached part uses GAFF2 atom types (the part shared with the parent residue keeps ff19SB types).

Inputs (both from parmchk2 on the same molecule):
  --primary   parmchk2 against GAFF2 (-s 2). Lists everything gaff2.dat lacks:
              the ff19SB backbone terms (already in ff19SB, dropped here), the
              junction terms between ff19SB and GAFF2 types, and GAFF2 gaps.
  --fallback  parmchk2 against ff19SB (parm19.dat + frcmod.ff19SB); only used to
              fill an entry the primary run marks "ATTN, need revision".

Output (-o): only terms that no loaded force field provides. Dropped:
  * any bond/angle/dihedral/mass/nonbond already in parm19.dat, frcmod.ff19SB or gaff2.dat
  * any dihedral made only of ff19SB types that tleap would already match with a generic X-B-C-X term
    in those files: a specific line would change that torsion for EVERY residue in the protein. A dihedral
    with a GAFF2 type (lower case) occurs only in the residue, so parmchk2's analogue ("same as ...") is
    kept: it is closer to the chemistry than the generic term tleap would otherwise fall back to
  * any improper tleap would already match (central atom third, outer atoms sorted,
    wildcards). parmchk2 writes some as specific lines (e.g. N -O -C -XC with a 1.1
    default instead of X -X -C -O 10.5); a specific line beats the generic one for
    EVERY residue in the protein.
Because nothing in the output overlaps ff19SB or GAFF2, tleap load order does not matter.

Usage
  python3 merge_frcmod.py --primary frcmod.NCA.gaff2 --fallback frcmod.NCA.ff19SB -o frcmod.NCA
"""
import argparse
import os
import re
import sys

SECTIONS = ("MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON")
NTYPES = {"BOND": 2, "ANGLE": 3, "DIHE": 4, "IMPROPER": 4}
PEN = re.compile(r"penalty score=\s*([-\d.]+)")


# ----------------------------------------------------------------- keys
def types_of(line, sec):
    if sec in ("MASS", "NONBON"):
        return (line.split()[0],)
    n = NTYPES[sec]
    return tuple(t.strip() for t in line[:3 * n - 1].split("-"))


def canon(sec, t):
    """Key under which tleap treats two entries as the same parameter."""
    if sec in ("BOND", "ANGLE", "DIHE"):
        return min(t, t[::-1])
    if sec == "IMPROPER":                      # central atom 3rd, others sorted, X first
        outer = sorted((t[0], t[1], t[3]), key=lambda s: (s != "X", s))
        return (outer[0], outer[1], t[2], outer[2])
    return t


def tleap_improper_match(ref, cand):
    """Port of zbParmSetMatchImproper (leap/parmSet.c); both keys canonical."""
    if ref[2] != cand[2]:
        return False
    s1, s2, s4 = cand[0], cand[1], cand[3]
    nwild = sum(1 for x in (ref[0], ref[1], ref[3]) if x == "X")
    if nwild == 0:
        return (ref[0], ref[1], ref[3]) == (s1, s2, s4)
    if nwild == 1:
        if ref[1] == s1:
            return ref[3] in (s2, s4)
        return ref[1] == s2 and ref[3] == s4
    if nwild == 2:
        return ref[3] in (s1, s2, s4)
    return True


# ------------------------------------------------------------- parsers
def read_frcmod(path):
    """Return {section: [entry]}; entry = dict(key, types, lines, text)."""
    out = {s: [] for s in SECTIONS}
    sec = None
    with open(path) as fh:
        lines = fh.read().splitlines()
    for line in lines[1:]:                       # first line is a title
        head = line.strip().upper()
        if head[:4] in ("MASS", "BOND", "ANGL", "DIHE", "IMPR", "NONB", "CMAP"):
            sec = {"ANGL": "ANGLE", "IMPR": "IMPROPER", "NONB": "NONBON"}.get(head[:4], head[:4])
            if sec == "CMAP":
                break
            continue
        if not line.strip() or sec is None:
            continue
        t = types_of(line, sec)
        key = canon(sec, t)
        # multi-term dihedrals: consecutive lines with the same key belong together
        if sec == "DIHE" and out[sec] and out[sec][-1]["key"] == key:
            out[sec][-1]["lines"].append(line)
            out[sec][-1]["text"] += " " + line
            continue
        out[sec].append({"key": key, "types": t, "lines": [line], "text": line})
    return out


def read_parm_dat(path):
    """Keys defined in a main Amber parameter file (parm19.dat / gaff2.dat style)."""
    keys = {s: set() for s in SECTIONS}
    order = ["MASS", "HYDRO", "BOND", "ANGLE", "DIHE", "IMPROPER", "HBOND", "EQUIV"]
    state = 0
    equiv = []
    with open(path) as fh:
        lines = fh.read().splitlines()[1:]
    in_nonb = False
    for line in lines:
        if line.startswith("END"):
            break
        if line.startswith("MOD4"):
            in_nonb = True
            continue
        if in_nonb:
            if line.strip():
                keys["NONBON"].add((line.split()[0],))
            continue
        name = order[state] if state < len(order) else None
        if name == "HYDRO":                      # single line, no blank after it
            state += 1
            continue
        if not line.strip():
            state += 1
            continue
        if name == "MASS":
            keys["MASS"].add((line.split()[0],))
        elif name in NTYPES:
            keys[name].add(canon(name, types_of(line, name)))
        elif name == "EQUIV":
            equiv.append(line.split())
    for grp in equiv:                            # e.g. "N NA N2 ..." share N's LJ
        for t in grp[1:]:
            keys["NONBON"].add((t,))
    return keys


def reference_keys(paths):
    ref = {s: set() for s in SECTIONS}
    for p in paths:
        with open(p) as fh:
            fh.readline()
            second = fh.readline()
        is_frcmod = any(second.upper().startswith(h) for h in ("MASS", "BOND", "ANGL", "DIHE",
                                                                "IMPR", "NONB"))
        if is_frcmod:
            for sec, ents in read_frcmod(p).items():
                ref[sec].update(e["key"] for e in ents)
        else:
            for sec, ks in read_parm_dat(p).items():
                ref[sec].update(ks)
    return ref


def already_defined(sec, e, ref, ref_impr):
    """Return the reference term tleap would use for this entry, or None."""
    if sec == "IMPROPER":
        hit = [r for r in ref_impr if tleap_improper_match(r, e["key"])]
        return "-".join(hit[0]) if hit else None
    if e["key"] in ref[sec]:
        return "-".join(e["key"])
    if sec == "DIHE" and not any(x[:1].islower() for x in e["types"]):
        t = e["types"]                           # ff19SB-only torsion: tleap falls back to X-B-C-X
        gen = canon("DIHE", ("X", t[1], t[2], "X"))
        if gen in ref[sec]:
            return "-".join(gen)
    return None


def penalty(entry):
    m = PEN.findall(entry["text"])
    return max(float(x) for x in m) if m else 0.0


# ----------------------------------------------------------------- main
def main():
    home = os.environ.get("AMBERHOME", "")
    dref = [os.path.join(home, "dat/leap/parm", f)
            for f in ("parm19.dat", "frcmod.ff19SB", "gaff2.dat")]
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--primary", required=True, help="parmchk2 output made against GAFF2 (-s 2)")
    ap.add_argument("--fallback", help="parmchk2 output made against ff19SB, for ATTN entries")
    ap.add_argument("--ref", nargs="+", default=dref,
                    help="force-field files tleap loads (default parm19.dat frcmod.ff19SB gaff2.dat)")
    ap.add_argument("--flag-above", type=float, default=10.0,
                    help="list entries whose parmchk2 penalty is above this for a manual check")
    ap.add_argument("-o", "--output", required=True)
    a = ap.parse_args()

    for p in a.ref:
        if not os.path.isfile(p):
            sys.exit("ERROR: reference file %s not found (is $AMBERHOME set?)" % p)
    ref = reference_keys(a.ref)
    ref_impr = list(ref["IMPROPER"])
    prim = read_frcmod(a.primary)
    fall = read_frcmod(a.fallback) if a.fallback else {s: [] for s in SECTIONS}
    fidx = {s: {e["key"]: e for e in fall[s]} for s in SECTIONS}

    final = {s: [] for s in SECTIONS}
    log, dropped = [], 0
    for sec in SECTIONS:
        for e in prim[sec]:
            label = "%-8s %s" % (sec, "-".join(e["types"]))
            used = already_defined(sec, e, ref, ref_impr)
            if used:
                dropped += 1
                log.append("drop     %-30s tleap uses %s from the loaded force fields" % (label, used))
                continue
            if "ATTN" in e["text"]:
                g = fidx[sec].get(e["key"])
                if g is None or "ATTN" in g["text"]:
                    log.append("FIX!!    %-30s no estimate from either run -- set by hand" % label)
                    final[sec].append(e["lines"])
                    continue
                final[sec].append(g["lines"])
                log.append("ff19SB   %-30s GAFF2 had no analogue; ff19SB analogue used" % label)
                continue
            final[sec].append(e["lines"])
            pen = penalty(e)
            note = "  <- check (penalty %.1f)" % pen if pen > a.flag_above else ""
            log.append("keep     %-30s%s" % (label, note))

    with open(a.output, "w") as fh:
        fh.write("frcmod for residue terms missing from ff19SB + GAFF2 (merge_frcmod.py)\n")
        for sec in SECTIONS:
            fh.write("%s\n" % sec)
            for lines in final[sec]:
                for l in lines:
                    fh.write(l.rstrip() + "\n")
            fh.write("\n")
        fh.write("\n")

    print("\n".join(l for l in log if not l.startswith("drop")))
    print("(%d entries dropped: already provided by ff19SB or GAFF2)" % dropped)
    left = sum(1 for s in SECTIONS for ls in final[s] for l in ls if "need revision" in l)
    print("\nWrote %s: %d entries, %d still marked ATTN"
          % (a.output, sum(len(v) for v in final.values()), left))
    if left:
        sys.exit(1)


if __name__ == "__main__":
    main()
