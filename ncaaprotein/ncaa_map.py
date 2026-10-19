"""Find the amino-acid part of a capped non-canonical residue and put it on ff19SB names, types and charges.

This renames and retypes the amino-acid part (step 6) with no atom numbers typed in by hand. The capped
molecule, ACE-X-NME, is read from antechamber's .ac file as a graph of bonds, and three things are located in it:

  the two caps         CH3-C(=O)-  and  -NH-CH3      (removed later by prepgen)
  the backbone         N, H, CA, HA, C, O
  the parent's atoms   the side-chain atoms X shares with the canonical residue it derives from
                       (for instance CB and SG when the parent is CYX)

Those atoms get the parent's names and ff19SB types from amino19.lib (CA becomes XC, which ties the residue to
its CMAP). Backbone atoms and caps get ff19SB charges; the leftover is spread over the rest of the residue so it
sums to the residue's net charge. Everything further out keeps antechamber's names, GAFF2 types and the charges
step 5 computed (charge_method).
"""
import json
import os
from pathlib import Path

import networkx as nx

# name, ff19SB type, charge -- aminont12.lib (ACE) and aminoct12.lib (NME)
ACE = {"CAY": ("CT", -0.3662), "HY": ("HC", 0.1123), "CY": ("C", 0.5972), "OY": ("O", -0.5679)}
NME = {"NT": ("N", -0.4157), "HNT": ("H", 0.2719), "CAT": ("CT", -0.1490), "HT": ("H1", 0.0976)}
BACKBONE = ("N", "H", "CA", "HA", "C", "O")


class MapError(RuntimeError):
    pass


def element(gaff_type):
    t = gaff_type.strip()
    return {"cl": "Cl", "br": "Br"}.get(t[:2].lower(), t[:1].upper())


def read_ac(path):
    lines = open(path).read().splitlines()
    g = nx.Graph()
    for k, line in enumerate(lines):
        if line.startswith("ATOM"):
            n = int(line[4:11])
            name, q, t = line[13:17].strip(), float(line[54:64]), line[64:].strip()
            if (q, t) != (float(line.split()[-2]), line.split()[-1]):
                raise MapError(f"unexpected .ac column layout on line {k + 1}")
            g.add_node(n, line=k, name=name, q=q, type=t, el=element(t))
        elif line.startswith("BOND"):
            f = line.split()
            g.add_edge(int(f[2]), int(f[3]))
    return lines, g


def parent_template(parent):
    """Names, types, charges and bonds of a canonical residue, from ff19SB's amino19.lib."""
    from parmed.amber import AmberOFFLibrary
    lib = AmberOFFLibrary.parse(os.path.join(os.environ["AMBERHOME"], "dat", "leap", "lib", "amino19.lib"))
    if parent not in lib:
        raise MapError(f"{parent} is not a residue of amino19.lib")
    res = lib[parent]
    p = nx.Graph()
    for a in res.atoms:
        p.add_node(a.name, type=a.type, q=a.charge, el="H" if a.atomic_number == 1 else a.name[0])
    for b in res.bonds:
        p.add_edge(b.atom1.name, b.atom2.name)
    return p


def _heavy(g, n):
    return [m for m in g[n] if g.nodes[m]["el"] != "H"]


def _hyd(g, n):
    return sorted(m for m in g[n] if g.nodes[m]["el"] == "H")


def find_backbone(g):
    """The one place where CH3-C(=O)-NH-CH(R)-C(=O)-NH-CH3 fits the molecule."""
    el = lambda n: g.nodes[n]["el"]
    methyl = lambda n, frm: el(n) == "C" and _heavy(g, n) == [frm] and len(_hyd(g, n)) == 3
    oxo = lambda n: el(n) == "O" and len(g[n]) == 1
    hits = []
    for N in g:
        if el(N) != "N" or len(_hyd(g, N)) != 1 or len(_heavy(g, N)) != 2:
            continue
        for CY, CA in (_heavy(g, N), _heavy(g, N)[::-1]):
            if el(CY) != "C" or el(CA) != "C":
                continue
            rest = [m for m in _heavy(g, CY) if m != N]
            OY = [m for m in rest if oxo(m)]
            CAY = [m for m in rest if methyl(m, CY)]
            if len(rest) != 2 or len(OY) != 1 or len(CAY) != 1:
                continue
            for C in [m for m in _heavy(g, CA) if m != N and el(m) == "C"]:
                r2 = [m for m in _heavy(g, C) if m != CA]
                O = [m for m in r2 if oxo(m)]
                NT = [m for m in r2 if el(m) == "N" and len(_hyd(g, m)) == 1 and len(_heavy(g, m)) == 2]
                if len(r2) != 2 or len(O) != 1 or len(NT) != 1:
                    continue
                CAT = [m for m in _heavy(g, NT[0]) if m != C and methyl(m, NT[0])]
                if len(CAT) != 1:
                    continue
                hits.append(dict(CAY=CAY[0], CY=CY, OY=OY[0], N=N, CA=CA, C=C, O=O[0], NT=NT[0], CAT=CAT[0]))
    if len(hits) != 1:
        raise MapError(f"found {len(hits)} places where ACE-NH-CH(R)-CO-NME fits, expected exactly 1: "
                       "the molecule must be drawn with an acetyl cap on the N and an N-methyl cap on the C")
    return hits[0]


def match_side_chain(g, bb, parent, keep=None):
    """Map the parent's side-chain heavy atoms onto the molecule, walking out from CA. Returns {parent name: atom}."""
    side = [m for m in _heavy(g, bb["CA"]) if m not in (bb["N"], bb["C"])]
    if len(side) != 1:
        raise MapError(f"CA carries {len(side)} side-chain atoms; one is supported (no Gly-like or disubstituted CA)")
    p_heavy = lambda a: [b for b in parent[a] if parent.nodes[b]["el"] != "H" and b not in ("N", "C", "O")]
    p_side = [b for b in p_heavy("CA")]
    if len(p_side) != 1:
        raise MapError("the parent residue has no single CB")
    mapped = {p_side[0]: side[0]}
    if g.nodes[side[0]]["el"] != parent.nodes[p_side[0]]["el"]:
        raise MapError("the atom on CA is not the parent's CB element")
    frontier = [p_side[0]]
    while frontier:
        pa = frontier.pop(0)
        for pb in p_heavy(pa):
            if pb in mapped or pb == "CA" or (keep and pb not in keep):
                continue
            free = [m for m in _heavy(g, mapped[pa]) if m not in mapped.values() and m != bb["CA"]
                    and g.nodes[m]["el"] == parent.nodes[pb]["el"]]
            if len(free) == 1:
                mapped[pb] = free[0]
                frontier.append(pb)
            elif len(free) > 1:
                raise MapError(f"cannot tell which atom is the parent's {pb}: set [ncaa] parent_atoms to stop before it")
    if keep:
        mapped = {k: v for k, v in mapped.items() if k in keep}
    return mapped


def build_table(g, parent_name, keep=None, fixed=BACKBONE):
    """{atom number: (new name, ff19SB type, charge or None)}, the cap atoms, and a description of the residue."""
    bb = find_backbone(g)
    parent = parent_template(parent_name)
    side = match_side_chain(g, bb, parent, keep)
    table, caps = {}, []

    def put(n, name, typ, q):
        table[n] = (name, typ, q)

    put(bb["CAY"], "CAY", *ACE["CAY"])
    for i, h in enumerate(_hyd(g, bb["CAY"]), 1):
        put(h, f"HY{i}", *ACE["HY"])
    put(bb["CY"], "CY", *ACE["CY"])
    put(bb["OY"], "OY", *ACE["OY"])
    caps += [bb["CAY"]] + _hyd(g, bb["CAY"]) + [bb["CY"], bb["OY"]]

    for pname, n in [("N", bb["N"]), ("CA", bb["CA"])] + list(side.items()) + [("C", bb["C"]), ("O", bb["O"])]:
        pq = lambda a: round(parent.nodes[a]["q"], 4) if a in fixed else None
        put(n, pname, parent.nodes[pname]["type"], pq(pname))
        p_h = sorted(b for b in parent[pname] if parent.nodes[b]["el"] == "H")
        hs = _hyd(g, n)
        if len(hs) > len(p_h):
            raise MapError(f"{pname} carries {len(hs)} hydrogens, the parent's {pname} only {len(p_h)}")
        for h, hname in zip(hs, p_h[len(p_h) - len(hs):] if pname != "N" else p_h):
            put(h, hname, parent.nodes[hname]["type"], pq(hname))

    put(bb["NT"], "NT", *NME["NT"])
    for h in _hyd(g, bb["NT"]):
        put(h, "HNT", *NME["HNT"])
    put(bb["CAT"], "CAT", *NME["CAT"])
    for i, h in enumerate(_hyd(g, bb["CAT"]), 1):
        put(h, f"HT{i}", *NME["HT"])
    caps += [bb["NT"]] + _hyd(g, bb["NT"]) + [bb["CAT"]] + _hyd(g, bb["CAT"])

    # the first two side-chain torsions, for the rotamer scan when the residue is grafted
    def onward(prev, cur):
        nxt = [m for m in _heavy(g, cur) if m != prev]
        if not nxt:
            return None
        h = g.copy()
        h.remove_node(cur)
        return max(nxt, key=lambda m: len(nx.node_connected_component(h, m)))
    cb = side[[k for k in side][0]]
    x1 = onward(bb["CA"], cb)
    x2 = onward(cb, x1) if x1 else None
    chi = [[bb["N"], bb["CA"], cb, x1]] if x1 else []
    if x2:
        chi.append([bb["CA"], cb, x1, x2])
    info = dict(backbone=bb, side=side, chi=chi)
    return table, caps, info


def retype(ac_in, ac_out, parent, net_charge=0, set_charges=True, keep=None, map_out=None, resname="NCA"):
    """Write ac_out with the residue's amino-acid part renamed, retyped and recharged; return the lines to print."""
    lines, g = read_ac(ac_in)
    table, caps, info = build_table(g, parent, keep)
    capset = set(caps)
    atoms = {n: [g.nodes[n]["line"], g.nodes[n]["name"], g.nodes[n]["q"], g.nodes[n]["type"]] for n in sorted(g)}
    out = [" no.  old name type    ->  new name type   charge"]
    for n in sorted(table, key=lambda n: (n in capset, n)):
        new, t, q = table[n]
        k, old, q_old, t_old = atoms[n]
        if set_charges and q is not None:
            atoms[n][2] = q
        out.append("%4d  %-5s %-4s   ->  %-5s    %-4s  %8.4f%s" % (n, old, t_old, new, t, atoms[n][2],
                   "   (cap, removed by prepgen)" if n in capset else ""))
        atoms[n][1], atoms[n][3] = new, t

    cb = info["side"][list(info["side"])[0]]
    if set_charges:
        fixed = {n for n, (_, _, q) in table.items() if q is not None}
        side = [n for n in atoms if n not in fixed]
        leftover = net_charge - sum(atoms[n][2] for n in atoms if n not in capset)
        for n in side:
            atoms[n][2] = round(atoms[n][2] + leftover / len(side), 6)
        atoms[cb][2] = round(atoms[cb][2] - (sum(atoms[n][2] for n in atoms if n not in capset) - net_charge), 6)
        out.append("\nleftover %+.4f e spread over %d side-chain atoms (%+.5f e each)"
                   % (leftover, len(side), leftover / len(side)))

    names = [a[1] for a in atoms.values()]
    dup = sorted({x for x in names if names.count(x) > 1})
    if dup:
        raise MapError(f"duplicate atom names after renaming: {dup}")

    for n, (k, name, q, t) in atoms.items():
        l = lines[k]
        lines[k] = l[:13] + "%-4s" % name + l[17:54] + "%10.6f" % q + "%10s" % t
    for k, line in enumerate(lines):
        if line.startswith("BOND"):
            f = line.split()
            i, j = int(f[2]), int(f[3])
            lines[k] = "BOND%5d%5d%5d%5s  %5s%5s" % (int(f[1]), i, j, f[4], atoms[i][1], atoms[j][1])
    Path(ac_out).write_text("\n".join(lines) + "\n")

    q_res = sum(atoms[n][2] for n in atoms if n not in capset)
    q_cap = sum(atoms[n][2] for n in capset)
    out.append("residue charge %+.6f   caps %+.6f   ->  wrote %s" % (q_res, q_cap, Path(ac_out).name))

    nm = lambda n: atoms[n][1]
    heavy = [n for n in atoms if g.nodes[n]["el"] != "H"]
    desc = dict(
        resname=resname, parent=parent, net_charge=net_charge,
        n_atoms=len(atoms), n_residue_atoms=len(atoms) - len(capset),
        cap_names=[nm(n) for n in caps],
        heavy_atoms={"ACE": 3, resname: sum(1 for n in heavy if n not in capset), "NME": 2},
        backbone_names=["N", "CA", "C", "O"],
        parent_side_names=list(info["side"]),
        chi=[[nm(n) for n in d] for d in info["chi"]],
        parent_types={nm(n): table[n][1] for n in table if n not in capset},
    )
    if map_out:
        Path(map_out).write_text(json.dumps(desc, indent=1) + "\n")
    return out, desc


def write_mc(path, desc):
    """prepgen's main-chain file: head N, tail C, main chain CA, the cap atoms to omit."""
    text = "HEAD_NAME N\nTAIL_NAME C\nMAIN_CHAIN CA\n"
    text += "".join(f"OMIT_NAME {n}\n" for n in desc["cap_names"])
    text += f"PRE_HEAD_TYPE C\nPOST_TAIL_TYPE N\nCHARGE {float(desc['net_charge']):.1f}\n"
    Path(path).write_text(text)
