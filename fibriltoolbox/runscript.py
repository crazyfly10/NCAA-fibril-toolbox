"""Run a script with some of the settings at its top replaced, without editing the script.

The scripts keep their settings as plain assignments near the top (SRC = "x.pdb", N_CHAINS = 11 ...).
This wrapper parses the script, swaps the value of every top-level assignment named in the overrides, and runs
the result as if the script had been started directly. It works under plain Python and inside PyMOL:

    FTB_SCRIPT=build_nmer.py FTB_OVERRIDES='{"N_CHAINS": 7}' pymol -cq runscript.py
    FTB_SCRIPT=check_npt.py FTB_ARGS='["npt.xvg", "1025", "1055"]' python runscript.py

FTB_OVERRIDES is JSON. A string value that starts with "py:" is taken as a Python expression instead, for values
JSON cannot hold (tuples, dicts with number keys, a glob): "py:{226: 'ASN'}", "py:sorted(glob.glob('m/*.cif'))".
FTB_STATUS names a file that receives "ok" or the error, because PyMOL's exit code does not report a failed script.
"""
import ast
import json
import os
import sys
import traceback


def _value(v):
    if isinstance(v, str) and v.startswith("py:"):
        return ast.parse(v[3:], mode="eval").body
    return ast.parse(repr(v), mode="eval").body


def run(script, overrides, argv):
    tree = ast.parse(open(script).read(), script)
    used = set()
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        t = node.targets[0]
        if isinstance(t, ast.Name) and t.id in overrides:
            node.value = _value(overrides[t.id])
            used.add(t.id)
        elif isinstance(t, ast.Tuple) and isinstance(node.value, ast.Tuple) and len(t.elts) == len(node.value.elts):
            for i, e in enumerate(t.elts):
                if isinstance(e, ast.Name) and e.id in overrides:
                    node.value.elts[i] = _value(overrides[e.id])
                    used.add(e.id)
    missing = sorted(set(overrides) - used)
    if missing:
        raise SystemExit(f"{os.path.basename(script)} has no top-level setting named {', '.join(missing)}")
    ast.fix_missing_locations(tree)
    sys.argv = [script] + list(argv)
    exec(compile(tree, script, "exec"), {"__name__": "__main__", "__file__": script})


def main():
    script = os.environ["FTB_SCRIPT"]
    overrides = json.loads(os.environ.get("FTB_OVERRIDES", "{}"))
    argv = json.loads(os.environ.get("FTB_ARGS", "[]"))
    status = os.environ.get("FTB_STATUS")
    state = "ok"
    try:
        run(script, overrides, argv)
    except SystemExit as e:
        if e.code not in (None, 0):
            state = f"exit: {e.code}"
            print(e.code if isinstance(e.code, str) else f"exit code {e.code}")
    except BaseException as e:                      # report it, whatever it is; PyMOL would swallow it
        state = f"{type(e).__name__}: {e}"
        traceback.print_exc()
    sys.stdout.flush()
    if status:
        with open(status, "w") as f:
            f.write(state.split("\n")[0][:500] + "\n")
    if state != "ok" and "pymol" not in sys.modules:
        sys.exit(1)


main()
