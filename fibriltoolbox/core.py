"""Shared machinery: the config file, running commands, step bookkeeping, the script wrapper, mdp editing."""
import json
import os
import re
import subprocess
import sys
import time
import tomllib
from pathlib import Path

PKG = Path(__file__).resolve().parent
SCRIPTS = PKG / "scripts"
MDP = PKG / "mdp"
VERDICT = re.compile(r"\[(PASS|FAIL|WARN|info)\]|VERDICT|ALL CHECKS PASSED|problem")


class StageError(RuntimeError):
    """A step failed; the message says which command and where its log is."""


def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def load_config(path):
    path = Path(path).resolve()
    with open(path, "rb") as f:
        cfg = tomllib.load(f)
    cfg.setdefault("project", {})
    cfg.setdefault("tools", {})
    cfg["_dir"] = path.parent
    return cfg


class Ctx:
    """One stage's working folder, log, list of finished steps and report."""

    def __init__(self, cfg, stage, folder):
        self.cfg = cfg
        t = cfg["tools"]
        self.gmx = t.get("gmx", "gmx_mpi")
        self.pymol = t.get("pymol", "pymol")
        self.python = t.get("python", sys.executable)
        self.mdrun_flags = t.get("mdrun_flags", "-nb gpu -pme gpu")
        root = Path(cfg["project"].get("workdir", "."))
        self.root = root if root.is_absolute() else (cfg["_dir"] / root).resolve()
        self.stage = stage
        self.dir = self.root / folder
        (self.dir / "logs").mkdir(parents=True, exist_ok=True)
        self._state = self.dir / "steps_done.json"
        self.done = json.loads(self._state.read_text()) if self._state.exists() else {}
        self._cur = "setup"
        self._n = 0
        self._report = []
        self.say(f"##### stage {stage} in {self.dir}")

    # ---- paths
    def input(self, p):
        """A path from the config: absolute, or relative to the config file."""
        p = Path(os.path.expanduser(str(p)))
        return p if p.is_absolute() else (self.cfg["_dir"] / p).resolve()

    def other(self, folder):
        return self.root / folder

    # ---- talking
    def say(self, msg):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        print(line, flush=True)
        with open(self.dir / "logs" / "stage.log", "a") as f:
            f.write(line + "\n")

    def note(self, text):
        """A result line: printed, and kept for REPORT.md."""
        self.say("   " + text)
        self._report.append("- " + text)

    def grab(self, out, pattern=VERDICT):
        lines = [l.strip() for l in out.split("\n") if pattern.search(l)]
        for l in lines:
            self.note(l)
        return lines

    # ---- steps
    def begin(self, key, title):
        if key in self.done:
            self.say(f"-- already done ({self.done[key]}): {title}")
            return False
        self._cur = key
        self.say(f"== {title}")
        self._report.append(f"\n### {title}\n")
        return True

    def end(self):
        self.done[self._cur] = now()
        self._state.write_text(json.dumps(self.done, indent=1))

    # ---- gates: a step's checks decide whether the stage goes on
    def _gates(self):
        f = self.dir / "gates.json"
        return json.loads(f.read_text()) if f.exists() else {}

    def gate(self, key, lines):
        """Record the verdict lines of a finished step. Anything marked FAIL, or a minimisation that did not converge,
        holds the stage at check_gate()."""
        fails = [l for l in lines if "[FAIL]" in l or "NOT converged" in l]
        g = self._gates()
        g[key] = fails
        (self.dir / "gates.json").write_text(json.dumps(g, indent=1))

    def check_gate(self, key, section, what, look_at, nxt):
        """Go on if the step passed, or if the config accepts its failure; otherwise stop before `nxt`.
        section: the config table that takes accept = [...], as text ("[linked]") and as the table itself."""
        where, table = section
        fails = self._gates().get(key) or []
        if not fails:
            return
        if key in table.get("accept", []):
            self.note(f"{what}: {len(fails)} check(s) did not pass; accepted in the config ({where} accept), going on to {nxt}")
            return
        raise StageError(
            f"{what} did not meet the protocol's criteria, so {nxt} was not started:\n" + "\n".join("    " + l for l in fails) +
            f"\n  Look at {look_at} and REPORT.md under {self.dir}." +
            f"\n  If this is acceptable for your system: add  accept = [\"{key}\"]  under {where} in the config and run the "
            f"same command again; it carries on from here."
            f"\n  If it is not: fix the cause, delete the entry of that step from {self._state.name} in that folder (the "
            f"finished steps, in order), and run again to repeat it.")

    def pause(self, section, nxt):
        where, table = section
        if table.get("stop_before") == nxt:
            raise StageError(f"paused before {nxt}, as {where} stop_before asks. Look at the checks and plots so far; to go on, "
                             f"remove that line from the config and run the same command again.")

    def finish(self, ok=True):
        if self._report:
            with open(self.dir / "REPORT.md", "a") as f:
                f.write(f"\n## {self.stage}, {now()}" + ("" if ok else " (stopped before the end)") + "\n" + "\n".join(self._report) + "\n")
        self.say(f"##### stage {self.stage} {'finished' if ok else 'STOPPED before the end'}; report: {self.dir / 'REPORT.md'}")

    # ---- running things
    def sh(self, cmd, cwd=None, stdin=None, log=None, check=True, env=None):
        """Run a shell command, keep everything it prints in logs/, return that text."""
        cwd = Path(cwd or self.dir)
        self._n += 1
        logf = self.dir / "logs" / f"{self._cur}_{self._n:03d}_{log or 'cmd'}.log"
        with open(self.dir / "commands.log", "a") as f:
            f.write(f"[{now()}] (cd {cwd})\n")
            for k, v in (env or {}).items():
                f.write(f"    export {k}='{v}'\n")
            f.write(f"    {cmd}" + (f"      <<< {stdin!r}" if stdin else "") + "\n")
        r = subprocess.run(["bash", "-c", cmd], cwd=cwd, input=stdin, text=True, errors="replace",
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env={**os.environ, **(env or {})})
        logf.write_text(r.stdout)
        if check and r.returncode != 0:
            tail = "\n".join(r.stdout.strip().split("\n")[-25:])
            raise StageError(f"command failed (exit {r.returncode}): {cmd}\n   full log: {logf}\n{tail}")
        return r.stdout

    def script(self, name, args=(), overrides=None, cwd=None, pymol=False, log=None, check=True):
        """Run one of the protocol's scripts. `overrides` replaces settings written at the top of the script
        (NAME = value) without editing the file; `pymol=True` runs it in headless PyMOL."""
        cwd = Path(cwd or self.dir)
        src = name if os.path.isabs(str(name)) else SCRIPTS / name
        self._n += 1
        status = self.dir / "logs" / f"{self._cur}_{self._n:03d}_{Path(name).stem}.status"
        env = {"FTB_SCRIPT": str(src), "FTB_ARGS": json.dumps([str(a) for a in args]),
               "FTB_OVERRIDES": json.dumps(overrides or {}), "FTB_STATUS": str(status)}
        runner = f"{self.pymol} -cq" if pymol else self.python
        out = self.sh(f"{runner} {PKG / 'runscript.py'}", cwd=cwd, log=log or Path(name).stem, check=False, env=env)
        state = status.read_text().strip() if status.exists() else "no status written (the interpreter did not start?)"
        if state != "ok" and check:
            tail = "\n".join(out.strip().split("\n")[-25:])
            raise StageError(f"{name} stopped: {state}\n{tail}")
        return out


def write_mdp(template, dest, changes=None):
    """Copy an mdp file from the package, changing only the keys named in `changes` (None removes a key)."""
    def norm(k):
        return k.strip().replace("-", "_").lower()
    todo = {norm(k): (k, v) for k, v in (changes or {}).items()}
    out = []
    for line in (MDP / template).read_text().rstrip("\n").split("\n"):
        body = line.split(";", 1)[0]
        if "=" in body and norm(body.split("=", 1)[0]) in todo:
            key, val = todo.pop(norm(body.split("=", 1)[0]))
            if val is None:
                out.append("; " + line)
            else:
                out.append(f"{body.split('=', 1)[0]}= {val}")
            continue
        out.append(line)
    out += [f"{k:<24}= {v}" for k, v in todo.values() if v is not None]
    Path(dest).write_text("\n".join(out) + "\n")
    return dest


def mdp_value(template, key):
    for line in (MDP / template).read_text().split("\n"):
        body = line.split(";", 1)[0]
        if "=" in body and body.split("=", 1)[0].strip().replace("-", "_") == key.replace("-", "_"):
            return body.split("=", 1)[1].strip()
    return None


def nsteps_for(template, ns):
    """Number of steps that makes `ns` nanoseconds at the time step of an mdp file."""
    return int(round(float(ns) * 1000.0 / float(mdp_value(template, "dt"))))


# ---------------------------------------------------------------- PDB helpers
def pdb_atoms(path):
    return [l for l in open(path) if l.startswith(("ATOM", "HETATM"))]


def pdb_chains(path):
    """Chains of a PDB in file order: [(chain id, [resid, ...]), ...]. A new chain starts where the chain ID
    changes, at a TER, or where the residue number drops."""
    chains, prev, ter = [], None, True
    for l in open(path):
        if l.startswith("TER"):
            ter = True
            continue
        if not l.startswith(("ATOM", "HETATM")):
            continue
        ch, resi = l[21], int(l[22:26])
        if ter or ch != prev[0] or resi < prev[1]:
            chains.append((ch, []))
        if not chains[-1][1] or chains[-1][1][-1] != resi:
            chains[-1][1].append(resi)
        prev, ter = (ch, resi), False
    return chains


def resindex_ranges(chains, spec):
    """'G' or 'L:8-73' -> (first, last) residue index counted straight through the structure from 1."""
    ch, _, rng = spec.partition(":")
    start = 1
    for cid, resids in chains:
        if cid == ch:
            if not rng:
                return start, start + len(resids) - 1
            lo, hi = (int(x) for x in rng.split("-"))
            idx = [i for i, r in enumerate(resids) if lo <= r <= hi]
            if not idx:
                raise StageError(f"chain {ch} has no residues {rng}")
            return start + idx[0], start + idx[-1]
        start += len(resids)
    raise StageError(f"no chain {ch!r} in the structure (chains: {''.join(c for c, _ in chains)})")


def chain_breaks(path, limit=4.2):
    """Consecutive C-alpha atoms of a chain further apart than `limit` A: [(chain, resid, distance), ...]."""
    out, prev = [], None
    for l in pdb_atoms(path):
        if l[12:16] != " CA ":
            continue
        cur = (l[21], float(l[30:38]), float(l[38:46]), float(l[46:54]))
        if prev and prev[0] == cur[0]:
            d = sum((a - b) ** 2 for a, b in zip(cur[1:], prev[1:])) ** 0.5
            if d > limit:
                out.append((cur[0], int(l[22:26]), d))
        prev = cur
    return out


# ---------------------------------------------------------------- stage registry
STAGES = {}


def stage(name, folder, steps):
    """Register a stage: its command name, its folder under the work directory, the protocol steps it covers."""
    def deco(fn):
        def run(cfg):
            ctx = Ctx(cfg, name, cfg.get("folders", {}).get(name, folder))
            try:
                fn(ctx)
            except BaseException:
                ctx.finish(ok=False)
                raise
            ctx.finish()
        run.steps, run.__doc__ = steps, fn.__doc__
        STAGES[name] = run
        return run
    return deco
