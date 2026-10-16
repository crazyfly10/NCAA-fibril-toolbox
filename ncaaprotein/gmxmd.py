"""GROMACS helpers: grompp, a detached mdrun that can be resumed, energies, the verdict scripts,
and the periodic-image check."""
import os
import re
import shlex
import subprocess
import time
from pathlib import Path

from .core import StageError, now


def grompp(ctx, d, mdp, c, p, o, r=None, t=None, n=None, maxwarn=0):
    cmd = f"{ctx.gmx} grompp -f {mdp} -c {c} -p {p} -o {o}"
    cmd += (f" -r {r}" if r else "") + (f" -t {t}" if t else "") + (f" -n {n}" if n else "")
    cmd += f" -maxwarn {maxwarn}" if maxwarn else ""
    return ctx.sh(cmd, cwd=d, log="grompp_" + Path(o).stem)


def _tail(path, n=4000):
    with open(path, "rb") as f:
        f.seek(0, 2)
        f.seek(max(0, f.tell() - n))
        return f.read().decode(errors="replace")


def finished(d, deffnm):
    log = Path(d) / f"{deffnm}.log"
    return log.exists() and "Finished mdrun" in _tail(log) and (Path(d) / f"{deffnm}.gro").exists()


def _alive(pid, d):
    try:
        os.kill(pid, 0)
        return os.path.realpath(f"/proc/{pid}/cwd") == os.path.realpath(d)     # ours, not a recycled PID
    except (OSError, ValueError):
        return False


def mdrun(ctx, d, deffnm, flags=None):
    """Run `gmx mdrun -deffnm X` in its own session, with its process ID in X.pid and its output in X.out.
    The run survives the toolbox (or the terminal) being closed. Called again, this waits for a run that is
    still going, or continues one that stopped from its checkpoint. It never stops a run."""
    d = Path(d)
    pidf = d / f"{deffnm}.pid"
    if finished(d, deffnm):
        return
    pid = int(pidf.read_text().split()[0]) if pidf.exists() and pidf.read_text().strip() else None
    if pid and _alive(pid, d):
        ctx.say(f"   {deffnm}: already running (process {pid}); waiting for it")
        while _alive(pid, d):
            time.sleep(20)
    else:
        cmd = [ctx.gmx, "mdrun", "-deffnm", deffnm] + shlex.split(ctx.mdrun_flags if flags is None else flags)
        out = d / f"{deffnm}.out"
        if (d / f"{deffnm}.cpt").exists():
            cmd += ["-cpi", f"{deffnm}.cpt"]
            k = 1
            while (d / f"{deffnm}_cont{k}.out").exists():
                k += 1
            out = d / f"{deffnm}_cont{k}.out"
            ctx.say(f"   {deffnm}: continuing from {deffnm}.cpt")
        with open(ctx.dir / "commands.log", "a") as f:
            pre = "".join(f"{k}={v} " for k, v in ctx.mdrun_env.items())
            f.write(f"[{now()}] (cd {d})\n    {pre}nohup {' '.join(cmd)} > {out.name} 2>&1 &   # PID in {pidf.name}\n")
        with open(out, "w") as fo:
            p = subprocess.Popen(cmd, cwd=d, stdout=fo, stderr=subprocess.STDOUT, start_new_session=True,
                                 env={**os.environ, **ctx.mdrun_env})
            pidf.write_text(f"{p.pid}\n")
            p.wait()
    if not finished(d, deffnm):
        outs = sorted(d.glob(f"{deffnm}*.out"), key=lambda q: q.stat().st_mtime)
        raise StageError(f"mdrun {deffnm} did not finish in {d}\n" + (_tail(outs[-1], 1500) if outs else ""))


def lincs(d, deffnm):
    return (Path(d) / f"{deffnm}.log").read_text(errors="replace").count("LINCS WARNING")


def performance(d, deffnm):
    for l in _tail(Path(d) / f"{deffnm}.log", 3000).split("\n"):
        if l.startswith("Performance:"):
            return float(l.split()[1])
    return None


def energy(ctx, d, edr, terms, xvg):
    return ctx.sh(f"{ctx.gmx} energy -f {edr} -o {xvg}", cwd=d, stdin="\n".join(terms) + "\n0\n", log="energy_" + Path(xvg).stem)


def draw(ctx, d, script, args=()):
    """Run a plotting script on the files an analysis left. A figure that fails is reported, and never stops a stage:
    the numbers and verdicts are already written."""
    out = ctx.script(script, list(args), cwd=d, check=False)
    made = [l.strip() for l in out.split("\n") if l.startswith("wrote ")]
    if not made:
        ctx.note(f"{script} drew no figure: " + (out.strip().split("\n")[-1][:200] if out.strip() else "no output"))
    return made


def run_em(ctx, d, mdp, c, p, n=None, maxwarn=0, second_pass=True):
    """Minimise; returns the name of the minimised structure. If the first pass ends NOT converged (it stopped at
    machine precision with Fmax above 1000), one more pass from em.gro starts the step size afresh."""
    grompp(ctx, d, mdp, c, p, "em.tpr", n=n, maxwarn=maxwarn)
    mdrun(ctx, d, "em", flags="")
    energy(ctx, d, "em.edr", ["Potential"], "potential.xvg")
    v = ctx.grab(ctx.script("check_em.py", ["potential.xvg", "em.log"], cwd=d))
    draw(ctx, d, "plot_em.py", ["em.csv", "em.png"])
    if second_pass and any("NOT converged" in l for l in v):
        ctx.note("first pass NOT converged: one more pass from em.gro")
        grompp(ctx, d, mdp, "em.gro", p, "em2.tpr", n=n, maxwarn=maxwarn)
        mdrun(ctx, d, "em2", flags="")
        energy(ctx, d, "em2.edr", ["Potential"], "potential2.xvg")
        v = ctx.grab(ctx.script("check_em.py", ["potential2.xvg", "em2.log"], cwd=d))
        draw(ctx, d, "plot_em.py", ["em.csv", "em.png"])
        ctx.gate("em", [l for l in v if l.startswith("VERDICT")])
        return "em2.gro"
    ctx.gate("em", [l for l in v if l.startswith("VERDICT")])
    return "em.gro"


def mdp_step(path):
    """dt (ps) and mass-repartition-factor set in an mdp file"""
    found = {}
    for line in open(path):
        body = line.split(";", 1)[0]
        if "=" in body:
            k, v = (x.strip() for x in body.split("=", 1))
            if k.replace("_", "-") in ("dt", "mass-repartition-factor") and v:
                found[k.replace("_", "-")] = float(v.split()[0])
    return found.get("dt", 0.001), found.get("mass-repartition-factor", 1.0)


def run_md(ctx, d, name, mdp, c, p, r=None, t=None, n=None, maxwarn=0, flags=None):
    if name in ("nvt", "npt"):
        dt, mrf = mdp_step(Path(d) / mdp)
        if dt > 0.0025 and mrf > 1:
            f = LARGE_STEP_MDP.get((ctx.stage, name), mdp)
            ctx.note(f"note: this {name.upper()} runs at {dt * 1000:g} fs with repartitioned hydrogen masses ({f}). The solute "
                     f"group then often reads 3-5 K above ref_t while the water stays at it, and fails the temperature check. "
                     f"To avoid that, run it at 2 fs: dt = 0.002 and mass-repartition-factor = 1 in your copy of {f} "
                     f"([project] mdp_dir; README, \"Your own mdp files\")")
    grompp(ctx, d, mdp, c, p, f"{name}.tpr", r=r, t=t, n=n, maxwarn=maxwarn)
    mdrun(ctx, d, name, flags=flags)
    ns = performance(d, name)
    ctx.note(f"{name}: finished, {lincs(d, name)} LINCS warnings" + (f", {ns:.0f} ns/day" if ns else ""))


def time_step(d, kind):
    """dt (ps) and mass-repartition-factor of a finished run, from its log"""
    found = {}
    for l in open(Path(d) / f"{kind}.log"):
        m = re.match(r"\s+(dt|mass-repartition-factor)\s+=\s+(\S+)", l)
        if m:
            found.setdefault(m.group(1), float(m.group(2)))
        if len(found) == 2:
            break
    return found.get("dt", 0.002), found.get("mass-repartition-factor", 1.0)


# the packaged mdp file behind each equilibration (free-run minimises and runs NVT with the linked files)
LARGE_STEP_MDP = {("protein-run", "nvt"): "protein_nvt.mdp", ("protein-run", "npt"): "protein_npt.mdp",
                  ("linked-run", "nvt"): "linked_nvt.mdp", ("linked-run", "npt"): "linked_npt.mdp",
                  ("free-run", "nvt"): "linked_nvt.mdp", ("free-run", "npt"): "free_npt.mdp"}


def _judge_once(ctx, d, kind, groups, ref_t, density):
    tg = [f"T-{g}" for g in groups if g != "System"]
    s = {"REF_T": ref_t}
    dt, mrf = time_step(d, kind)
    if dt > 0.0025 and mrf > 1:
        # at 4 fs with repartitioned hydrogen masses the solute can read a few K warm while the water stays at ref_t:
        # the check then names the mdp file to run at 2 fs instead
        s["LARGE_STEP_MDP"] = LARGE_STEP_MDP.get((ctx.stage, kind), f"the {kind} mdp file")
    if kind == "nvt":
        energy(ctx, d, "nvt.edr", ["Temperature"] + tg + ["Potential", "Pressure"], "nvt.xvg")
        v = ctx.grab(ctx.script("check_nvt.py", ["nvt.xvg"], s, cwd=d))
        draw(ctx, d, "plot_nvt.py")
    else:
        energy(ctx, d, "npt.edr", ["Temperature"] + tg + ["Pressure", "Density", "Potential"], "npt.xvg")
        v = ctx.grab(ctx.script("check_npt.py", ["npt.xvg"] + [str(x) for x in (density or [])], s, cwd=d))
        draw(ctx, d, "plot_npt.py")
    ctx.gate(kind, v)
    return v


def _run_ps(d, kind):
    """Length of a finished run in ps: the last time in its energy trace."""
    last = [l for l in open(Path(d) / f"{kind}.xvg") if l.strip() and l[0] not in "#@"][-1]
    return float(last.split()[0])


def judge(ctx, d, kind, groups, ref_t=298.0, density=None, table=None, flags=None):
    """The verdicts of check_nvt.py / check_npt.py for a finished stage; plot_nvt.py / plot_npt.py draw the figure.
    `groups` are the thermostat groups.

    If a line fails, the run is continued for its own length again and checked afresh, up to `extend_max` times
    (2 unless the stage's config table says otherwise; 0 switches this off). The checks look at the last half of
    the run, so each time they judge new data: a failure that was chance passes, a system that is not
    equilibrated keeps failing and holds the stage."""
    table = table or {}
    v = _judge_once(ctx, d, kind, groups, ref_t, density)
    n, k = int(table.get("extend_max", 2)), 0
    failed = lambda lines: [l for l in lines if "[FAIL]" in l]
    if not failed(v) or kind in table.get("accept", []):
        return v
    first = _run_ps(d, kind)
    while failed(v) and k < n:
        k += 1
        ctx.note(f"{kind.upper()}: {len(failed(v))} check(s) did not pass after {_run_ps(d, kind):g} ps; continuing for another "
                 f"{first:g} ps and checking again (extension {k} of at most {n})")
        gro = Path(d) / f"{kind}.gro"
        gro.rename(Path(d) / f"{kind}_before_extension{k}.gro")
        ctx.sh(f"{ctx.gmx} convert-tpr -s {kind}.tpr -extend {first:g} -o {kind}.tpr", cwd=d, log=f"extend_{kind}")
        mdrun(ctx, d, kind, flags=flags)
        ctx.note(f"{kind}: now {_run_ps_log(d, kind)}, {lincs(d, kind)} LINCS warnings")
        v = _judge_once(ctx, d, kind, groups, ref_t, density)
    if k:
        ctx.note(f"{kind.upper()}: " + (f"every check passes after {k} extension(s), {_run_ps(d, kind):g} ps in all" if not failed(v)
                                          else f"still {len(failed(v))} check(s) failing after {k} extension(s), {_run_ps(d, kind):g} ps in all"))
    return v


def _run_ps_log(d, kind):
    ns = performance(d, kind)
    return "extended" + (f" at {ns:.0f} ns/day" if ns else "")


def box_of(gro):
    return _tail(gro, 300).strip().split("\n")[-1].split()


# what a box of solute, OPC water and salt weighs per volume: water with 0.15 M NaCl as measured in this pipeline
# (ncaa-relax, 1003 kg/m3), and a typical partial specific volume for protein (0.73 cm3/g)
RHO_WATER, RHO_SOLUTE = 1003.0, 1370.0
_MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06, "P": 30.974, "F": 18.998}
_ION = {"NA": 22.990, "CL": 35.453, "K": 39.098, "MG": 24.305, "CA": 40.078, "ZN": 65.38}
_WATER = {"SOL", "WAT", "HOH"}


def expected_density(gro, window=10.0):
    """The density range to expect after NPT, worked out from what is in the box: more water around the same solute
    means a lower density, so no fixed range fits every box. Used when the config gives no density_range."""
    lines = open(gro).read().split("\n")
    m_wat = m_sol = 0.0
    for l in lines[2:2 + int(lines[1])]:
        res, name = l[5:10].strip(), l[10:15].strip()
        if res in _WATER:
            m_wat += _MASS.get(name[0], 0.0)              # the fourth site of OPC has no mass
        elif res.upper().rstrip("+-") in _ION:
            m_wat += _ION[res.upper().rstrip("+-")]
        else:
            m_sol += _MASS.get(name.lstrip("0123456789")[:1], 12.011)
    rho = (m_wat + m_sol) / (m_wat / RHO_WATER + m_sol / RHO_SOLUTE)
    return [round(rho - window), round(rho + window)]


def density_range(ctx, given, gro):
    """The range from the config, or the one worked out from the box, with a line in the report saying which."""
    if given:
        return list(given)
    r = expected_density(gro)
    ctx.note(f"expected density {r[0]}-{r[1]} kg/m3, worked out from the water, salt and solute in the box "
             f"(no density_range in the config)")
    return r


def image_verdict(xvg):
    """Smallest periodic-image distance in a `gmx mindist -pi` file, judged against 2.0 nm (PASS) and 1.2 nm (FAIL)."""
    rows = [(float(l.split()[0]), float(l.split()[1])) for l in open(xvg) if l.strip() and l[0] not in "#@"]
    t, dmin = min(rows, key=lambda x: x[1])
    b = sum(1 for _, x in rows if x < 2.0)
    c = sum(1 for _, x in rows if x < 1.2)
    text = f"smallest periodic image distance {dmin:.3f} nm at {t:.0f} ps; under 2.0 nm in {b} of {len(rows)} frames, under 1.2 nm in {c}"
    verdict = ("FAIL: within reach of its own image (cutoff 1.0 nm) -- rebuild with a larger box" if c else
               "OK, but it dips under 2.0 nm -- see which part reaches out before going on" if b else
               "PASS: never under 2.0 nm")
    return dmin, text, verdict


def image_clean_stretch(xvg, floor=1.2):
    """After a FAIL: the stretch at the end of the run in which the distance stays above the floor: the part
    of the run that can still be analysed. Returns None when the distance never goes under the floor."""
    rows = [(float(l.split()[0]), float(l.split()[1])) for l in open(xvg) if l.strip() and l[0] not in "#@"]
    bad = [t for t, x in rows if x < floor]
    if not bad:
        return None
    tail = [(t, x) for t, x in rows if t > bad[-1]]
    if not tail:
        return f"the distance is under {floor:g} nm in the last frame measured: no clean stretch at the end of the run"
    return (f"last frame under {floor:g} nm at {bad[-1]:.0f} ps; from {tail[0][0]:.0f} ps to the end "
            f"({(tail[-1][0] - tail[0][0]) / 1000:.1f} ns) it never goes under {min(x for _, x in tail):.2f} nm")


def mindist_pi(ctx, d, tpr, xtc, group, dt=None, out="mindist_pi.xvg"):
    """gmx mindist -pi: the smallest distance between the solute and its own periodic images, frame by frame.
    It runs on one CPU core. With [tools] cpu_jobs above 1 the trajectory is cut into that many stretches of time,
    measured side by side and joined: each frame is measured on its own, so the result is the same."""
    step = f" -dt {int(dt)}" if dt else ""
    jobs = ctx.cpu_jobs
    if jobs > 1:
        last = re.findall(r"Last frame\s+\d+\s+time\s+([0-9.]+)", ctx.sh(f"{ctx.gmx} check -f {xtc}", cwd=d, log="check_traj"))
        end = float(last[-1]) if last else 0.0
        grain = float(dt or 1)
        cuts = [round(end * i / jobs / grain) * grain for i in range(jobs + 1)]
        cuts[-1] = end
        if end <= 0 or len(set(cuts)) < 3:
            jobs = 1
    if jobs == 1:
        ctx.sh(f"{ctx.gmx} mindist -s {tpr} -f {xtc} -pi -od {out}{step}", cwd=d, stdin=f"{group}\n", log="mindist")
        return
    parts = [f"{Path(out).stem}_part{i}.xvg" for i in range(jobs)]
    runs = " ".join(f"(echo {shlex.quote(group)} | {ctx.gmx} mindist -s {tpr} -f {xtc} -pi -od {parts[i]}{step} "
                    f"-b {cuts[i]:g} -e {cuts[i + 1]:g} > {Path(parts[i]).stem}.log 2>&1) &" for i in range(jobs))
    ctx.sh(f"{runs} wait", cwd=d, log="mindist_parallel")
    head, rows, seen = [], [], -1.0
    for i, p in enumerate(parts):
        for l in open(Path(d) / p):
            if l[:1] in "#@":
                if i == 0:
                    head.append(l)
            elif l.strip() and float(l.split()[0]) > seen:
                rows.append(l)
                seen = float(l.split()[0])
    head.insert(0, f"# joined from {jobs} gmx mindist runs over consecutive stretches of time ([tools] cpu_jobs); "
                   f"the commands are in commands.log\n")
    (Path(d) / out).write_text("".join(head + rows))
    ctx.note(f"image distance measured in {jobs} parallel parts ([tools] cpu_jobs = {ctx.cpu_jobs}): {len(rows)} frames")


def image_check_single(ctx, d, tpr, xtc, ndx, group, dt=100):
    """Solute that is ONE molecule (tleap/acpype systems): make it whole, then measure."""
    ctx.sh(f"{ctx.gmx} convert-tpr -s {tpr} -n {ndx} -o solute_only.tpr", cwd=d, stdin=f"{group}\n", log="convert_tpr")
    ctx.sh(f"{ctx.gmx} trjconv -s {tpr} -f {xtc} -n {ndx} -o whole_solute.xtc -pbc whole -dt {dt}", cwd=d, stdin=f"{group}\n", log="trjconv_whole")
    mindist_pi(ctx, d, "solute_only.tpr", "whole_solute.xtc", "0")
    dmin, text, verdict = image_verdict(Path(d) / "mindist_pi.xvg")
    ctx.note(text)
    ctx.note(verdict)
    clean = image_clean_stretch(Path(d) / "mindist_pi.xvg")
    if clean:
        ctx.note(clean)
    return dmin
