"""GROMACS helpers: grompp, a detached mdrun that can be resumed, energies, the verdict scripts,
and the periodic-image check."""
import os
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
            f.write(f"[{now()}] (cd {d})\n    nohup {' '.join(cmd)} > {out.name} 2>&1 &   # PID in {pidf.name}\n")
        with open(out, "w") as fo:
            p = subprocess.Popen(cmd, cwd=d, stdout=fo, stderr=subprocess.STDOUT, start_new_session=True)
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


def run_md(ctx, d, name, mdp, c, p, r=None, t=None, n=None, maxwarn=0, flags=None):
    grompp(ctx, d, mdp, c, p, f"{name}.tpr", r=r, t=t, n=n, maxwarn=maxwarn)
    mdrun(ctx, d, name, flags=flags)
    ns = performance(d, name)
    ctx.note(f"{name}: finished, {lincs(d, name)} LINCS warnings" + (f", {ns:.0f} ns/day" if ns else ""))


def _judge_once(ctx, d, kind, groups, ref_t, density):
    tg = [f"T-{g}" for g in groups if g != "System"]
    if kind == "nvt":
        energy(ctx, d, "nvt.edr", ["Temperature"] + tg + ["Potential", "Pressure"], "nvt.xvg")
        v = ctx.grab(ctx.script("check_nvt.py", ["nvt.xvg"], {"REF_T": ref_t}, cwd=d))
        draw(ctx, d, "plot_nvt.py")
    else:
        energy(ctx, d, "npt.edr", ["Temperature"] + tg + ["Pressure", "Density", "Potential"], "npt.xvg")
        v = ctx.grab(ctx.script("check_npt.py", ["npt.xvg"] + [str(x) for x in (density or [])], {"REF_T": ref_t}, cwd=d))
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


def image_check_single(ctx, d, tpr, xtc, ndx, group, dt=100):
    """Solute that is ONE molecule (tleap/acpype systems): make it whole, then measure."""
    ctx.sh(f"{ctx.gmx} convert-tpr -s {tpr} -n {ndx} -o solute_only.tpr", cwd=d, stdin=f"{group}\n", log="convert_tpr")
    ctx.sh(f"{ctx.gmx} trjconv -s {tpr} -f {xtc} -n {ndx} -o whole_solute.xtc -pbc whole -dt {dt}", cwd=d, stdin=f"{group}\n", log="trjconv_whole")
    ctx.sh(f"{ctx.gmx} mindist -s solute_only.tpr -f whole_solute.xtc -pi -od mindist_pi.xvg", cwd=d, stdin="0\n", log="mindist")
    dmin, text, verdict = image_verdict(Path(d) / "mindist_pi.xvg")
    ctx.note(text)
    ctx.note(verdict)
    clean = image_clean_stretch(Path(d) / "mindist_pi.xvg")
    if clean:
        ctx.note(clean)
    return dmin
