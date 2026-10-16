"""Command line: python -m ncaaprotein <stage> config.toml"""
import argparse
import sys

from .core import STAGES, StageError, load_config
from .stages import relax, ncaa, linked, free  # noqa: F401  (importing registers the stages)


def copy_mdp(dest):
    """Copy the package's mdp templates into a folder of your own, to edit and name as [project] mdp_dir."""
    import shutil
    from pathlib import Path
    from .core import MDP
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    made, kept = [], []
    for f in sorted(MDP.glob("*.mdp")):
        if (dest / f.name).exists():
            kept.append(f.name)
        else:
            shutil.copy(f, dest / f.name)
            made.append(f.name)
    print(f"copied {len(made)} mdp files to {dest}" + (f"; left {len(kept)} that were already there untouched" if kept else ""))
    print(f'edit them, then add  mdp_dir = "{dest}"  under [project] in your config. A file you delete from that folder falls back\n'
          "to the package's. The toolbox still sets these keys itself, from the config: nsteps (run lengths, worked out from\n"
          "your dt; equilibration lengths only when the config gives them), define (restraints), tc-grps and compressed-x-grps (the solute group), and ref_t / gen_temp when the config\n"
          "gives a temperature.")


def main():
    if len(sys.argv) >= 2 and sys.argv[1] == "copy-mdp":
        if len(sys.argv) != 3:
            sys.exit("usage: python -m ncaaprotein copy-mdp FOLDER")
        return copy_mdp(sys.argv[2])
    listing = "\n".join(f"  {k:<15} {('steps ' if '-' in v.steps else 'step ') + v.steps:<12} {(v.__doc__ or '').strip().splitlines()[0]}"
                        for k, v in STAGES.items())
    ap = argparse.ArgumentParser(prog="python -m ncaaprotein", formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description="Run one stage of the pipeline from a config file.",
                                 epilog="stages, in pipeline order:\n" + listing +
                                 "\n\nalso:  python -m ncaaprotein copy-mdp FOLDER   copies the mdp templates for you to edit ([project] mdp_dir)")
    ap.add_argument("stage", nargs="+", choices=list(STAGES), metavar="stage", help="one or more stages, run in the order given")
    ap.add_argument("config", help="the TOML config file")
    a = ap.parse_args()
    try:
        cfg = load_config(a.config)
    except StageError as e:
        print(f"\nSTOPPED: {e}", file=sys.stderr)
        sys.exit(1)
    for s in a.stage:
        try:
            STAGES[s](cfg)
        except StageError as e:
            print(f"\nSTOPPED in {s}: {e}", file=sys.stderr)
            sys.exit(1)
