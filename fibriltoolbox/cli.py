"""Command line: python -m fibriltoolbox <stage> config.toml"""
import argparse
import sys

from .core import STAGES, StageError, load_config
from .stages import ncaa, relax, linked, free  # noqa: F401  (importing registers the stages)


def main():
    listing = "\n".join(f"  {k:<15} protocol steps {v.steps:<10} {(v.__doc__ or '').strip().splitlines()[0]}"
                        for k, v in STAGES.items())
    ap = argparse.ArgumentParser(prog="python -m fibriltoolbox", formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description="Run one stage of the pipeline from a config file.",
                                 epilog="stages, in pipeline order:\n" + listing)
    ap.add_argument("stage", nargs="+", choices=list(STAGES), metavar="stage", help="one or more stages, run in the order given")
    ap.add_argument("config", help="the TOML config file")
    a = ap.parse_args()
    cfg = load_config(a.config)
    for s in a.stage:
        try:
            STAGES[s](cfg)
        except StageError as e:
            print(f"\nSTOPPED in {s}: {e}", file=sys.stderr)
            sys.exit(1)
