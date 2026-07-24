"""Command-line interface required by the SPLAS development contract."""

from __future__ import annotations

import argparse
from pathlib import Path

from .analysis import analyze_results
from .figures import make_run_figures
from .reproduce import reproduce_thesis
from .runner import run, run_ensemble


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(prog="splas")
    commands = root.add_subparsers(dest="command", required=True)

    single = commands.add_parser("run", help="run one validated SPLAS configuration")
    single.add_argument("config")
    single.add_argument("--output", required=True)

    ensemble = commands.add_parser("ensemble", help="run a versioned ensemble manifest")
    ensemble.add_argument("manifest")
    ensemble.add_argument("--jobs", type=int, default=1)
    ensemble.add_argument("--output", required=True)

    analyze = commands.add_parser("analyze", help="post-process a run dataset")
    analyze.add_argument("results_dir")
    analyze.add_argument("--output", required=True)

    figures = commands.add_parser("figures", help="make non-interactive figures for one run")
    figures.add_argument("results_dir")
    figures.add_argument("--output", required=True)

    reproduce = commands.add_parser("reproduce-thesis", help="run or resume the complete thesis suite")
    reproduce.add_argument("--jobs", type=int, default=1)
    reproduce.add_argument("--output", default="reproduction")
    return root


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.command == "run":
        run(args.config, args.output)
    elif args.command == "ensemble":
        run_ensemble(args.manifest, args.jobs, args.output)
    elif args.command == "analyze":
        analyze_results(args.results_dir, args.output)
    elif args.command == "figures":
        make_run_figures(args.results_dir, args.output)
    elif args.command == "reproduce-thesis":
        reproduce_thesis(args.output, args.jobs)
    return 0
