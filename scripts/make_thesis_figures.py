#!/usr/bin/env python3
from __future__ import annotations

import argparse

from splas.figures import make_thesis_figures, write_figure_manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("reproduction_dir", nargs="?", default="reproduction")
    args = parser.parse_args()
    created = make_thesis_figures(args.reproduction_dir)
    write_figure_manifest(f"{args.reproduction_dir}/figure_manifest.csv", created)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
