#!/usr/bin/env python3
from __future__ import annotations

import argparse

from splas.reproduce import reproduce_thesis


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="reproduction")
    parser.add_argument("--jobs", type=int, default=1)
    args = parser.parse_args()
    reproduce_thesis(args.output, args.jobs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
