#!/usr/bin/env python3
"""Move the one freshly built f2py extension into the package directory."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--module", required=True)
    parser.add_argument("--destination", required=True)
    args = parser.parse_args()
    destination = Path(args.destination).resolve()
    candidates = sorted(
        (
            path
            for path in Path.cwd().rglob(args.module)
            if destination not in path.resolve().parents
        ),
        key=lambda path: path.stat().st_mtime,
    )
    if not candidates:
        raise SystemExit(f"No extension matched {args.module}")
    destination.mkdir(parents=True, exist_ok=True)
    for stale in destination.glob("_core*.so"):
        stale.unlink()
    shutil.move(str(candidates[-1]), destination / candidates[-1].name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
