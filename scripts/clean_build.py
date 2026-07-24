#!/usr/bin/env python3
"""Remove only generated build products from this repository."""

from __future__ import annotations

import shutil
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    for directory in (root / "build", root / "dist"):
        if directory.exists():
            shutil.rmtree(directory)
    for pattern in ("*.o", "*.mod", "_core*.so"):
        for path in root.glob(pattern):
            path.unlink()
        for path in (root / "src" / "splas").glob(pattern):
            path.unlink()
    for cache in root.rglob("__pycache__"):
        shutil.rmtree(cache)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
