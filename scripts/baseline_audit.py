#!/usr/bin/env python3
"""Machine-readable checks that supplement docs/baseline_audit.md."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    historical = [
        root / "functions.f90",
        root / "solver.f90",
        root / "solver_closed_system.f90",
        root / "splas.f90",
        root / "splas_py_wrapper.f90",
        root / "legacy" / "splas.py",
        root / "legacy" / "main.py",
        root / "input_default.json",
    ]
    payload = {
        "files": {
            str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in historical
        },
        "environment": {
            "gfortran": shutil.which("gfortran"),
            "python": shutil.which("python3"),
        },
        "known_missing_original_directories": [
            "input",
            "inputIC",
            "input_closed_world_tv",
            "input_test_1.5",
            "sensitivity_json",
            "SPLAS_PICKLE",
        ],
    }
    destination = root / "docs" / "baseline_audit.json"
    destination.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(destination)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
