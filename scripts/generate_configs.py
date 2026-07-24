#!/usr/bin/env python3
"""Generate the versioned thesis experiment configs and manifests."""

from __future__ import annotations

import json
from pathlib import Path

from splas.config import config_hash, load_config


LEVELS = [
    ("I_m30", 0.70),
    ("I_m20", 0.80),
    ("I_m10", 0.90),
    ("I_default", 1.00),
    ("I_p10", 1.10),
    ("I_p20", 1.20),
    ("I_p30", 1.30),
]


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def make_config(path: Path, extends: str, name: str, multiplier: float, overrides: dict | None = None) -> None:
    payload: dict = {
        "extends": extends,
        "name": name,
        "forcing": {"multiplier": multiplier},
    }
    if overrides:
        for key, value in overrides.items():
            payload[key] = value
    write_json(path, payload)


def make_manifest(root: Path, name: str, config_dir: str, entries: list[tuple[str, Path]]) -> None:
    runs = []
    for run_name, config in entries:
        resolved = load_config(config)
        runs.append(
            {
                "name": run_name,
                "config": str(Path("..") / config_dir / config.name),
                "duration_days": resolved["numerics"]["duration_days"],
                "timestep_days": resolved["numerics"]["timestep_days"],
                "output_interval_days": resolved["numerics"]["output_interval_days"],
                "seed": resolved["numerics"]["seed"],
                "scientific_profile": resolved["scientific_profile"],
                "config_hash": config_hash(resolved),
            }
        )
    write_json(root / "configs" / "manifests" / f"{name}.json", {"name": name, "runs": runs})


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    groups: dict[str, list[tuple[str, Path]]] = {}
    for label, multiplier in LEVELS:
        constant_open = root / "configs" / "constant_light" / f"{label}_open.json"
        make_config(
            constant_open,
            "../base/thesis_reported.json",
            f"constant_open_{label}",
            multiplier,
        )
        groups.setdefault("constant_open", []).append((f"constant_open_{label}", constant_open))

        constant_closed = root / "configs" / "constant_light" / f"{label}_closed.json"
        make_config(
            constant_closed,
            "../reference/reference_closed_constant_default.json",
            f"constant_closed_{label}",
            multiplier,
        )
        groups.setdefault("constant_closed", []).append((f"constant_closed_{label}", constant_closed))

        seasonal = root / "configs" / "seasonal_light" / f"{label}.json"
        make_config(
            seasonal,
            "../base/thesis_reported.json",
            f"seasonal_200m_{label}",
            multiplier,
            {"forcing": {"mode": "seasonal", "multiplier": multiplier}},
        )
        groups.setdefault("seasonal_200m", []).append((f"seasonal_200m_{label}", seasonal))

        critical = root / "configs" / "critical_depth_100m" / f"{label}.json"
        make_config(
            critical,
            "../base/thesis_reported.json",
            f"seasonal_100m_{label}",
            multiplier,
            {
                "forcing": {"mode": "seasonal", "multiplier": multiplier},
                "grid": {"layers": 21, "dz_m": 5.0},
            },
        )
        groups.setdefault("critical_depth_100m", []).append((f"seasonal_100m_{label}", critical))

    for name, entries in groups.items():
        make_manifest(root, name, entries[0][1].parent.name, entries)
    reference_entries = [
        ("reference_open_constant_default", root / "configs" / "reference" / "reference_open_constant_default.json"),
        ("reference_closed_constant_default", root / "configs" / "reference" / "reference_closed_constant_default.json"),
    ]
    make_manifest(root, "reference", "reference", reference_entries)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
