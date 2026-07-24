"""Single-run and ensemble orchestration with deterministic metadata."""

from __future__ import annotations

import concurrent.futures
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

from . import __version__
from .analysis import carbon_from_phyto, shannon_metrics
from .backend import core_source_hash, run_core
from .config import canonical_config, config_hash, load_config, write_resolved_config
from .dataset import write_dataset
from .traits import log_size_grid, phyto_traits


def _command_output(command: list[str]) -> str | None:
    try:
        return subprocess.run(command, check=True, text=True, capture_output=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _repository_root(config: dict[str, Any]) -> Path:
    return Path(config["_source"]).resolve().parents[2] if "configs" in Path(config["_source"]).parts else Path.cwd()


def run(config_path: str | Path, output_dir: str | Path) -> Path:
    started = time.time()
    config = load_config(config_path)
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    resolved_hash = config_hash(config)
    resolved_path = output / "resolved_config.json"
    write_resolved_config(config, resolved_path)
    raw = run_core(config)

    phyto_size = log_size_grid(1.0, 200.0, 70)
    zoo_size = log_size_grid(2.0, 1000.0, 70)
    traits = phyto_traits(phyto_size, config["community"]["phyto_sinking_coefficient_m_d"])
    cell_number, carbon_biomass, _ = carbon_from_phyto(raw["P"], raw["quota_P"], traits["carbon_pg_cell"])
    net_flux = (
        raw["phyto_growth_rate"] * raw["light_limitation"] * raw["P"]
        - raw["grazing_loss"]
    )
    net_growth_rate = np.divide(
        net_flux,
        raw["P"],
        out=np.full_like(net_flux, np.nan),
        where=raw["P"] > 1e-20,
    )
    potential_flux = traits["mu_max_d"][None, None, :] * raw["P"]
    total_growth_reduction = 1.0 - np.divide(
        net_flux,
        potential_flux,
        out=np.full_like(net_flux, np.nan),
        where=potential_flux > 1e-20,
    )
    entropy, evenness, _, diversity_valid = shannon_metrics(
        raw["P"],
        absolute_threshold=float(config.get("analysis", {}).get("diversity_absolute_threshold", 1e-12)),
        relative_threshold=float(config.get("analysis", {}).get("diversity_relative_threshold", 1e-8)),
        community_threshold=float(config.get("analysis", {}).get("diversity_biomass_threshold", 1e-10)),
    )
    numerics = config["numerics"]
    time_days = np.arange(raw["P"].shape[0], dtype=float) * float(numerics["output_interval_days"])
    depth = (np.arange(int(config["grid"]["layers"]), dtype=float) + 0.5) * float(config["grid"]["dz_m"])
    coords = {
        "time": time_days,
        "depth": depth,
        "phyto_size": phyto_size,
        "zoo_size": zoo_size,
        "nutrient": np.arange(3, dtype=float),
    }
    data = {
        "P": raw["P"],
        "Z": raw["Z"],
        "quota_P": raw["quota_P"],
        "cell_number_P": cell_number,
        "nutrient_pool": raw["nutrient"],
        "detritus": raw["detritus"],
        "surface_irradiance": raw["surface_irradiance"],
        "irradiance": raw["irradiance"],
        "nutrient_limitation": raw["nutrient_limitation"],
        "grazing_limitation": raw["grazing_limitation"],
        "light_limitation": raw["light_limitation"],
        "phyto_growth_rate": raw["phyto_growth_rate"],
        "grazing_loss": raw["grazing_loss"],
        "net_phyto_growth_rate": net_growth_rate,
        "total_growth_reduction": total_growth_reduction,
        "phyto_biomass": raw["P"].sum(axis=-1),
        "zoo_biomass": raw["Z"].sum(axis=-1),
        "shannon_entropy_P": entropy,
        "shannon_evenness_P": evenness,
        "diversity_valid_P": diversity_valid.astype(float),
        "carbon_per_cell_P": traits["carbon_pg_cell"],
        "carbon_biomass_P": carbon_biomass,
        "export_phyto": raw["export_phyto"],
        "export_detritus": raw["export_detritus"],
        "export_total": raw["export_total"],
        "mass_balance_error": raw["mass_balance_error"],
        "nonnegative_correction": raw["nonnegative_correction"],
    }
    root = _repository_root(config)
    source_hash = core_source_hash(root)
    metadata = {
        "schema_version": 1,
        "splas_version": __version__,
        "config_hash": resolved_hash,
        "fortran_source_sha256": source_hash,
        "git_commit": _command_output(["git", "-C", str(root), "rev-parse", "HEAD"]),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "compiler": _command_output(["gfortran", "--version"]),
        "seed": int(numerics["seed"]),
        "elapsed_seconds": time.time() - started,
        "backend": "Fortran f2py",
    }
    write_dataset(output / "results.nc", coords, data, canonical_config(config), metadata)
    (output / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    run_manifest = {
        "name": config["name"],
        "config": str(resolved_path),
        "config_hash": resolved_hash,
        "dataset": str(output / "results.nc"),
        "complete": True,
    }
    (output / "manifest.json").write_text(
        json.dumps(run_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return output


def _run_item(item: tuple[str, str]) -> str:
    config, output = item
    run(config, output)
    return output


def run_ensemble(manifest_path: str | Path, jobs: int, output_dir: str | Path) -> list[Path]:
    manifest_source = Path(manifest_path).resolve()
    manifest = json.loads(manifest_source.read_text(encoding="utf-8"))
    base_output = Path(output_dir).resolve()
    base_output.mkdir(parents=True, exist_ok=True)
    items: list[tuple[str, str]] = []
    for entry in manifest["runs"]:
        config = (manifest_source.parent / entry["config"]).resolve()
        target = base_output / entry["name"]
        if (target / "manifest.json").exists():
            existing = json.loads((target / "manifest.json").read_text(encoding="utf-8"))
            if existing.get("complete") and existing.get("config_hash") == config_hash(load_config(config)):
                continue
        items.append((str(config), str(target)))
    if jobs <= 1:
        completed = [_run_item(item) for item in items]
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as executor:
            completed = list(executor.map(_run_item, items))
    return [Path(path) for path in completed]
