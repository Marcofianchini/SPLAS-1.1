"""Configuration loading, inheritance, validation, and hashing."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from typing import Any


class ConfigError(ValueError):
    """Raised when a SPLAS configuration is incomplete or inconsistent."""


def _merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(base)
    for key, value in override.items():
        if key == "extends":
            continue
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def load_config(path: str | Path) -> dict[str, Any]:
    """Load a JSON config, recursively resolving ``extends`` relative to the file."""

    source = Path(path).resolve()
    try:
        data = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ConfigError(f"Cannot read configuration {source}: {exc}") from exc
    if not isinstance(data, dict):
        raise ConfigError(f"Configuration {source} must contain a JSON object")
    if "extends" in data:
        parent = (source.parent / data["extends"]).resolve()
        data = _merge(load_config(parent), data)
    data["_source"] = str(source)
    validate_config(data)
    return data


def _need(mapping: dict[str, Any], path: str) -> Any:
    current: Any = mapping
    for part in path.split("."):
        if not isinstance(current, dict) or part not in current:
            raise ConfigError(f"Missing required configuration field: {path}")
        current = current[part]
    return current


def validate_config(config: dict[str, Any]) -> None:
    """Validate the fixed 70 P / 70 Z / 3-pool SPLAS contract."""

    required = [
        "name",
        "scientific_profile",
        "system.mode",
        "system.nutrient_supply_uM_N_d",
        "system.diffusion_m2_s",
        "system.remineralization_d",
        "system.detritus_sinking_m_d",
        "system.mortality_recycling_fraction",
        "grid.layers",
        "grid.dz_m",
        "community.phyto_classes",
        "community.zoo_classes",
        "community.nutrient_pools",
        "community.phyto_size_min_um",
        "community.phyto_size_max_um",
        "community.zoo_size_min_um",
        "community.zoo_size_max_um",
        "community.initial_phyto_uM_N",
        "community.initial_zoo_uM_N",
        "community.initial_nutrient_uM_N",
        "community.phyto_mortality_d",
        "community.phyto_sinking_coefficient_m_d",
        "community.zoo_growth_efficiency",
        "community.zoo_egestion_fraction",
        "community.zoo_mortality_uM_N_d",
        "community.zoo_half_saturation_uM_N",
        "community.prey_tolerance_log10_um",
        "optics.kw_m",
        "optics.kp_m_uM_N",
        "optics.theta_chl_c",
        "optics.alpha_m2_W",
        "forcing.mode",
        "forcing.reference_irradiance_W_m2",
        "forcing.multiplier",
        "numerics.duration_days",
        "numerics.timestep_days",
        "numerics.output_interval_days",
        "numerics.seed",
    ]
    for field in required:
        _need(config, field)
    if config["scientific_profile"] not in {"thesis_consistent", "legacy_reconstructed"}:
        raise ConfigError("scientific_profile must be thesis_consistent or legacy_reconstructed")
    if config["system"]["mode"] not in {"open", "closed"}:
        raise ConfigError("system.mode must be open or closed")
    if config["forcing"]["mode"] not in {"constant", "seasonal", "daily", "linear", "normal"}:
        raise ConfigError("forcing.mode must be constant, seasonal, daily, linear, or normal")
    expected = {"phyto_classes": 70, "zoo_classes": 70, "nutrient_pools": 3}
    for key, value in expected.items():
        if int(config["community"][key]) != value:
            raise ConfigError(f"{key} is fixed at {value} by the current Fortran core")
    fixed_ranges = {
        "phyto_size_min_um": 1.0,
        "phyto_size_max_um": 200.0,
        "zoo_size_min_um": 2.0,
        "zoo_size_max_um": 1000.0,
    }
    for key, value in fixed_ranges.items():
        if float(config["community"][key]) != value:
            raise ConfigError(f"{key} is fixed at {value} by the current Fortran core")
    layers = int(config["grid"]["layers"])
    if layers < 3:
        raise ConfigError("grid.layers must be at least 3")
    for field in ("duration_days", "timestep_days", "output_interval_days"):
        if float(config["numerics"][field]) <= 0:
            raise ConfigError(f"numerics.{field} must be positive")
    dt = float(config["numerics"]["timestep_days"])
    output = float(config["numerics"]["output_interval_days"])
    duration = float(config["numerics"]["duration_days"])
    steps_per_output = output / dt
    total_steps = duration / dt
    if abs(steps_per_output - round(steps_per_output)) > 1e-9:
        raise ConfigError("output_interval_days must be an integer multiple of timestep_days")
    if abs(total_steps - round(total_steps)) > 1e-9:
        raise ConfigError("duration_days must be an integer multiple of timestep_days")
    if int(round(total_steps)) % int(round(steps_per_output)):
        raise ConfigError("duration_days must contain a whole number of output intervals")
    fractions = (
        float(config["community"]["zoo_growth_efficiency"]),
        float(config["community"]["zoo_egestion_fraction"]),
    )
    if sum(fractions) > 1.0:
        raise ConfigError("zoo growth efficiency plus egestion fraction cannot exceed one")
    if config["system"]["mode"] == "closed":
        nonzero = [
            ("nutrient_supply_uM_N_d", config["system"]["nutrient_supply_uM_N_d"]),
            ("diffusion_m2_s", config["system"]["diffusion_m2_s"]),
            ("phyto_sinking_coefficient_m_d", config["community"]["phyto_sinking_coefficient_m_d"]),
            ("detritus_sinking_m_d", config["system"]["detritus_sinking_m_d"]),
        ]
        if any(abs(float(value)) > 0 for _, value in nonzero):
            raise ConfigError("closed mode requires zero supply, diffusion, and sinking")


def canonical_config(config: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in config.items() if not key.startswith("_")}


def config_hash(config: dict[str, Any]) -> str:
    payload = json.dumps(canonical_config(config), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def write_resolved_config(config: dict[str, Any], path: str | Path) -> None:
    Path(path).write_text(
        json.dumps(canonical_config(config), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
