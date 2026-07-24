"""Thin, explicit bridge from validated Python configuration to the Fortran core."""

from __future__ import annotations

import hashlib
import importlib
from pathlib import Path
from typing import Any

import numpy as np


class CoreUnavailable(RuntimeError):
    """Raised when the compiled Fortran extension has not been built."""


def require_core() -> Any:
    try:
        module = importlib.import_module("splas._core")
    except ImportError as exc:
        raise CoreUnavailable(
            "The SPLAS Fortran core is unavailable. Run `make build` with gfortran installed."
        ) from exc
    return getattr(module, "splas_core", module)


def core_source_hash(repository_root: str | Path) -> str:
    source = Path(repository_root) / "src" / "fortran" / "splas_core.f90"
    return hashlib.sha256(source.read_bytes()).hexdigest()


def parameter_vector(config: dict[str, Any]) -> np.ndarray:
    forcing = config["forcing"]
    system = config["system"]
    community = config["community"]
    optics = config["optics"]
    values = [
        forcing["reference_irradiance_W_m2"],
        forcing["multiplier"],
        forcing.get("total_delta_W_m2", 0.0),
        forcing.get("normal_sd_W_m2", 0.0),
        system["nutrient_supply_uM_N_d"],
        system["diffusion_m2_s"],
        community["phyto_sinking_coefficient_m_d"],
        system["detritus_sinking_m_d"],
        system["remineralization_d"],
        system["mortality_recycling_fraction"],
        optics["kw_m"],
        optics["kp_m_uM_N"],
        optics["theta_chl_c"],
        optics["alpha_m2_W"],
        community["phyto_mortality_d"],
        community["zoo_growth_efficiency"],
        community["zoo_egestion_fraction"],
        community["zoo_mortality_uM_N_d"],
        community["zoo_half_saturation_uM_N"],
        community["prey_tolerance_log10_um"],
        community["initial_phyto_uM_N"],
        community["initial_zoo_uM_N"],
        community["initial_nutrient_uM_N"],
        community.get("initial_detritus_uM_N", 0.0),
        community.get("initial_quota_qmin_multiplier", 1.1),
        system.get("interior_supply_fraction", 0.1),
    ]
    return np.asarray(values, dtype=np.float64)


def run_core(config: dict[str, Any]) -> dict[str, np.ndarray]:
    api = require_core()
    numerics = config["numerics"]
    grid = config["grid"]
    dt = float(numerics["timestep_days"])
    duration = float(numerics["duration_days"])
    output_interval = float(numerics["output_interval_days"])
    n_steps = int(round(duration / dt))
    output_every = int(round(output_interval / dt))
    nout = n_steps // output_every + 1
    mode_codes = {"closed": 0, "open": 1}
    forcing_codes = {"constant": 0, "seasonal": 1, "daily": 2, "linear": 3, "normal": 4}
    profile_codes = {"thesis_consistent": 0, "legacy_reconstructed": 1}
    outputs = api.run_model(
        n_steps,
        output_every,
        nout,
        int(grid["layers"]),
        dt,
        float(grid["dz_m"]),
        mode_codes[config["system"]["mode"]],
        forcing_codes[config["forcing"]["mode"]],
        profile_codes[config["scientific_profile"]],
        int(numerics["seed"]),
        parameter_vector(config),
    )
    names = [
        "P_raw",
        "Z_raw",
        "quota_raw",
        "nutrient_raw",
        "detritus",
        "surface_irradiance",
        "irradiance",
        "nutrient_limitation_raw",
        "grazing_limitation_raw",
        "light_limitation_raw",
        "phyto_growth_rate_raw",
        "grazing_loss_raw",
        "export_phyto",
        "export_detritus",
        "mass_balance_error",
        "nonnegative_correction",
        "status",
    ]
    if len(outputs) != len(names):
        raise RuntimeError(f"Unexpected Fortran output count: {len(outputs)} (expected {len(names)})")
    result = dict(zip(names, outputs, strict=True))
    status = int(np.asarray(result.pop("status")).item())
    if status != 0:
        meanings = {
            10: "invalid dimensions",
            11: "inconsistent output cadence",
            12: "non-positive grid step",
            20: "closed-system mass-balance tolerance exceeded",
            30: "non-finite state",
            31: "negative state beyond roundoff tolerance",
        }
        raise RuntimeError(f"Fortran core stopped with status {status}: {meanings.get(status, 'unknown')}")
    for raw_name, final_name in [
        ("P_raw", "P"),
        ("Z_raw", "Z"),
        ("quota_raw", "quota_P"),
        ("nutrient_raw", "nutrient"),
        ("nutrient_limitation_raw", "nutrient_limitation"),
        ("grazing_limitation_raw", "grazing_limitation"),
        ("light_limitation_raw", "light_limitation"),
        ("phyto_growth_rate_raw", "phyto_growth_rate"),
        ("grazing_loss_raw", "grazing_loss"),
    ]:
        result[final_name] = np.transpose(np.asarray(result.pop(raw_name)), (1, 2, 0))
    for name, values in list(result.items()):
        result[name] = np.asarray(values, dtype=np.float64)
    result["export_total"] = result["export_phyto"] + result["export_detritus"]
    return result
