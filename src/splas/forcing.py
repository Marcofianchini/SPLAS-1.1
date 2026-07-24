"""Deterministic light-forcing functions mirrored by the Fortran core."""

from __future__ import annotations

import numpy as np


def surface_irradiance(
    time_days: np.ndarray | float,
    mode: str,
    reference: float = 959.825439,
    multiplier: float = 1.0,
    total_delta: float = 0.0,
    sd: float = 0.0,
    seed: int = 0,
    duration_days: float | None = None,
) -> np.ndarray:
    time = np.asarray(time_days, dtype=float)
    baseline = reference * multiplier
    if mode == "constant":
        return np.full_like(time, baseline)
    if mode == "seasonal":
        return baseline * (0.5 + 0.5 * np.sin(2.0 * np.pi * time / 365.0) ** 2)
    if mode == "daily":
        return baseline * np.sin(np.pi * time) ** 2
    if mode == "linear":
        duration = float(duration_days if duration_days is not None else np.max(time))
        fraction = np.zeros_like(time) if duration == 0 else time / duration
        return baseline + total_delta * fraction
    if mode == "normal":
        state = max(1, abs(int(seed)))
        draws = np.empty(time.size, dtype=float)
        for index in range(time.size):
            state = (16807 * state) % 2147483647
            u1 = state / 2147483647.0
            state = (16807 * state) % 2147483647
            u2 = state / 2147483647.0
            radius = np.sqrt(-2.0 * np.log(max(u1, np.finfo(float).tiny)))
            draws[index] = radius * np.cos(2.0 * np.pi * u2)
        return np.maximum(0.0, baseline + sd * draws.reshape(time.shape))
    raise ValueError(f"Unknown forcing mode: {mode}")


def attenuated_irradiance(
    surface: float,
    depth_m: np.ndarray,
    phyto_biomass_by_depth: np.ndarray,
    dz_m: float,
    kw_m: float = 0.03,
    kp_m_uM_N: float = 0.01,
) -> np.ndarray:
    depth = np.asarray(depth_m, dtype=float)
    biomass = np.asarray(phyto_biomass_by_depth, dtype=float)
    if depth.shape != biomass.shape:
        raise ValueError("depth and phyto_biomass_by_depth must have the same shape")
    cumulative = np.cumsum(np.maximum(biomass, 0.0)) * float(dz_m)
    return float(surface) * np.exp(-float(kw_m) * depth - float(kp_m_uM_N) * cumulative)


def light_limitation(
    irradiance_W_m2: np.ndarray | float,
    mu_max_d: np.ndarray,
    theta_chl_c: float = 0.02,
    alpha_m2_W: float = 0.32,
) -> np.ndarray:
    irradiance = np.asarray(irradiance_W_m2, dtype=float)
    mu = np.asarray(mu_max_d, dtype=float)
    return 1.0 - np.exp(
        -np.expand_dims(irradiance, -1) * theta_chl_c * alpha_m2_W / mu
    )
