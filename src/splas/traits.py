"""Pure size-grid and allometric functions used for validation and analysis."""

from __future__ import annotations

import numpy as np


def log_size_grid(minimum_um: float, maximum_um: float, classes: int = 70) -> np.ndarray:
    if minimum_um <= 0 or maximum_um <= minimum_um or classes < 2:
        raise ValueError("invalid logarithmic size-grid bounds")
    return np.geomspace(minimum_um, maximum_um, classes)


def phyto_traits(size_um: np.ndarray, sinking_coefficient: float = 0.0041) -> dict[str, np.ndarray]:
    size = np.asarray(size_um, dtype=float)
    volume = (np.pi / 6.0) * size**3
    return {
        "cell_volume_um3": volume,
        "mu_max_d": 2.6 * size**-0.45,
        "vmax_uM_N_cell_d": 9.10e-12 * volume**0.67,
        "qmin_uM_N_cell": 1.36e-12 * volume**0.77,
        "half_saturation_uM_N": 0.1 * size,
        "sinking_m_d": sinking_coefficient * size**1.17,
        "carbon_pg_cell": 0.0109 * volume**0.991,
    }


def zoo_traits(size_um: np.ndarray) -> dict[str, np.ndarray]:
    size = np.asarray(size_um, dtype=float)
    return {
        "optimal_prey_um": 0.16 * size * np.exp(-0.02 * np.log(size) ** 2),
        "ingestion_d": 26.0 * size**-0.4,
    }


def preference_matrix(
    phyto_size_um: np.ndarray,
    optimal_prey_um: np.ndarray,
    tolerance_log10_um: float = 0.25,
) -> np.ndarray:
    phyto = np.asarray(phyto_size_um, dtype=float)[:, None]
    optimum = np.asarray(optimal_prey_um, dtype=float)[None, :]
    return np.exp(-((np.log10(phyto) - np.log10(optimum)) / tolerance_log10_um) ** 2)


def monod(nutrient: np.ndarray | float, half_saturation: np.ndarray) -> np.ndarray:
    nutrient_nonnegative = np.maximum(np.asarray(nutrient, dtype=float), 0.0)
    return nutrient_nonnegative / (nutrient_nonnegative + np.asarray(half_saturation))


def droop(mu_max: np.ndarray, qmin: np.ndarray, quota: np.ndarray) -> np.ndarray:
    q = np.maximum(np.asarray(quota, dtype=float), np.asarray(qmin, dtype=float))
    return np.asarray(mu_max, dtype=float) * (1.0 - np.asarray(qmin, dtype=float) / q)
