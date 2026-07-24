"""Pure scientific post-processing for SPLAS results."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .dataset import read_dataset


def shannon_metrics(
    biomass: np.ndarray,
    absolute_threshold: float = 1e-12,
    relative_threshold: float = 1e-8,
    community_threshold: float = 1e-10,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    values = np.maximum(np.asarray(biomass, dtype=float), 0.0)
    total = values.sum(axis=-1)
    cutoff = np.maximum(absolute_threshold, relative_threshold * total[..., None])
    present = values > cutoff
    count = present.sum(axis=-1)
    safe_total = np.where(total > 0.0, total, 1.0)
    proportions = np.where(present, values / safe_total[..., None], 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        entropy = -np.sum(np.where(proportions > 0.0, proportions * np.log(proportions), 0.0), axis=-1)
    normalizer = np.where(count > 1, np.log(np.maximum(count, 2)), 1.0)
    evenness = np.where(count > 1, entropy / normalizer, np.where(count == 1, 0.0, np.nan))
    valid = total >= community_threshold
    entropy = np.where(valid, entropy, np.nan)
    evenness = np.where(valid, evenness, np.nan)
    return entropy, evenness, count, valid


def weighted_mean_size(biomass: np.ndarray, size_um: np.ndarray) -> np.ndarray:
    values = np.maximum(np.asarray(biomass, dtype=float), 0.0)
    denominator = values.sum(axis=-1)
    numerator = np.sum(values * np.asarray(size_um, dtype=float), axis=-1)
    return np.divide(numerator, denominator, out=np.full_like(numerator, np.nan), where=denominator > 0)


def biological_maximum_depth(
    profile: np.ndarray,
    depth_m: np.ndarray,
    relative_prominence: float = 0.05,
) -> list[np.ndarray]:
    profiles = np.atleast_2d(np.asarray(profile, dtype=float))
    depth = np.asarray(depth_m, dtype=float)
    maxima: list[np.ndarray] = []
    for values in profiles:
        candidates = []
        peak_floor = relative_prominence * max(float(np.nanmax(values)), 0.0)
        for index in range(1, len(values) - 1):
            local_prominence = values[index] - max(values[index - 1], values[index + 1])
            if values[index] >= values[index - 1] and values[index] >= values[index + 1] and local_prominence >= peak_floor:
                candidates.append(depth[index])
        if not candidates and np.any(np.isfinite(values)):
            candidates = [depth[int(np.nanargmax(values))]]
        maxima.append(np.asarray(candidates, dtype=float))
    return maxima


def euphotic_depth(
    irradiance: np.ndarray,
    surface_irradiance: np.ndarray,
    depth_m: np.ndarray,
    fraction: float = 0.01,
) -> np.ndarray:
    light = np.asarray(irradiance, dtype=float)
    surface = np.asarray(surface_irradiance, dtype=float)
    depth = np.asarray(depth_m, dtype=float)
    result = np.full(light.shape[0], np.nan)
    for index in range(light.shape[0]):
        valid = np.flatnonzero(light[index] > fraction * surface[index])
        if valid.size:
            result[index] = depth[valid[-1]]
    return result


def critical_depth(
    phyto_biomass: np.ndarray,
    time_days: np.ndarray,
    depth_m: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    biomass = np.maximum(np.asarray(phyto_biomass, dtype=float), np.finfo(float).tiny)
    time = np.asarray(time_days, dtype=float)
    depth = np.asarray(depth_m, dtype=float)
    cumulative = np.cumsum(biomass, axis=1)
    growth = np.full_like(cumulative, np.nan)
    growth[1:] = np.diff(np.log(cumulative), axis=0) / np.diff(time)[:, None]
    crossing = np.full(len(time), np.nan)
    for ti in range(1, len(time)):
        values = growth[ti]
        indices = np.flatnonzero(values[:-1] * values[1:] <= 0.0)
        if indices.size:
            k = int(indices[-1])
            y0, y1 = values[k], values[k + 1]
            if y1 == y0:
                crossing[ti] = depth[k]
            else:
                crossing[ti] = depth[k] - y0 * (depth[k + 1] - depth[k]) / (y1 - y0)
    return crossing, growth


def compensation_light(
    critical_depth_m: np.ndarray,
    irradiance: np.ndarray,
    depth_m: np.ndarray,
) -> np.ndarray:
    result = np.full(len(critical_depth_m), np.nan)
    for index, value in enumerate(np.asarray(critical_depth_m, dtype=float)):
        if np.isfinite(value):
            result[index] = np.interp(value, depth_m, irradiance[index])
    return result


def detect_blooms(
    biomass_time_depth: np.ndarray,
    time_days: np.ndarray,
    depth_m: np.ndarray,
    minimum_biomass: float = 1e-3,
) -> list[dict[str, float]]:
    biomass = np.asarray(biomass_time_depth, dtype=float)
    time = np.asarray(time_days, dtype=float)
    depth = np.asarray(depth_m, dtype=float)
    integrated = np.trapezoid(biomass, depth, axis=1)
    rate = np.full_like(integrated, np.nan)
    rate[1:] = np.diff(np.log(np.maximum(integrated, np.finfo(float).tiny))) / np.diff(time)
    events: list[dict[str, float]] = []
    active: int | None = None
    for index in range(1, len(time)):
        onset = rate[index - 1] <= 0.0 < rate[index] and integrated[index] >= minimum_biomass
        if onset and active is None:
            active = index
        if active is not None and (rate[index] <= 0.0 or index == len(time) - 1):
            stop = max(index, active + 1)
            peak = active + int(np.argmax(integrated[active : stop + 1]))
            depth_peak = int(np.argmax(biomass[peak]))
            events.append(
                {
                    "onset_day": float(time[active]),
                    "peak_day": float(time[peak]),
                    "duration_days": float(time[stop] - time[active]),
                    "amplitude": float(integrated[peak] - integrated[active]),
                    "peak_depth_m": float(depth[depth_peak]),
                }
            )
            active = None
    return events


def carbon_from_phyto(
    phyto_uM_N: np.ndarray,
    quota_uM_N_cell: np.ndarray,
    carbon_pg_cell: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    cells = np.divide(
        phyto_uM_N,
        quota_uM_N_cell,
        out=np.zeros_like(phyto_uM_N, dtype=float),
        where=np.asarray(quota_uM_N_cell) > 0,
    )
    carbon_mmol = cells * np.asarray(carbon_pg_cell) * (1e-9 / 12.0)
    # Numerically, 1 uM in water is 1 mmol m-3.
    nitrogen_mmol = np.asarray(phyto_uM_N)
    cn = np.divide(
        carbon_mmol,
        nitrogen_mmol,
        out=np.full_like(carbon_mmol, np.nan),
        where=nitrogen_mmol > 0,
    )
    return cells, carbon_mmol, cn


def export_efficiencies(
    export_total: np.ndarray,
    phyto_biomass: np.ndarray,
    surface_light: np.ndarray,
    time_days: np.ndarray,
    depth_m: np.ndarray,
    transient_days: float = 365.0,
) -> dict[str, float]:
    time = np.asarray(time_days)
    mask = time >= transient_days
    if not np.any(mask):
        mask = np.ones_like(time, dtype=bool)
    export = np.asarray(export_total)[mask]
    biomass = np.trapezoid(np.asarray(phyto_biomass)[mask], depth_m, axis=1)
    light = np.asarray(surface_light)[mask]
    mean_export = float(np.mean(export))
    mean_biomass = float(np.mean(biomass))
    mean_light = float(np.mean(light))
    return {
        "export_mean_mg_C_m2_d": mean_export,
        "export_integral_mg_C_m2": float(np.trapezoid(export, time[mask])),
        "export_per_biomass": mean_export / mean_biomass if mean_biomass else float("nan"),
        "export_per_light": mean_export / mean_light if mean_light else float("nan"),
        "export_per_biomass_light": mean_export / (mean_biomass * mean_light) if mean_biomass and mean_light else float("nan"),
    }


def analyze_results(results_dir: str | Path, output_dir: str | Path) -> dict[str, Any]:
    source = Path(results_dir) / "results.nc"
    arrays, attrs = read_dataset(source)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    depth = arrays["depth"]
    time = arrays["time"]
    p_total = arrays["P"].sum(axis=-1)
    z_total = arrays["Z"].sum(axis=-1)
    euphotic = euphotic_depth(arrays["irradiance"], arrays["surface_irradiance"], depth)
    dcrit, cumulative_growth = critical_depth(p_total, time, depth)
    icrit = compensation_light(dcrit, arrays["irradiance"], depth)
    mean_size = weighted_mean_size(arrays["P"], arrays["phyto_size"])
    blooms = detect_blooms(p_total, time, depth)
    efficiencies = export_efficiencies(
        arrays["export_total"],
        p_total,
        arrays["surface_irradiance"],
        time,
        depth,
        transient_days=float(attrs["config"].get("analysis", {}).get("transient_days", 365.0)),
    )
    summary: dict[str, Any] = {
        "run_name": attrs["config"]["name"],
        "config_hash": attrs["metadata"]["config_hash"],
        "mean_integrated_phyto_uM_N_m": float(np.mean(np.trapezoid(p_total, depth, axis=1))),
        "mean_integrated_zoo_uM_N_m": float(np.mean(np.trapezoid(z_total, depth, axis=1))),
        "mean_euphotic_depth_m": float(np.nanmean(euphotic)),
        "mean_critical_depth_m": float(np.nanmean(dcrit)) if np.any(np.isfinite(dcrit)) else None,
        "mean_compensation_light_W_m2": float(np.nanmean(icrit)) if np.any(np.isfinite(icrit)) else None,
        "bloom_count": len(blooms),
        "blooms": blooms,
        **efficiencies,
    }
    np.savez_compressed(
        output / "analysis.npz",
        euphotic_depth=euphotic,
        critical_depth=dcrit,
        compensation_light=icrit,
        cumulative_growth=cumulative_growth,
        weighted_mean_phyto_size=mean_size,
    )
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with (output / "summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["metric", "value"])
        for key, value in summary.items():
            if not isinstance(value, (dict, list)):
                writer.writerow([key, value])
    return summary
