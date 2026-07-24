"""Non-interactive figures for individual runs and the thesis reproduction suite."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .analysis import critical_depth, euphotic_depth, weighted_mean_size
from .dataset import read_dataset


def _finish(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _load(run_dir: Path) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    return read_dataset(run_dir / "results.nc")


def make_run_figures(results_dir: str | Path, output_dir: str | Path) -> list[Path]:
    arrays, attrs = _load(Path(results_dir))
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    created: list[Path] = []
    time, depth = arrays["time"], arrays["depth"]
    ptotal, ztotal = arrays["P"].sum(axis=-1), arrays["Z"].sum(axis=-1)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=True)
    axes[0].plot(arrays["nutrient_pool"][-1, :, 0], depth, label="N")
    axes[0].plot(arrays["irradiance"][-1] / max(arrays["irradiance"][-1, 0], 1e-30), depth, label="I/I(surface)")
    axes[0].set_xlabel("normalized chemical/physical state")
    axes[0].legend()
    axes[1].plot(ptotal[-1], depth, label="P")
    axes[1].plot(ztotal[-1], depth, label="Z")
    axes[1].set_xlabel("biomass (uM N)")
    axes[1].legend()
    axes[0].set_ylabel("depth (m)")
    axes[0].invert_yaxis()
    path = output / "profiles.png"
    _finish(fig, path)
    created.append(path)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    pmean = arrays["P"].mean(axis=0).T
    zmean = arrays["Z"].mean(axis=0).T
    axes[0].pcolormesh(depth, arrays["phyto_size"], pmean, shading="auto")
    axes[1].pcolormesh(depth, arrays["zoo_size"], zmean, shading="auto")
    for ax, title in zip(axes, ("phytoplankton", "zooplankton")):
        ax.set_yscale("log")
        ax.set_xlabel("depth (m)")
        ax.set_title(title)
    axes[0].set_ylabel("ESD (um)")
    path = output / "size_depth_spectra.png"
    _finish(fig, path)
    created.append(path)

    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    axes[0].plot(time, np.trapezoid(ptotal, depth, axis=1), label="P")
    axes[0].plot(time, np.trapezoid(ztotal, depth, axis=1), label="Z")
    axes[0].set_ylabel("integrated biomass")
    axes[0].legend()
    axes[1].plot(time, arrays["surface_irradiance"], color="goldenrod", label="surface light")
    axes[1].plot(time, arrays["export_total"], color="tab:brown", label="C export")
    axes[1].set_xlabel("time (days)")
    axes[1].legend()
    path = output / "time_series.png"
    _finish(fig, path)
    created.append(path)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
    for ax, variable, title in zip(
        axes,
        ("nutrient_limitation", "grazing_limitation", "light_limitation"),
        ("nutrient limitation", "grazing limitation", "light growth factor"),
    ):
        field = arrays[variable].mean(axis=(0, 2))
        ax.plot(field, depth)
        ax.set_title(title)
        ax.set_xlabel("fraction")
    axes[0].set_ylabel("depth (m)")
    axes[0].invert_yaxis()
    path = output / "limitations.png"
    _finish(fig, path)
    created.append(path)

    dcrit, growth = critical_depth(ptotal, time, depth)
    eup = euphotic_depth(arrays["irradiance"], arrays["surface_irradiance"], depth)
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(time, dcrit, label="critical depth")
    ax.plot(time, eup, label="euphotic depth / MLD*")
    ax.set_xlabel("time (days)")
    ax.set_ylabel("depth (m)")
    ax.invert_yaxis()
    ax.legend()
    path = output / "critical_depth.png"
    _finish(fig, path)
    created.append(path)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    axes[0].plot(time, np.nanmean(arrays["shannon_evenness_P"], axis=1))
    axes[0].set_xlabel("time (days)")
    axes[0].set_ylabel("Shannon evenness")
    axes[1].scatter(arrays["surface_irradiance"], arrays["export_total"], s=5, alpha=0.5)
    axes[1].set_xlabel("surface irradiance (W m-2)")
    axes[1].set_ylabel("export (mg C m-2 d-1)")
    path = output / "diversity_export.png"
    _finish(fig, path)
    created.append(path)
    return created


def _ensemble(group_dir: Path) -> list[tuple[dict[str, np.ndarray], dict[str, Any]]]:
    runs = []
    if not group_dir.exists():
        return runs
    for dataset in sorted(group_dir.glob("*/results.nc")):
        runs.append(read_dataset(dataset))
    return sorted(runs, key=lambda item: float(item[1]["config"]["forcing"]["multiplier"]))


def _intensity(item: tuple[dict[str, np.ndarray], dict[str, Any]]) -> float:
    arrays, _ = item
    return float(np.mean(arrays["surface_irradiance"]))


def _metric_plot(
    runs: list[tuple[dict[str, np.ndarray], dict[str, Any]]],
    metric,
    title: str,
    ylabel: str,
    path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    x = [_intensity(item) for item in runs]
    y = [metric(item[0]) for item in runs]
    ax.plot(x, y, "o-")
    ax.set_title(title)
    ax.set_xlabel("mean surface irradiance (W m-2 PAR)")
    ax.set_ylabel(ylabel)
    _finish(fig, path)


def _timeseries_plot(
    runs: list[tuple[dict[str, np.ndarray], dict[str, Any]]],
    variable,
    title: str,
    ylabel: str,
    path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 4))
    for arrays, attrs in runs:
        label = f"{100 * (attrs['config']['forcing']['multiplier'] - 1):+.0f}%"
        ax.plot(arrays["time"], variable(arrays), label=label, linewidth=0.9)
    ax.set_title(title)
    ax.set_xlabel("time (days)")
    ax.set_ylabel(ylabel)
    ax.legend(ncol=2, fontsize=8)
    _finish(fig, path)


def make_thesis_figures(reproduction_dir: str | Path) -> list[Path]:
    root = Path(reproduction_dir)
    output = root / "figures"
    output.mkdir(parents=True, exist_ok=True)
    created: list[Path] = []
    reference = _ensemble(root / "results" / "reference")
    constant = _ensemble(root / "results" / "constant_open")
    seasonal = _ensemble(root / "results" / "seasonal_200m")
    critical = _ensemble(root / "results" / "critical_depth_100m")

    if len(reference) >= 2:
        fig, axes = plt.subplots(2, 2, figsize=(9, 8), sharey=True)
        for row, (arrays, attrs) in enumerate(reference[:2]):
            depth = arrays["depth"]
            axes[row, 0].plot(arrays["nutrient_pool"][-1, :, 0], depth, label="N")
            axes[row, 0].plot(arrays["detritus"][-1], depth, label="OM")
            axes[row, 1].plot(arrays["P"][-1].sum(axis=-1), depth, label="P")
            axes[row, 1].plot(arrays["Z"][-1].sum(axis=-1), depth, label="Z")
            axes[row, 0].set_title(attrs["config"]["system"]["mode"])
            axes[row, 0].legend()
            axes[row, 1].legend()
        axes[0, 0].invert_yaxis()
        for ax in axes.ravel():
            ax.set_ylabel("depth (m)")
        path = output / "fig_3_1_profiles.png"
        _finish(fig, path)
        created.append(path)

        for number, variable, title in [
            ("3_2", "P", "mean size-depth phytoplankton spectrum"),
            ("3_3", "shannon_evenness_P", "Shannon evenness profiles"),
            ("3_4", "carbon_biomass_P", "phytoplankton C:N profiles"),
            ("3_5", "export_total", "biomass and carbon export in time"),
        ]:
            fig, axes = plt.subplots(1, 2, figsize=(10, 4))
            for ax, (arrays, attrs) in zip(axes, reference[:2]):
                if variable == "P":
                    field = arrays["P"].mean(axis=0).T
                    ax.pcolormesh(arrays["depth"], arrays["phyto_size"], field, shading="auto")
                    ax.set_yscale("log")
                    ax.set_ylabel("ESD (um)")
                elif variable == "shannon_evenness_P":
                    ax.plot(np.nanmean(arrays[variable], axis=0), arrays["depth"])
                    ax.invert_yaxis()
                    ax.set_ylabel("depth (m)")
                elif variable == "carbon_biomass_P":
                    carbon = arrays[variable].sum(axis=-1)
                    nitrogen = arrays["P"].sum(axis=-1)
                    ratio = np.divide(carbon, nitrogen, out=np.full_like(carbon, np.nan), where=nitrogen > 0)
                    ax.plot(np.nanmean(ratio, axis=0), arrays["depth"])
                    ax.invert_yaxis()
                    ax.set_ylabel("depth (m)")
                else:
                    ax.plot(arrays["time"], arrays["export_total"], label="export")
                    ax2 = ax.twinx()
                    ax2.plot(arrays["time"], np.trapezoid(arrays["P"].sum(axis=-1), arrays["depth"], axis=1), color="tab:green", alpha=0.7)
                ax.set_title(attrs["config"]["system"]["mode"])
            path = output / f"fig_{number}_{title.replace(' ', '_')}.png"
            _finish(fig, path)
            created.append(path)

    def integrated_p(arrays):
        return float(np.mean(np.trapezoid(arrays["P"].sum(axis=-1), arrays["depth"], axis=1)))

    if constant:
        constant_specs = [
            ("4_2", lambda a: np.trapezoid(a["P"].sum(axis=-1), a["depth"], axis=1), "phytoplankton biomass", "uM N m", "time"),
            ("4_3", lambda a: integrated_p(a), "mean phytoplankton biomass", "uM N m", "metric"),
            ("4_4", lambda a: float(np.mean(a["Z"].sum() / max(a["P"].sum(), 1e-30))), "Z:P / TTE proxy", "1", "metric"),
            ("4_5", lambda a: float(np.nanmean(weighted_mean_size(a["P"], a["phyto_size"]))), "mean phytoplankton size", "um ESD", "metric"),
            ("4_6", lambda a: float(np.nanmean(weighted_mean_size(a["P"], a["phyto_size"]) ** 3)), "mean cell volume proxy", "um3", "metric"),
            ("4_7", lambda a: np.nanmean(weighted_mean_size(a["P"], a["phyto_size"]), axis=1), "size structure in time", "um ESD", "time"),
            ("4_8", lambda a: float(np.nanmean(1.0 - a["light_limitation"])), "mean light limitation", "1", "metric"),
            ("4_9", lambda a: float(np.nanmean(a["total_growth_reduction"])), "total growth reduction", "fraction of mu_max", "metric"),
            ("4_10", lambda a: float(np.mean(a["nutrient_pool"][..., 0])), "mean dissolved nitrogen", "uM N", "metric"),
            ("4_11", lambda a: float(np.nanmean(a["shannon_evenness_P"])), "mean diversity", "Shannon evenness", "metric"),
            ("4_12", lambda a: float(np.mean(a["export_total"])), "mean carbon export", "mg C m-2 d-1", "metric"),
        ]
        for number, metric, title, ylabel, kind in constant_specs:
            path = output / f"fig_{number}_{title.replace(' ', '_').replace('/', '_')}.png"
            if kind == "time":
                _timeseries_plot(constant, metric, title, ylabel, path)
            else:
                _metric_plot(constant, metric, title, ylabel, path)
            created.append(path)

    if seasonal:
        seasonal_specs = [
            ("4_13", lambda a: a["surface_irradiance"], "seasonal surface light", "W m-2 PAR"),
            ("4_14", lambda a: np.nanmean(weighted_mean_size(a["P"], a["phyto_size"]), axis=1), "seasonal mean size", "um ESD"),
            ("4_15", lambda a: np.trapezoid(a["P"].sum(axis=-1), a["depth"], axis=1), "seasonal phytoplankton biomass", "uM N m"),
            ("4_16", lambda a: np.trapezoid(a["P"].sum(axis=-1), a["depth"], axis=1), "biomass regime transitions", "uM N m"),
            ("4_19", lambda a: np.nanmean(a["shannon_evenness_P"], axis=1), "diversity in time", "Shannon evenness"),
            ("4_20", lambda a: np.nansum(a["carbon_biomass_P"], axis=(1, 2)) / np.maximum(np.nansum(a["P"], axis=(1, 2)), 1e-30), "community C:N", "mol C mol N-1"),
            ("4_21", lambda a: a["export_total"], "carbon export in time", "mg C m-2 d-1"),
            ("4_22", lambda a: a["export_total"] / np.maximum(a["surface_irradiance"], 1e-30), "export efficiency per light", "mg C W-1"),
        ]
        for number, variable, title, ylabel in seasonal_specs:
            path = output / f"fig_{number}_{title.replace(' ', '_').replace(':', '')}.png"
            _timeseries_plot(seasonal, variable, title, ylabel, path)
            created.append(path)
        for number, metric, title, ylabel in [
            ("4_17", lambda a: float(np.nanmean(critical_depth(a["P"].sum(axis=-1), a["time"], a["depth"])[0])), "mean critical depth", "m"),
            ("4_18", lambda a: float(np.nanmean(euphotic_depth(a["irradiance"], a["surface_irradiance"], a["depth"]))), "mean euphotic depth / MLD*", "m"),
        ]:
            path = output / f"fig_{number}_{title.replace(' ', '_').replace('/', '_')}.png"
            _metric_plot(seasonal, metric, title, ylabel, path)
            created.append(path)

    if critical:
        critical_specs = [
            ("5_2", lambda a: np.trapezoid(a["P"].sum(axis=-1), a["depth"], axis=1), "100m biomass", "uM N m"),
            ("5_3", lambda a: np.nanmean(1.0 - a["light_limitation"], axis=(1, 2)), "100m light limitation", "1"),
            ("5_4", lambda a: critical_depth(a["P"].sum(axis=-1), a["time"], a["depth"])[0], "critical-depth crossing", "m"),
            ("5_5", lambda a: euphotic_depth(a["irradiance"], a["surface_irradiance"], a["depth"]), "euphotic-depth crossing", "m"),
            ("5_6", lambda a: a["export_total"], "100m carbon export", "mg C m-2 d-1"),
            ("5_7", lambda a: a["export_total"] / np.maximum(a["surface_irradiance"], 1e-30), "100m export efficiency", "mg C W-1"),
        ]
        for number, variable, title, ylabel in critical_specs:
            path = output / f"fig_{number}_{title.replace(' ', '_')}.png"
            _timeseries_plot(critical, variable, title, ylabel, path)
            created.append(path)

    return created


def write_figure_manifest(path: str | Path, created: list[Path], status: str = "generated_unverified") -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["thesis_figure", "file", "source_run", "metrics", "status", "explanation"])
        for figure in created:
            stem = figure.stem
            number = stem.split("_", 3)[1:3]
            writer.writerow([
                ".".join(number),
                str(figure),
                "reproduction/results",
                stem,
                status,
                "Generated from reconstructed SPLAS outputs; original numerical arrays are unavailable.",
            ])
