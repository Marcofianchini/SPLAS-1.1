"""End-to-end thesis suite, ecological regression checks, and report generation."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .analysis import critical_depth, detect_blooms, euphotic_depth
from .dataset import read_dataset
from .figures import make_thesis_figures, write_figure_manifest
from .runner import run_ensemble


def _runs(group: Path) -> list[tuple[dict[str, np.ndarray], dict[str, Any]]]:
    values = [read_dataset(path) for path in group.glob("*/results.nc")]
    return sorted(values, key=lambda item: float(item[1]["config"]["forcing"]["multiplier"]))


def _monotonic(values: list[float], increasing: bool = True, tolerance: float = 0.0) -> bool:
    differences = np.diff(np.asarray(values, dtype=float))
    return bool(np.all(differences >= -tolerance)) if increasing else bool(np.all(differences <= tolerance))


def _mean_integrated_p(arrays: dict[str, np.ndarray]) -> float:
    return float(np.mean(np.trapezoid(arrays["P"].sum(axis=-1), arrays["depth"], axis=1)))


def _dbm(arrays: dict[str, np.ndarray]) -> float:
    profile = arrays["P"].mean(axis=(0, 2))
    return float(arrays["depth"][int(np.argmax(profile))])


def _critical(arrays: dict[str, np.ndarray]) -> float:
    values, _ = critical_depth(arrays["P"].sum(axis=-1), arrays["time"], arrays["depth"])
    return float(np.nanmean(values))


def _compensation(arrays: dict[str, np.ndarray]) -> float:
    dcrit, _ = critical_depth(arrays["P"].sum(axis=-1), arrays["time"], arrays["depth"])
    values = [
        np.interp(depth, arrays["depth"], arrays["irradiance"][index])
        for index, depth in enumerate(dcrit)
        if np.isfinite(depth)
    ]
    return float(np.mean(values)) if values else float("nan")


def _efficiency(arrays: dict[str, np.ndarray]) -> float:
    return float(np.mean(arrays["export_total"] / np.maximum(arrays["surface_irradiance"], 1e-30)))


def _write_summary_table(root: Path, groups: dict[str, list[tuple[dict[str, np.ndarray], dict[str, Any]]]]) -> None:
    path = root / "tables" / "ensemble_summary.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "ensemble",
            "run",
            "light_multiplier",
            "mean_surface_light_W_m2",
            "mean_integrated_P_uM_N_m",
            "dbm_depth_m",
            "mean_export_mg_C_m2_d",
            "export_per_light_mg_C_W",
            "mean_diversity",
            "mean_critical_depth_m",
            "mean_compensation_light_W_m2",
        ])
        for group_name, runs in groups.items():
            for arrays, attrs in runs:
                writer.writerow([
                    group_name,
                    attrs["config"]["name"],
                    attrs["config"]["forcing"]["multiplier"],
                    float(np.mean(arrays["surface_irradiance"])),
                    _mean_integrated_p(arrays),
                    _dbm(arrays),
                    float(np.mean(arrays["export_total"])),
                    _efficiency(arrays),
                    float(np.nanmean(arrays["shannon_evenness_P"])),
                    _critical(arrays) if len(arrays["time"]) > 1 else "",
                    _compensation(arrays) if len(arrays["time"]) > 1 else "",
                ])


def ecological_checks(groups: dict[str, list[tuple[dict[str, np.ndarray], dict[str, Any]]]]) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []

    def add(name: str, expected: str, passed: bool, observed: str) -> None:
        checks.append(
            {
                "name": name,
                "expected": expected,
                "passed": bool(passed),
                "level": "R1" if passed else "R0",
                "observed": observed,
            }
        )

    constant = groups.get("constant_open", [])
    if constant:
        biomass = [_mean_integrated_p(arrays) for arrays, _ in constant]
        dbm = [_dbm(arrays) for arrays, _ in constant]
        export = [float(np.mean(arrays["export_total"])) for arrays, _ in constant]
        efficiency = [_efficiency(arrays) for arrays, _ in constant]
        add("constant_biomass", "biomass increases with I0", _monotonic(biomass, True), repr(biomass))
        add("constant_dbm", "DBM tends to deepen with I0", dbm[-1] >= dbm[0], repr(dbm))
        add("constant_export", "total export increases with I0", _monotonic(export, True), repr(export))
        add("constant_export_efficiency", "export per light decreases with I0", _monotonic(efficiency, False), repr(efficiency))

    closed = groups.get("constant_closed", [])
    if closed:
        diversity = [float(np.nanmean(arrays["shannon_evenness_P"])) for arrays, _ in closed]
        add("closed_diversity", "direct light effect lowers diversity", diversity[-1] <= diversity[0], repr(diversity))

    seasonal = groups.get("seasonal_200m", [])
    if seasonal:
        bloom_counts = [
            len(detect_blooms(arrays["P"].sum(axis=-1), arrays["time"], arrays["depth"]))
            for arrays, _ in seasonal
        ]
        critical_values = [_critical(arrays) for arrays, _ in seasonal]
        compensation = [_compensation(arrays) for arrays, _ in seasonal]
        add("seasonal_blooms", "recurring blooms occur under seasonal forcing", max(bloom_counts) >= 2, repr(bloom_counts))
        add("critical_depth", "mean critical depth decreases with I0", _monotonic(critical_values, False), repr(critical_values))
        add("compensation_light", "compensation light increases with I0", _monotonic(compensation, True), repr(compensation))

    critical = groups.get("critical_depth_100m", [])
    if critical:
        second_bloom = []
        for arrays, _ in critical:
            events = detect_blooms(arrays["P"].sum(axis=-1), arrays["time"], arrays["depth"])
            second_bloom.append(len(events) >= 2)
        add("second_bloom_100m", "a second seasonal peak appears when Dcr crosses 100 m", any(second_bloom), repr(second_bloom))
    return checks


def _write_checks(root: Path, checks: list[dict[str, Any]]) -> None:
    path = root / "tables" / "ecological_checks.csv"
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["name", "expected", "passed", "level", "observed"])
        writer.writeheader()
        writer.writerows(checks)


def _write_report(root: Path, checks: list[dict[str, Any]], figures: list[Path]) -> None:
    passed = sum(check["passed"] for check in checks)
    lines = [
        "# SPLAS 1.1 thesis reproduction report",
        "",
        "## Scope and evidence",
        "",
        "The numerical pipeline was reconstructed from the Fianchini thesis and the supplied SPLAS-1.1 source archive. "
        "No original simulation arrays were supplied, so R2 quantitative equivalence cannot be claimed. "
        "R1 denotes a reproduced qualitative sign, ordering, timing, or pattern; R0 denotes a failed or unverifiable result.",
        "",
        "## Automated ecological checks",
        "",
        "| Check | Expected | Observed | Level |",
        "|---|---|---|---|",
    ]
    for check in checks:
        observed = str(check["observed"]).replace("|", "/")
        lines.append(f"| {check['name']} | {check['expected']} | {observed} | {check['level']} |")
    lines.extend(
        [
            "",
            f"{passed} of {len(checks)} qualitative checks passed.",
            "",
            "## Figure status",
            "",
            f"{len(figures)} thesis-numbered figures were generated from the reconstructed outputs. "
            "They remain R1 only where the associated ecological check passes; visual or pixel-level identity is R0 "
            "because the original arrays and plotting scripts were not supplied.",
            "",
            "## Main unresolved scientific limits",
            "",
            "- Table 3 reports open-system S=1 uM N d-1, whereas section 3.2 states S=4; the recommended profile follows the text and an alternative remains documented.",
            "- Equation 6 uses cellular uptake while equation 4 evolves nitrogen biomass by Droop growth. The conservative thesis profile couples dissolved-N loss to realized P production; the literal historical uptake form remains available as legacy_reconstructed.",
            "- Initial cellular quota is not specified by the thesis. The reconstruction starts at 1.1 times Qmin and records that choice.",
            "- The thesis export-efficiency wording is ambiguous. All candidate ratios are saved; the figure script uses export divided by contemporaneous surface light for the mgC:W panel.",
            "",
            "## Reproducibility level",
            "",
            "- R2: unavailable without original numerical outputs or digitized target values.",
            "- R1: assigned per automated ecological check above.",
            "- R0: failed checks and visual/numeric details that cannot be verified.",
            "",
        ]
    )
    (root / "reproduction_report.md").write_text("\n".join(lines), encoding="utf-8")


def reproduce_thesis(output_dir: str | Path, jobs: int = 1) -> Path:
    root = Path(output_dir).resolve()
    repository = Path(__file__).resolve().parents[2]
    manifests = repository / "configs" / "manifests"
    groups: dict[str, list[tuple[dict[str, np.ndarray], dict[str, Any]]]] = {}
    for group in ("reference", "constant_open", "constant_closed", "seasonal_200m", "critical_depth_100m"):
        run_ensemble(manifests / f"{group}.json", jobs, root / "results" / group)
        groups[group] = _runs(root / "results" / group)
    _write_summary_table(root, groups)
    checks = ecological_checks(groups)
    _write_checks(root, checks)
    figures = make_thesis_figures(root)
    write_figure_manifest(root / "figure_manifest.csv", figures)
    _write_report(root, checks, figures)
    manifest = {
        "schema_version": 1,
        "result_groups": {name: len(runs) for name, runs in groups.items()},
        "figures": [str(path) for path in figures],
        "ecological_checks": checks,
        "complete": True,
    }
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return root
