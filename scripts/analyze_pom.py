"""Reinterpret frozen results for a thesis-grounded POM study; never run the ODEs.

python scripts/analyze_pom.py --freeze-baseline  # once, before further development
python scripts/analyze_pom.py                   # verify, diagnose, draw figures
python scripts/analyze_pom.py --check-only      # read-only checkpoint verification
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import subprocess
import zipfile

import numpy as np
import scipy
from scipy.signal import find_peaks
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "research/pom"
FIG = ROOT / "results/pom"
BASE = "ae55985e0f05ae05f0930c1bdfc1cca6f04052b6"
TAG = "splas-cn-0.2-pom-baseline-20261016"
REGISTRIES = ["experiment_registry.json", "supplement_registry.json",
              "carbon_boundary_registry.json", "recommended_registry.json",
              "recommended_refinement_registry.json", "export_intervention_registry.json"]


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def load(path):
    return json.loads(path.read_text())


def save(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def freeze():
    target = DOC / "BASELINE.json"
    if target.exists():
        raise RuntimeError("Baseline already exists; do not silently replace a checkpoint")
    current = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    if current != BASE:
        raise RuntimeError("Freeze must start at the recorded pre-POM commit")
    files = sorted(p for p in (ROOT / "results").rglob("*")
                   if p.is_file() and p.suffix in {".json", ".npz", ".png", ".pdf"}
                   and "pom" not in p.relative_to(ROOT / "results").parts)
    files += [ROOT / "research" / n for n in REGISTRIES + [
        "validation_summary.json", "acceptance_assessment.json", "experiment_summary.csv",
        "calibrated_physiology.json", "requirements-reference.txt", "SOURCE_PROVENANCE.json"]]
    files += sorted((ROOT / "splas2").glob("*.py"))
    files += sorted((ROOT / "tests").glob("*.py"))
    report = ROOT / "output/pdf/SPLAS_Sviluppo_Verifica_2026-10-16.pdf"
    if report.exists():
        files.append(report)
    DOC.mkdir(parents=True, exist_ok=True)
    records = [{"path": str(p.relative_to(ROOT)), "bytes": p.stat().st_size,
                "sha256": sha(p)} for p in sorted(set(files))]
    save(target, {
        "schema_version": 1, "date": "2026-10-16", "purpose": "Frozen pre-POM evidence, not a new model release",
        "baseline_commit": BASE, "baseline_tree": "730e010891828e29f404f260aa9d2eb5d1457633",
        "local_tag": TAG, "model_sha256": sha(ROOT / "splas2/model.py"),
        "community_scenarios": 117, "new_ecological_simulations_in_this_revision": 0,
        "source_provenance": load(ROOT / "research/SOURCE_PROVENANCE.json"),
        "files": records,
        "availability": "Large generated results are excluded from Git; the preceding distributed PDF embeds its numerical archive. This manifest records the available checkpoint bytes."
    })


def verify():
    baseline = load(DOC / "BASELINE.json")
    missing, changed = [], []
    for rec in baseline["files"]:
        path = ROOT / rec["path"]
        if not path.exists():
            missing.append(rec["path"])
        elif sha(path) != rec["sha256"]:
            changed.append(rec["path"])
    if missing or changed:
        raise RuntimeError(json.dumps({"missing": missing, "changed": changed}))
    source = (ROOT / "splas2/model.py").read_text()
    registry = sum([load(ROOT / "research" / n) for n in REGISTRIES], [])
    assert len(registry) == len({c["id"] for c in registry}) == 117
    cases = {}
    for c in registry:
        path = ROOT / "results/thesis_replay" / (c["id"] + ".json")
        record = load(path)
        sig = hashlib.sha256((json.dumps(c, sort_keys=True) + source).encode()).hexdigest()
        if sig != record["run_signature"]:
            raise RuntimeError("Stale model/case signature: " + c["id"])
        archive = path.with_suffix(".npz")
        with zipfile.ZipFile(archive) as z:
            if z.testzip():
                raise RuntimeError("Corrupt archive: " + c["id"])
        with np.load(archive) as a:
            s = record["summary"]
            assert np.isclose(a["spectrum_N"].sum(), s["phyto_N_mmol_m2"], rtol=1e-11)
            assert np.isclose(a["spectrum_C"].sum(), s["phyto_C_mmol_m2"], rtol=1e-11)
            assert np.isclose(a["spectrum_Z_C"].sum(), s["zoo_C_mmol_m2"], rtol=1e-11)
        cases[c["id"]] = record
    return baseline, cases


def association(x, y):
    x, y = np.asarray(x), np.asarray(y)
    return {"Pearson_raw": float(np.corrcoef(x, y)[0, 1]),
            "Spearman_raw": float(spearmanr(x, y).statistic)}


def rebin_log_mass(sizes, mass, n):
    """Conservative overlap integral, uniform within each native log-ESD interval.

    This explicit measurement convention avoids counting different native class
    numbers as different biological diversity. It is not a new ecological process.
    """
    x = np.log10(sizes)
    old = np.r_[x[0], .5 * (x[1:] + x[:-1]), x[-1]]
    new = np.linspace(x[0], x[-1], n + 1)
    overlap = np.maximum(0, np.minimum(new[1:, None], old[None, 1:])
                         - np.maximum(new[:-1, None], old[None, :-1]))
    rebinned = overlap @ (mass / np.diff(old))
    assert np.isclose(rebinned.sum(), mass.sum(), rtol=1e-12)
    return rebinned


def entropy_evenness(mass):
    p = mass[mass > 0] / mass.sum()
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)) / np.log(len(mass)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze-baseline", action="store_true")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    if args.freeze_baseline:
        freeze()
    baseline, cases = verify()
    if args.check_only:
        print(json.dumps({"frozen_files_verified": len(baseline["files"]),
                          "community_signatures_verified": len(cases),
                          "model_unchanged": True, "new_simulations": 0}))
        return
    rows = {k: v["summary"] for k, v in cases.items()}
    cohort = [k for k in cases if k == "D00_reference" or
              k.startswith(("D01_", "D02_", "D04_", "D05_", "D06_"))]
    supply = sorted([k for k in cohort if k == "D00_reference" or k.startswith("D01_")],
                    key=lambda k: cases[k]["config"]["nutrient_supply"])
    mortality = sorted([k for k in cohort if k == "D00_reference" or k.startswith("D02_")],
                       key=lambda k: cases[k]["physiology"]["zoo_mortality"])
    def stock_N(k):
        s = rows[k]
        ans = s["zoo_C_mmol_m2"] / cases[k]["physiology"]["zoo_CN"] / s["phyto_N_mmol_m2"]
        # Independent identity connecting both measurement currencies.
        assert np.isclose(ans, s["Z_to_P_stock_C"] * s["phyto_CN_inventory_ratio"] /
                          cases[k]["physiology"]["zoo_CN"], rtol=1e-12)
        return ans
    out = {
        "schema_version": 1, "date": "2026-10-16", "baseline_commit": BASE,
        "baseline_manifest_sha256": sha(DOC / "BASELINE.json"),
        "analysis_script_sha256": sha(Path(__file__)), "new_simulations": 0,
        "analysis_environment": {"Python": platform.python_version(), "NumPy": np.__version__, "SciPy": scipy.__version__},
        "analysis_role": "Retrospective exploratory POM diagnosis; not independent validation or preregistration",
        "cohort_rule": "D00,D01,D02,D04,D05,D06: steady default-closure perturbations corresponding to De Pasquale3.4.4. No claim of an identical historical ensemble.",
        "cohort_ids": cohort, "associations": {}, "supply_hump_diagnostics": {},
        "grid_observables": [], "vertical_profiles": [],
        "aggregation_warning": "Correlations across deterministic scenarios depend on the perturbation distribution. No p-values, causal identification, or population inference."
    }
    for name, ids in [("pooled_24", cohort), ("N_supply_9", supply), ("Z_mortality_5", mortality)]:
        ep = [rows[k]["phyto_evenness_N"] for k in ids]
        ez = [rows[k]["zoo_evenness_C"] for k in ids]
        block = {"case_ids": ids, "metrics": {}}
        for metric in ["stock_ZN_PN", "Z_to_P_stock_C", "TTE_flux_C"]:
            vals = [stock_N(k) if metric == "stock_ZN_PN" else rows[k][metric] for k in ids]
            block["metrics"][metric] = {"prey_evenness_N": association(ep, vals),
                                       "predator_evenness": association(ez, vals)}
        out["associations"][name] = block
    for name, mass, ev in [("phyto", "phyto_N_mmol_m2", "phyto_evenness_N"),
                           ("zoo", "zoo_C_mmol_m2", "zoo_evenness_C")]:
        x = np.log10([rows[k][mass] for k in supply])
        y = np.array([rows[k][ev] for k in supply])
        coef = np.polyfit(x, y, 2)
        fit = np.polyval(coef, x)
        out["supply_hump_diagnostics"][name] = {
            "case_ids": supply, "log10_mass_field": mass, "evenness_field": ev,
            "quadratic_coefficients": coef.tolist(),
            "R2_descriptive": float(1 - np.sum((y-fit)**2) / np.sum((y-y.mean())**2)),
            "interior_observed_maximum": bool(np.argmax(y) not in (np.argmin(x), np.argmax(x))),
            "quadratic_vertex_log10_mass": float(-coef[1]/(2*coef[0])),
            "observed_log10_mass_range": [float(x.min()), float(x.max())],
            "interpretation": "Interior diversity maximum is compatible with a hump. Rebounds remain; the quadratic is a descriptive diagnostic, not a selection criterion, causal proof or new calibration."
        }
    for k in ["D07_size_bins_35", "D00_reference", "D07_size_bins_140"]:
        with np.load(ROOT / "results/thesis_replay" / (k + ".npz")) as a:
            rec = {"id": k, "n_classes": len(a["sizes"]),
                   "dominant_ESD_um": float(a["sizes"][a["spectrum_N"].argmax()]),
                   "native_evenness_N": rows[k]["phyto_evenness_N"],
                   "common_log_bin_evenness_N": {
                       str(n): entropy_evenness(rebin_log_mass(a["sizes"], a["spectrum_N"], n))
                       for n in [7, 14, 28]}}
            out["grid_observables"].append(rec)
    for k in ["F00_open_reference", "F02_semiannual", "F05_no_shading", "F05_no_grazers",
              "F04_mixing_0.1", "F04_mixing_10", "R00_buffered_reference"]:
        with np.load(ROOT / "results/thesis_replay" / (k + ".npz")) as a:
            p, z = a["profile_N"], a["depth"]
            j = int(p.argmax())
            out["vertical_profiles"].append({"id": k, "maximum_depth_m": float(z[j]),
                "maximum_type": "surface_cell" if j == 0 else ("bottom_cell" if j == len(z)-1 else "interior"),
                "interior_peak_depths_by_relative_prominence": {
                    str(q): [float(z[i]) for i in find_peaks(p, prominence=q*p.max())[0]]
                    for q in [.02, .05, .1]},
                "note": "Local maxima of the saved mean N-biomass profile. Chlorophyll maxima, instantaneous multiple maxima and optical causation are not inferred."})
    mechanisms = load(ROOT / "results/mechanisms/mechanisms.json")
    out["fixed_CN_assembly"] = mechanisms["fixed_CN_assembly"]
    out["dynamic_assembly"] = []
    for k in cases:
        if not k.startswith("D11_"):
            continue
        with np.load(ROOT / "results/thesis_replay" / (k + ".npz")) as a:
            out["dynamic_assembly"].append({"id": k, "n_zoo": cases[k]["config"]["n_zoo"],
                "dominant_ESD_um": float(a["sizes"][a["spectrum_N"].argmax()]),
                "geomean_ESD_um": rows[k]["N_weighted_geomean_ESD_um"],
                "ESD95_um": rows[k]["N_weighted_ESD95_um"]})
    out["metric_crosswalk"] = {
        "historical_ZP": "N-stock Z:N / P:N; compute Z_C / zoo_CN / P_N. C-stock ratio and C-production ratio are additional different observables.",
        "Fianchini_export_efficiency": "Export weighted by biomass then normalized by irradiation (pp39-40,43-46); exact weighting requires historical code/run provenance. Not export/NPP.",
        "time_windows": "Theses often average10years including transients; current spectra/profiles average last364days of3652day runs. Full ten-year size-resolved trajectories cannot be reconstructed from the reduced archives.",
        "seasonality": "Fianchini p33 deliberately specifies two light maxima per year; emergent targets are response phase, amplitude, succession and persistence, not the imposed count of light peaks."
    }
    save(DOC / "current_assessment.json", out)
    draw_figures(cases, cohort, supply, stock_N, mechanisms)
    print(json.dumps({"frozen_files_verified": len(baseline["files"]), "scenarios": len(cases),
                      "POM_cohort": len(cohort), "new_simulations": 0, "outputs": str(DOC)}))


def draw_figures(cases, cohort, supply, stock_N, mechanisms):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    FIG.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.spines.top": False, "axes.spines.right": False})
    def figsave(fig, name):
        fig.savefig(FIG / (name + ".png"), dpi=180, bbox_inches="tight", facecolor="white")
        fig.savefig(FIG / (name + ".pdf"), bbox_inches="tight", facecolor="white")
        plt.close(fig)
    fig = plt.figure(figsize=(9.2, 6.8), constrained_layout=True)
    gs = fig.add_gridspec(2, 2)
    ax = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])]
    fixed = mechanisms["fixed_CN_assembly"]
    ax[0].plot([r["n_zoo"] for r in fixed], [r["winning_phyto_ESD"] for r in fixed], "o-", label="Monod, C:N fisso")
    dyn = [v for k, v in cases.items() if k.startswith("D11_")]
    modes = []
    for c in dyn:
        with np.load(ROOT / "results/thesis_replay" / (c["id"]+".npz")) as a:
            modes.append(float(a["sizes"][a["spectrum_N"].argmax()]))
    ax[0].plot([c["config"]["n_zoo"] for c in dyn], modes, "s--", label="Riserve, C:N dinamico")
    ax[0].set(xlabel="Classi di predatori introdotte", ylabel="ESD dominante P (µm)", yscale="log", title="A  Assemblaggio progressivo")
    ax[0].legend(frameon=False, fontsize=9)
    bath = sorted([v for k, v in cases.items() if k.startswith("D14_Cbath_supply_")], key=lambda c: c["config"]["nutrient_supply"])
    for ids, label in [(supply, "Scambio DIC 0,01/d")]:
        ax[1].plot([cases[k]["config"]["nutrient_supply"] for k in ids], [cases[k]["summary"]["N_weighted_ESD95_um"] for k in ids], "o-", label=label)
    ax[1].plot([c["config"]["nutrient_supply"] for c in bath], [c["summary"]["N_weighted_ESD95_um"] for c in bath], "s--", label="Serbatoio DIC 1/d")
    ax[1].set(xscale="log", yscale="log", xlabel="Apporto N (mmol m⁻³ d⁻¹)", ylabel="ESD95 ponderato N (µm)", title="B  Condizioni di risorsa")
    ax[1].legend(frameon=False, fontsize=9)
    refs = ["D07_size_bins_35", "D00_reference", "D07_size_bins_140"]
    for k in refs:
        with np.load(ROOT / "results/thesis_replay" / (k+".npz")) as a:
            b = rebin_log_mass(a["sizes"], a["spectrum_N"], 14)
            edges = np.linspace(np.log10(a["sizes"][0]), np.log10(a["sizes"][-1]), 15)
            ax[2].stairs(b/b.sum(), 10**edges, label=str(len(a["sizes"]))+" classi")
    ax[2].set(xscale="log", xlabel="ESD (µm), 14 intervalli comuni", ylabel="Frazione biomassa N", title="C  Dominanza e distribuzione")
    ax[2].legend(frameon=False, fontsize=9)
    figsave(fig, "01_assembly_resource_resolution")
    fig, ax = plt.subplots(2, 2, figsize=(9.2, 7.2), constrained_layout=True)
    for axis, mass, ev, title in [(ax[0,0], "phyto_N_mmol_m2", "phyto_evenness_N", "A  Fitoplancton: diversità-biomassa"),
                                  (ax[0,1], "zoo_C_mmol_m2", "zoo_evenness_C", "B  Zooplancton: diversità-biomassa")]:
        factor = 5.5 if mass.startswith("zoo") else 1.
        axis.plot([cases[k]["summary"][mass]/factor for k in supply], [cases[k]["summary"][ev] for k in supply], "o-", color="#147d92")
        axis.set(xscale="log", xlabel="Biomassa N (mmol m⁻²)", ylabel="Shannon evenness", title=title)
    families = [("D00", "Riferimento", "#35495c"), ("D01", "Apporto N", "#147d92"), ("D02", "Mortalità Z", "#d07830"),
                ("D04", "Selettività", "#7b5ba5"), ("D05", "Mortalità P", "#359968"), ("D06", "Perdite sinking P", "#b65362")]
    for prefix, label, color in families:
        ids = [k for k in cohort if k.startswith(prefix)]
        for axis, ev in [(ax[1,0], "phyto_evenness_N"), (ax[1,1], "zoo_evenness_C")]:
            axis.scatter([cases[k]["summary"][ev] for k in ids], [stock_N(k) for k in ids], s=30, color=color, label=label)
    for axis, title in [(ax[1,0], "C  Rapporto Z:P e diversità delle prede"), (ax[1,1], "D  Rapporto Z:P e diversità dei predatori")]:
        axis.set(xlabel="Shannon evenness", ylabel="Rapporto stock Z_N / P_N", title=title)
    ax[1,1].legend(frameon=False, fontsize=9, ncol=2)
    figsave(fig, "02_diversity_trophic_patterns")


if __name__ == "__main__":
    main()
