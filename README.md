# SPLAS-CN 0.2 - an ecological toy model with dynamic C and N

A new conservative implementation alongside the original SPLAS-1.1 Fortran.
It retains the size-structured ecological questions of De Pasquale (2018) and
Fianchini (2021), while resolving carbon, nitrogen, cellular reserves, energetic
costs, actual prey stoichiometry, detrital C/N and conservative vertical transport.

The scientific aim is **Pattern Oriented Modelling**: explain recognizable
ecosystem organization with relatively simple rules, following the progressive
mechanism experiments of the two theses. Success means jointly generating and
explaining size spectra, resource/grazer control, diversity, vertical organization
and temporal succession within declared toy regimes. Site forecasting and matching
an absolute ocean-export value are not the primary objective.

The [thesis-grounded POM protocol](research/pom/STUDY_PROTOCOL.md) reconstructs
that reasoning, separates imposed rules from emergent patterns, and orders the
next experiments. The [current POM diagnosis](research/pom/CURRENT_ASSESSMENT.md)
reuses existing results: competitive release is recognizable, while the joint
diversity/transfer pattern remains unresolved and diversity needs further
resolution checks. [The source argument review](research/pom/argument-review.md)
includes page-specific evidence and corrects unequal historical comparisons.

The pre-POM checkpoint is commit `ae55985`; its
[manifest](research/pom/BASELINE.json) records333 file hashes and117 community
cases. The methodological consolidation at `a8192dc` changed no model equations
and added no ecological simulations. Reproduce that checkpoint diagnosis with:

```bash
python scripts/analyze_pom.py --check-only
python scripts/analyze_pom.py
```

The subsequent [dynamic-stoichiometry study](research/stoichiometry/RESULTS.md)
executes **166 new 0D configurations**: resource-quality planes, matched closed
assembly, selective/uniform/no-grazer controls, physiology/storage alternatives,
step/size refinement and pulse sequences. It records full class histories and
process observations without changing the frozen numerical model. Lower consumer
N assimilation consistently contracts phytoplankton size and evenness in the
paired checks. Bulk C:N can conceal a much richer ingested diet; oscillatory
regimes remain sensitive in phase and fine composition. A robust stoichiometric
memory effect is **not established** by the pulse experiment.

```bash
python scripts/run_stoichiometry.py --workers 4
python scripts/run_stoichiometry_controls.py --workers 2
python scripts/run_stoichiometry_robustness.py --workers 3
python scripts/analyze_stoichiometry.py
```

All cases, including counterexamples, are in the
[study summary](research/stoichiometry/summary.csv). See the
[pattern assessment](research/stoichiometry/PATTERNS.json) and
[current primary-literature comparison](research/stoichiometry/LITERATURE.md).

The [independent 0D POM challenge](research/pom_validation_0d/RESULTS.md)
adds **221 runs** with internally registered pattern criteria and explicit
structural alternatives, followed by declared discriminating controls. The whole
reference plane at a finer step retains enrichment/size expansion and size
groups. The quantity/quality pattern remains conditional and fails its declared
generality criterion; the opposite diversity/Z:P signature is not reproduced.
Rare-trait redistribution can reorganize apparently stationary low-N communities
at unchanged C/N inventories. Thus finite-window size groups do not establish
invasion-resistant coexistence. This is **partial POM support**, not general
ecological validation.

Common-state pulse controls now support a phase-dependent interaction of event
order and a physiological component of memory in reserve common gardens.
Lower N retention still decreases grazer production and producer mean size in
the new recycling contrasts, but the earlier evenness direction is **not
universal**. See the [protocols and evidence](research/pom_validation_0d/PROTOCOL.md),
[complete summary](research/pom_validation_0d/summary.csv),
[fine-plane decisions](research/pom_validation_0d/summary.json) and
[verification](research/pom_validation_0d/VERIFICATION.json). Equations,
integrator and previous data remain unchanged.

```bash
python scripts/run_pom_validation.py --workers 4
python scripts/run_pom_validation.py --stage memory --workers 4
python scripts/run_pom_validation.py --stage garden --workers 4
python scripts/run_pom_followup.py
python scripts/run_pom_refinement.py
python scripts/run_pom_rare_pool.py
python scripts/analyze_pom_invasion.py
python scripts/analyze_pom_validation.py
python scripts/plot_pom_validation.py
```

**Status: research prototype.** The following implementation and auxiliary
physiology checks support the POM study; they are not its ecological success
score. Units, element closure and robustness of the claimed patterns remain
necessary in an idealized model.

| Evidence | Assessment |
|---|---|
|26 scientific implementation/observation tests | Pass locally; the workflow repeats the same suite. The latest three checks reconstruct the five-pool N network from the independent RHS, match all new initial inventories and conserve physical-bin mass in both trophic size domains. |
|Carbon/nitrogen inventories and positivity | Verified for the registered experiments; inspect exact residuals in the validation summary. |
|Annual time-step and16/32 particle-bin comparisons | Pass the declared numerical thresholds from a common state. |
|Vertical resolution |40/80 fails. Between80/160, export changes2.18%, but N stock6.62% and NPP5.44% exceed the5% bulk criterion. Overall convergence is not demonstrated. |
|Four N-starvation cultures | Partial exploratory skill; initial-value persistence outperforms the shared model for Micromonas. |

Implementation evidence is in [VALIDATION.md](research/VALIDATION.md). Ecological
pattern assessment and the new study design are in [research/pom](research/pom).
Closed element ledgers do not alone identify the ecological mechanisms.

## Run

Python 3.10+; the recorded reference environment is Python3.12, NumPy2.3.5,
SciPy1.17, Numba0.68. Numba is strongly recommended for ensembles.
`research/requirements-reference.txt` records exact scientific-package versions
as install constraints; Python3.12.14/Linux was used for the reported runs.

```bash
python -m pip install -c research/requirements-reference.txt -e '.[fast,research]'
python -m unittest discover -s tests -v
python -m splas2 --days 30 --output results/example.npz
# Revised idealized column (a short run, not a spin-up):
python -m splas2 --config research/recommended_idealized.json --days 30 --output results/column.npz
```

The default CLI uses20 phytoplankton and20 grazer classes. The thesis
registry uses70+70 classes unless explicitly testing another resolution.
Both C and N are in mmol element/m3 and all C:N ratios are **molar**.

```bash
python scripts/prepare_cultures.py
python scripts/calibrate_cultures.py
python scripts/run_experiments.py --workers 4
python scripts/run_mechanisms.py
python scripts/run_physiological_protocols.py
python scripts/run_supplement.py
python scripts/run_carbon_boundary.py
python scripts/run_recommended.py
python scripts/run_recommended_refinement.py
python scripts/run_export_interventions.py
python scripts/analyze_results.py
python scripts/assess_validation.py
python scripts/make_validation_report.py
python scripts/package_report_data.py
```

`run_experiments.py --filter F02_,F03_ --workers 2` selects a family. Every case
saves actual configuration, initial inventory, diagnostics and reduced
trajectories. Cached runs are keyed by configuration and model source hash.
Large generated trajectories are excluded from Git. The scientific report
and aggregate validation results document the tested revision and failures.

## What changes scientifically

- Independent C-only and organic-N pools, with carbon explicitly accompanying
  organic N. Carbon fixation, energetic assimilation, division and maintenance
  use the same conservative element ledgers. Consumer homeostasis remains an
  assumption; phytoplankton and detrital C:N vary.
- Grazing is normalized by `phi.T @ prey_C`, consumes actual prey C:N, obeys
  its maximum ingestion and releases excess assimilated elements.
- Finite-volume face fluxes conserve both elements at variable diffusivity;
  diffusivity in m2/s is converted to the day timebase.
- Explicit SSPRK3 with embedded error control is the reference integrator.
  Long ensembles use an IMEX option: positive SSPRK3 with a transport CFL bound for nonstiff processes
  plus an element-conservative implicit N-assimilation solve. The combined
  method is **first order**; step refinement must be checked. No state clipping.
- Stock Z:P, flux TTE, gravitational export and higher-trophic export are
  different reported quantities. Euphotic depth is not a proxy for mixing depth.

See [mathematical specification](research/MODEL.md),
[data provenance](data/README.md), [experiment registry](research/experiment_registry.json)
[thesis experiment map](research/THESIS_REPLAY.md),
and [validation conclusions](research/VALIDATION.md).

## Relation to the old code and Moloney-Field

The original Fortran and Python controller are retained unchanged. Their
instructions are in [legacy README](research/README_legacy.md). No compiled
legacy execution is claimed. The attached MoloneyField0D v0.1 provided a
conservative developmental benchmark, not an archival Moloney-Field release.
Moloney's1988 thesis used fixed group C:N mass ratios; this implementation is
a new dynamic-stoichiometry hypothesis. The original1991 article was not
fully available for source-level equivalence verification.

Montagnes etal1994 report `0.109*V^0.991 pgC/cell`, whereas the historical
SPLAS carbon conversion uses0.0109. This correction changes reference
biomass by a factor10 and must not be hidden inside a fitted parameter.

## Literature and limits

Culture observations: Liefer etal2019, DOI10.3389/fmicb.2019.00763;
raw dataset Liefer etal2026, DOI10.1594/PANGAEA.900248 (CC-BY4.0).
The archive explicitly lacks final formal author approval; published quota
means were cross-checked before use. Ot/Tw estimate shared parameters;
Tp/Micromonas observations are excluded from the optimizer. This is an
exploratory withheld-species comparison, conditional on observed initial
quotas/densities, and not an untouched prospective test set.

[Current literature matrix](research/LITERATURE.md) explains model decisions
against2024-2026 primary studies. P, Fe, Si, O2, carbonate chemistry, bacteria,
mixotrophy, aggregation and explicit macromolecular allocation are absent.
More physiological states do not themselves establish greater predictive skill.

The revised idealized control configuration is in
`research/recommended_idealized.json`: buffered DIC with an external ledger,
16 particle bins,70+70 plankton classes,40 cells of5m and max_step0.0125day.
The40-cell vertical gate failed; this configuration is a traceable control,
not a spatially converged export estimate. The nested80/160-cell result is
reported separately. Temporal/particle refinement does not repair a failed
vertical-resolution check. The historical irradiation is retained
for traceability and does not describe a measured ocean site.

Element conservation does not imply stationarity: an openN input can exceed
exports while phytoplankton biomass appears stable. The result summary reports
net N retention and DIN accumulation separately from biomass drift.

The nominally open-N legacy counterpart has no matching C input and can
become strongly DIC limited. `run_carbon_boundary.py` tests that confound
with declared C reservoirs. More N does not mean a pure bottom-up N test
when another modeled resource changes its limiting role.

The distributed PDF embeds generated numerical evidence as an attached ZIP.
Extract it with `python scripts/extract_report_results.py report.pdf --output results`.
Checksums are verified automatically; `--output results` restores the directory expected
by the analysis scripts. Repository source and tracked files stay
in Git; the report attachment contains generated outputs.
