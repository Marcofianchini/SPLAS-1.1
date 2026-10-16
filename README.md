# SPLAS-CN 0.2 - dynamic carbon and nitrogen

A new conservative implementation alongside the original SPLAS-1.1 Fortran.
It retains the size-structured ecological questions of De Pasquale (2018) and
Fianchini (2021), while resolving carbon, nitrogen, cellular reserves, energetic
costs, real prey stoichiometry, detrital C/N and conservative vertical transport.

**Status: research prototype.** Numerical verification and culture-starvation
assessment are distinguished from field validation. Do not use this release to
claim validated ocean export, sequestration, or full C:N:P dynamics.

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

The improved idealized research configuration is in
`research/recommended_idealized.json`: buffered DIC with an external ledger,
16 particle bins,70+70 plankton classes,40 cells of5m and max_step0.0125day.
Its status must be read with the same-state step-refinement results, not
assumed to establish field validity. The historical irradiation is retained
for traceability and does not describe a measured ocean site.

The nominally open-N legacy counterpart has no matching C input and can
become strongly DIC limited. `run_carbon_boundary.py` tests that confound
with declared C reservoirs. More N does not mean a pure bottom-up N test
when another modeled resource changes its limiting role.

The distributed PDF embeds generated numerical evidence as an attached ZIP.
Extract it with `python scripts/extract_report_results.py report.pdf --output results`.
Checksums are verified automatically; `--output results` restores the directory expected
by the analysis scripts. Repository source and tracked files stay
in Git; the report attachment contains generated outputs.
