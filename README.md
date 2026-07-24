# SPLAS 1.1 - reproducible thesis pipeline

This repository reconstructs the 1-D, size-structured NPZD model used in Marco
Fianchini's 2021 thesis on light, plankton community structure, biodiversity,
stoichiometry, blooms, critical depth, and carbon export.

The numerical core remains Fortran. Python provides validated configuration,
serial/parallel ensemble orchestration, self-describing NetCDF output,
post-processing, figures, and reproducibility reports. The historical sources
are retained for archaeology; the active implementation is under `src/`.

## Quick start

Requirements: Python 3.10+, GNU Fortran, Make, NumPy, SciPy, and Matplotlib.

```bash
python3 -m pip install -e .
make build
make test
make smoke
```

Run one case:

```bash
splas run configs/reference/reference_open_constant_default.json \
  --output reproduction/results/reference_open
splas analyze reproduction/results/reference_open \
  --output reproduction/tables/reference_open
splas figures reproduction/results/reference_open \
  --output reproduction/figures/reference_open
```

Run or resume all thesis experiments:

```bash
make reproduce-thesis
# or:
splas reproduce-thesis --jobs 4 --output reproduction
```

The full suite contains 30 ten-year simulations and is intentionally much
larger than the two-day smoke run. Existing completed runs are reused only when
their resolved configuration hash matches.

## Scientific profiles

- `thesis_consistent` is recommended. It follows the thesis parameters and
  supply profile, uses a complete coupled Heun step, and uses a conservative
  dissolved-N/realized-production coupling required by the closed-system mass
  contract.
- `legacy_reconstructed` retains the historical uptake term and legacy timing
  choices without undefined memory access.

Every run writes `resolved_config.json`, `run_metadata.json`, `manifest.json`,
and `results.nc`. See [docs/user_guide.md](docs/user_guide.md) for the schema
and [docs/source_discrepancies.md](docs/source_discrepancies.md) for choices
that cannot be resolved silently.

## Attribution

SPLAS was developed by Marco De Pasquale and extended for light experiments by
Marco Fianchini. The historical README, code, and thesis remain the scientific
provenance for this reconstruction.
