# User guide

## Environment

Native build:

```bash
python3 -m pip install -e .
make build
make test
```

Debug build:

```bash
make build-debug
SPLAS_REQUIRE_CORE=1 python3 -m unittest discover -s tests -v
```

Container build:

```bash
docker build -t splas-1.1 .
docker run --rm -v "$PWD/reproduction:/opt/splas/reproduction" \
  splas-1.1 reproduce-thesis --output reproduction
```

## Commands

```bash
splas run CONFIG --output DIR
splas ensemble MANIFEST --jobs N --output DIR
splas analyze RESULTS_DIR --output DIR
splas figures RESULTS_DIR --output DIR
splas reproduce-thesis --jobs N --output reproduction
```

Equivalent targets:

```bash
make build
make test
make smoke
make reproduce-thesis
```

`make smoke` runs a two-day, seven-layer configuration with all 70 P and 70 Z
classes, then performs analysis and plotting. It is a functional check, not a
scientific substitute for the ten-year suite.

## Configuration

JSON supports recursive `extends` relative to the child file. Resolved configs
are validated and saved with a SHA-256 hash. Active dimensions are fixed at 70
P classes, 70 Z classes, and three output nutrient coordinates.

Important selectors:

- `scientific_profile`: `thesis_consistent` or `legacy_reconstructed`;
- `system.mode`: `open` or `closed`;
- `forcing.mode`: `constant`, `seasonal`, `daily`, `linear`, or `normal`;
- `numerics.seed`: required even for deterministic forcing;
- `numerics.output_interval_days`: physical days, not step counts.

The closed-mode validator rejects any non-zero supply, diffusion, P sinking, or
OM sinking.

## Canonical output

Each run directory contains:

- `results.nc`: NetCDF3 64-bit-offset dataset;
- `resolved_config.json`: complete inherited configuration;
- `run_metadata.json`: source/config hashes, versions, platform, compiler,
  seed, and execution time;
- `manifest.json`: completion flag and dataset identity.

Coordinates are `time`, `depth`, `phyto_size`, `zoo_size`, and `nutrient`.
Variables include P, Z, quota, cell number, N/diagnostic P/C pools, OM,
surface and attenuated light, three limitation diagnostics, growth/grazing,
net growth and total reduction relative to mu-max, biomass,
entropy/evenness/validity, cell carbon, carbon biomass, separate
exports, mass error, and non-negative correction.

NetCDF uses layer-center depth, positive downward. Nominal MLD/domain depths
are 200 m for M=41 and 100 m for M=21, while the last centers are 202.5 m and
102.5 m.

## Resume and parallelism

Ensembles run serially with `--jobs 1`; larger values use independent Python
processes. MPI is not required. A completed run is reused only if the saved
resolved-config hash equals the current hash.

## Interpretation safeguards

- `light_limitation` is the multiplicative fraction of potential growth
  retained under light; the limiting fraction is `1-light_limitation`.
- `total_growth_reduction` is computed from realized Droop-light production
  minus grazing relative to `mu_max*P`. It is not an unexplained sum of
  incompatible component metrics and may exceed one when grazing makes net
  growth negative.
- Shannon entropy and evenness are separate variables. Low-biomass cells are
  masked by `diversity_valid_P`.
- The 1% euphotic depth is an optical MLD* proxy, not the physical fixed MLD.
- Export efficiency is not a single hidden formula. Analysis saves export per
  biomass, per light, and per biomass-light.
