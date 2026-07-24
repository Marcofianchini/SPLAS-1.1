# Baseline audit

## Sources and integrity

The immutable supplied archive was `SPLAS-1.1-main.zip`, SHA-256
`bf8ed3f6584b365a3af470defa8164f568845f046dfa96ba5db25e15d5ef6449`.
The usable re-uploaded thesis was
`Fianchini_tesilaurea_2307_corretta(2).pdf`, SHA-256
`19b0ec52dd06d35422e168dc7af506e85fb42d6c130717064bee2cf7faa2922f`
(51 pages). The first local thesis copy was empty and was not used.

Original source hashes:

| Original file | SHA-256 |
|---|---|
| `Makefile` | `4f93348aa5486f880034c3c3442ac175f9bd7eede228983d0e040c1a67961d7d` |
| `README.md` | `18555f4750458b7ad2eb048864d5a59f9801822600eba327e29c3bc01147e846` |
| `functions.f90` | `50668ad994360326c5d3d00dbe00009f9d8b4300d1abae1423fde2ed707b8196` |
| `heuns_method.f90` | `be7d16142c48ec19e2b7f0c4c9be4abd73346cd9c31c89549f050f5ec283215b` |
| `solver.f90` | `a1886cff0154ad9e4616519c02c0dedec77cf268c6905804c4ac031d4add45a7` |
| `solver_closed_system.f90` | `10cfffe37b0e69338b87374fd434c1961d68e4ad68fc2bb340b7ee127881d82f` |
| `splas.f90` | `52f1dffd496f945f872783edfe61e8356952641fe8af1dd2154d78fcbf7a0e0f` |
| `splas_py_wrapper.f90` | `4220526f3e9400249f1a9a7a699b48d2af0b13bb08159008eeae486e931fbe9b` |
| `splas.py` | `25890f067c9db970bd4c099fa2959e27257a754bb3ea9c6c71c08250091b35a5` |
| `splasIO.py` | `c38ffe4b462f600c4877e28902e7a8fcab9af5b08b7c832f670dd2b167849f3c` |
| `splasplot.py` | `c39c599f21ea456b430ae72e604fd87bb6101469b2ad116683b2b076eae4cecb` |
| `main.py` | `68e371adb09614dc745bbdb43777682ba3492a7812e2832d418145df15716229` |
| `input_default.json` | `f13d2c7021c1966821881934f06a33e453aeb0555bf043c042e76dd466a6b6f0` |

## Baseline build and import

- Original `make all` failed immediately because `gfortran` is absent in the
  present execution image. Independently, the original Makefile hard-codes
  `f2py3.5`, which is not available under Python 3.12.
- No compiled `splaspy` extension is included in the ZIP.
- Importing the historical controller therefore fails at
  `from splaspy import wrapper_py`.
- `main.py` imports `mpi4py` unconditionally, so a serial run also requires MPI.
- The supplied environment has Python 3.12.13, NumPy 2.3.5, SciPy 1.17.0, and
  Matplotlib 3.10.8, but no GNU Fortran compiler, Docker, pytest, xarray, or
  mpi4py.

## Missing inputs and implicit workflow

The archive references but does not contain:

- `input/`, `inputIC/`, `input_closed_world_tv/`, `input_test_1.5/`, and
  `sensitivity_json/`;
- dozens of named JSON configurations such as `input/input3.json`;
- `SPLAS_PICKLE/` and output directories;
- original numerical results and the scripts/configurations that produced
  every thesis panel.

Light mode is selected by commenting one of two lines in `solver.f90`.
Open/closed mode is selected by compiling a different file named `solver`,
not by `splas_simulation` in JSON. The JSON value is used only as a run name.

## Hard-coded scientific parameters

| Parameter/choice | Historical location | Value/behavior |
|---|---|---|
| P/Z/nutrient dimensions | `functions.f90:9-10`, `heuns_method.f90:8-9` | 70/70/3 |
| Detritus remineralization | `solver.f90:115` | 0.15 d-1 |
| Detritus sinking | `solver.f90:116` | 8 m d-1 |
| Zoo mortality recycling | `solver.f90:117` | 0.5 |
| Closed alternatives | `solver_closed_system.f90:115-117` | 0.40, 0, 0.8 |
| Water/self-shading | `functions.f90:611-612` | 0.03, 0.01 |
| Chl:C and PI slope | `functions.f90:647-651` | 0.02, 0.32 |
| Seasonal vs constant light | `solver.f90:165-166` | source edit |
| Interior supply | `solver.f90:209,237` | `k * 0.1 * S`, contrary to thesis |
| Initial quota | `functions.f90:420` | JSON value multiplied by 1e-11 |

## Numerical defects in the baseline

1. `det_old` and `det_pred` are read before complete initialization in the open
   solver.
2. The closed solver calls `dcopy(dim, det_temp(i), 1, det_old(i), 1)` on
   scalars, causing out-of-bounds memory access.
3. The predictor is built one vertical layer at a time, yet its derivative
   reads `i-1` and `i+1`; neighboring predictor values are therefore stale or
   undefined.
4. Boundary values are applied only after the interior corrector, so they are
   absent from the predicted coupled state.
5. Only N is read in uptake, but a scalar right-hand side is assigned to all
   three nutrient pools.
6. Non-negativity is applied only to nutrients. NaN/negative checks are
   partial and output-time dependent.
7. Printed states are averages over step counts, while the `days` coordinate
   can have a different length and is interpreted ambiguously.
8. The historical "Shannon index" routine actually returns normalized Shannon
   evenness and divides by total biomass without a low-biomass mask.
9. Pickles have no stable schema, units, complete configuration, seed, source
   hash, or compiler metadata.

## Baseline conclusion

The archive expresses the intended ecological model but cannot reproduce the
thesis as supplied. A parametrized core, complete coupled Heun step, generated
configurations, schema-bearing outputs, and independent analysis/figure
pipeline are necessary rather than optional cleanup.
