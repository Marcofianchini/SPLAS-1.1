# SPLAS 1.1 thesis reproduction report

The reproduction pipeline and experiment manifests are present, but the current
execution image has no Fortran compiler. Therefore no model simulation was
executed here and all thesis results remain R0 (not yet verified) in this
checkout.

Compiler-independent status:

- 34 versioned configurations validate;
- 18 scientific/configuration/NetCDF/plotting tests pass;
- 4 core tests are correctly skipped because the extension is absent;
- f2py parses the new Fortran public interface;
- the re-uploaded 51-page thesis was inspected and mapped to the core.

After `make build && make test && make reproduce-thesis` in a GNU Fortran
environment, this file is replaced automatically by a result-specific R1/R0
report. R2 remains unavailable without the original numerical outputs or
digitized target tolerances.
