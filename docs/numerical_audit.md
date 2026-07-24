# Numerical audit

## Integrator reconstruction

The historical vertical Heun loop calculated one predicted layer and
immediately read its not-yet-calculated neighbors. The new core uses this
ordering for every time step:

1. copy the complete old state;
2. compute old-state light, biological rates, and vertical derivatives for all layers;
3. build the complete predicted state;
4. validate/floor only roundoff-scale negative predictor values;
5. recompute light, rates, and derivatives on the complete predicted state;
6. apply the trapezoidal corrector once to every state variable;
7. validate and save only on physical output times.

This ordering is visible in `run_model` lines 145-178. There is no source-file
swap and no use of partially populated neighboring arrays.

## Precision, bounds, and invalid values

- All public numerical arrays and internal state use double precision/real64.
- `make build-debug` enables bounds checks, floating-point traps, backtraces,
  and warnings with GNU Fortran.
- Non-finite state returns status 30.
- A negative value below -1e-9 returns status 31 rather than being silently
  clipped.
- Roundoff-scale negatives are set to zero and their cumulative nitrogen-mass
  correction is saved as `nonnegative_correction`.
- Quota is bounded at Qmin, consistent with the Droop definition; the choice
  is explicit and does not alter P/N inventories directly.

## Transport and boundaries

Diffusion uses centered second differences and zero-gradient ghost values at
surface and bottom. P and OM sinking use a conservative downward upwind flux;
the domain-integrated flux telescopes to the bottom export. Closed mode sets
diffusion, P sinking, OM sinking, and supply to zero through validated config.

Diffusivity is configured in m2 s-1 and converted once to m2 d-1 inside the
daily-rate right-hand side.

## Mass and export budgets

For `thesis_consistent`, dissolved-N loss is coupled to realized P production.
With `gamma=1`, no transport, and no supply, the inventory
\[
\Delta z\sum_k(N_k+\sum_iP_{ik}+\sum_jZ_{jk}+OM_k)
\]
must remain within 1e-6 relative tolerance or the core returns status 20.

Open runs save phytoplankton and detrital lower-boundary export separately.
The P flux uses Montagnes cell carbon; OM uses Redfield C:N, matching the
thesis's stated approximation.

## Determinism

Constant, seasonal, daily, and linear forcing contain no random state. Normal
forcing uses an explicit Park-Miller generator and Box-Muller transform seeded
from config. The complete seed and configuration hash are stored per run.

## Verification status in the supplied execution image

The f2py parser successfully generated the expected public signature from the
new Fortran source. The present image has no Fortran compiler, so compilation,
debug bounds execution, smoke integration, mass-budget execution, time-step
convergence, and output-frequency invariance could not be executed here.

The compiler-independent suite currently passes 18 scientific/config/I/O/plotting tests;
four core tests are skipped only when the extension is absent. `make test`
first builds the extension and sets `SPLAS_REQUIRE_CORE=1`, converting those
tests into mandatory checks in a compiler-equipped environment.

No R2 or R1 model result is claimed until those tests and at least the smoke
run complete.
