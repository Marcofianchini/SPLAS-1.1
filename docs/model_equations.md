# Model equations and implementation map

## State and grid

SPLAS resolves 70 logarithmic phytoplankton ESD classes (1-200 um), 70
zooplankton ESD classes (2-1000 um), one prognostic dissolved-N pool, and one
detrital-N pool in M vertical layers. The fixed-size contract is validated in
`src/splas/config.py`; it is not advertised as dynamic.

Depth coordinates follow thesis Eq.1,

\[
z_k=(k-\tfrac12)\Delta z,
\]

while nominal domain/MLD depth follows Eq.2,
\[
h=(M-1)\Delta z.
\]

The source profile implements Eq.3: zero at the surface, \(0.1S\) in each
interior layer, and \(S\) at the bottom.

## Equations 4-7: NPZD dynamics

The coupled right-hand side is `src/fortran/splas_core.f90:188-331`.

| Thesis equation | Process | Active implementation |
|---|---|---|
| Eq.4 | P growth, grazing, mortality, sinking, diffusion | realized production at line 237; grazing at 244-247; P transport at 299-316 |
| Eq.5 | Z assimilation, density-dependent mortality, diffusion | lines 253-255 and 318-328 |
| Eq.6 | N supply, uptake/production loss, unassimilated grazing, remineralization, diffusion | profile switch at 257-263; supply at 267-274; diffusion at 291 |
| Eq.7 | OM remineralization, P mortality, fraction γ of Z mortality, egestion, sinking, diffusion | lines 264-265 and 292-296 |

`thesis_consistent` uses realized P production as the dissolved-N loss so that
the stated closed-system mass constraint is enforceable. `legacy_reconstructed`
uses the literal historical cellular uptake term. Both choices are recorded in
the resolved configuration and output metadata.

## Equations 8-12: light, quota, cells, and carbon

Eq.8 is the product of the Droop factor and the saturating PI factor:

\[
\mu_i(z)=\mu_i^0(1-Q_i^{min}/Q_i)
\left[1-\exp\left(-I(z)\vartheta\alpha/\mu_i^0\right)\right].
\]

The factors are computed independently at lines 233-237 so nutrient/quota and
light limitations can be saved separately.

Eq.9 is evaluated at line 232 with a cumulative discrete approximation to
\(\int_0^z\sum_i P_i(\xi)d\xi\). Parameters `kw`, `kp`, `theta`, and `alpha`
come from config and metadata.

Eq.10, \(dQ_i/dt=V_i^{max}N/(k_i+N)-\mu_iQ_i\), is line 256.
Eq.11, \(P_i=v_iQ_i\), is used in Python to derive cell number as P/Q.
Eq.12, \(C_i=0.0109V_i^{0.991}\) pg C cell-1, is in
`src/splas/traits.py` and the Fortran trait routine at line 37.

## Equations 13-15: grazing and preference

Holling-II grazing (Eq.13) is evaluated at lines 239-248. The log-Gaussian
preference (Eq.14) is line 42, centered on the Wirtz non-isometric optimum
(Eq.15), line 40. Half saturation and preference breadth are configuration
parameters.

## Allometries

`model_traits` at lines 17-46 reproduces the historical volume-form
allometries, algebraically equivalent to the size expressions on thesis
PDF pp.11-12:

- \(V_i^{max}=9.10\,10^{-12}V_i^{0.67}\);
- \(Q_i^{min}=1.36\,10^{-12}V_i^{0.77}\);
- \(k_i=0.1x_i\);
- \(w_i=0.0041x_i^{1.17}\);
- \(\mu_i^0=2.6x_i^{-0.45}\);
- \(I_j=26x_j^{-0.4}\).

## Equations 16-20: forcing and critical depth

Seasonal Eq.16 is implemented in `forcing_value` and
`src/splas/forcing.py` as
\[
I_s(t)=I_0[0.5+0.5\sin^2(2\pi t/365)].
\]

Post-processing implements:

- Eq.17: deepest layer with \(I(t,k)>0.01I(t,0)\);
- Eq.18: \(\log(B_{0:k,t+1}/B_{0:k,t})/\Delta t\);
- Eq.19: interpolated depth where that cumulative rate crosses zero;
- Eq.20: interpolated irradiance at critical depth.

These functions are pure and tested in `tests/test_analysis.py`. Optical
euphotic depth/MLD* remains distinct from the fixed 200 m or 100 m nominal MLD.
