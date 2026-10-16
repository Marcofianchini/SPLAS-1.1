# SPLAS-CN 0.2: mathematical and epistemic specification

This is a new Python implementation beside the unchanged SPLAS-1.1 Fortran,
base commit `76f992c9924420c715a660d017387882efe992dd`. It is a conservative
reinterpretation of thesis experiments, not a bitwise reproduction or a claim
that historical figures have been independently validated. The attached
MoloneyField 0D v0.1 proposal is the developmental starting point; Moloney's
1988 thesis itself used fixed group C:N in **mass**, not independent reserves.

## Units and state

Time day; depth m; ESD micrometres; cell volume micrometres^3; PAR micromol
photons m^-2 s^-1. Chemical concentrations mmol element m^-3. `a_i` is in
10^8 cells m^-3. Per-abundance quotas are mmol element per10^8 cells.

Each phytoplankton class has abundance a, C-only reserve R_C, and an organic
N reserve R_N carrying carbon at rho_R=4 molC/molN. Structural quotas c0,n0
are constant per cell. Actual biomass is

    B_C = c0*a + R_C + rho_R*R_N
    B_N = n0*a + R_N.

The organic N reserve is a reduced functional/storage pool, not explicit
protein/RNA biochemistry. With c0/n0>=rho_R and nonnegative state, C:N>=rho_R
is an algebraic invariant. This prevents the carbon-free N accumulation
possible when a free-standing N pool is mistaken for carbon-containing biomass.
An inorganic-storage alternative rho_R=0 is available only as an explicit
hypothesis; it changes this lower-bound interpretation.

Other pools: DIN, DIC, DOC, DON, carbon biomass of each zooplankton class,
and independently transported/remineralized detrital C,N in particle bins.
Zooplankton N=Z_C/rho_Z, rho_Z=5.5 by default. Thus consumer homeostasis is an
assumption, not dynamic consumer stoichiometry. P, Fe, Si and O2 are not states.
No C:N:P or iron/silicate limitation validation is claimed.

One mmolC=12mgC; one mmolN=14mgN. A mass C:N is multiplied by14/12 to become
molar. The 10^8 abundance conversion gives q_C=pgC/120 and q_N=pgN/140.
Reference carbon is Montagnes etal1994 `0.109*V^0.991 pgC/cell`; the historical
0.0109 coefficient is incorrect relative to that source. Applying the relation
to1-200um extrapolates outside the original2-60um sample and is uncertain.
Reference N=q_C/6.625 in community experiments, measured replete N in cultures.
This reference ratio initializes traits; it does not constrain future C:N.

## Physiology

The following are hypotheses, with numerical defaults in `Physiology`:

    mu_cap = mu_scale * 2.6 * ESD^(-0.45)
    Vmax_N = uptake_scale * 9.1e-4 * V^0.67
    n0 = Nref*clip(nmin_fraction*(ESD/5)^nmin_exponent, .05, .8)
    c0 = min(structural_CN*n0, .85*Cref)
    Cmax = Cref*clip(cmax_ratio*(ESD/5)^cmax_exponent,1.05,8)
    nstatus = min(1, B_N/(Nref*a))
    theta = theta_min+(theta_max-theta_min)*nstatus/(1+PAR/acclimation_PAR)
    p_specific = photo_scale*mu_cap
    fT = Q10^((temperature-18)/10)
    F = fT*a*Cref*p_specific*nstatus*(1-exp(-alpha*theta*PAR/(p_specific*nstatus)))
        * DIC/(K_DIC+DIC) * max(0,1-(R_C/(Cmax*a-c0*a-rho_R*R_N))^2).

At18C the low-PAR slope per reference carbon is alpha*theta; at other
temperatures it is fT*alpha*theta. The universal Q10 thus changes both
photosynthetic maximum and low-light slope in this formulation. This is an
explicit temperature hypothesis, not a demonstrated acclimation law.
nstatus must occur in the saturation
**denominator** as well as maximum, otherwise low-PAR N limitation is counted
twice. Theta is diagnostic instantaneous acclimation, not a prognostic pigment
pool; photoinhibition and acclimation lags are absent.

Potential N uptake is a*Vmax_N*DIN/(0.1*ESD+DIN), reduced by luxury-uptake
feedback max(0,1-(B_N/(2*Nref*a))^2). Actual uptake U is capped by
kappa*R_C/(cost_N+rho_R): assimilation both transfers rho_R*U carbon into the
organic N pool and respires cost_N*U carbon. Thus dark uptake can use existing
energy reserves but cannot assimilate N at zero available carbon reserve.

Division J=a*mu_cap*min(R_C/B_C,R_N/B_N), or a declared product-colimitation
alternative. J is capped by reserve mobilization for both required elements.
An extra cost_growth*c0*J goes to DIC. Carbon needed from R_C for structural
synthesis is (c0-rho_R*n0)*J; the rest already accompanies the N reserve.

Maintenance demand is maintenance*c0*a. The paid fraction draws R_C. A deficit
induces whole-cell attrition at deficit/B_C, respiring all C of that fraction
and returning its N to DIN. This famine closure is an explicit reduced
hypothesis, not a claim of resolved autophagy. Death/exudation release the
instantaneous actual C:N, including carbon accompanying released organic N.
All physiological capacities except fixed exudation rates scale by the declared
Q10 temperature factor. One universal Q10 is a prior, not a temperature fit.

## Grazing and recycling

For predator j, available carbon A_j=sum_i(phi_ij*B_Ci). Correct matrix
contraction is phi.T @ B_C. Carbon ingestion of i is

    G_ij = Imax_j*Z_Cj*phi_ij*B_Ci/(K_C+A_j), Imax=26*ESD_Z^(-0.4).

Total ingestion cannot exceed Imax*Z_C. Prey abundance and both reserves are
removed in the same fractional amount, so N ingestion follows actual prey C:N.
Carbon and nitrogen assimilation efficiencies are distinct. Consumer production
is min(eta_C*I_C/(1+cost_growth), rho_Z*eta_N*I_N). Assimilated excess C is
respired and excess N excreted; unassimilated C,N enters detritus.
Basal consumer respiration removes homeostatic biomass C and returns its
paired N to DIN; consumer C reserves are not separately represented.
Community-dependent mortality xi*(sum Z_N)*Z_C and own quadratic mortality
xi*Z_Nj*Z_Cj are separate hypotheses. Recycled fraction goes to detritus;
nonrecycled fraction is an explicit higher-trophic export ledger. This export
is not gravitational carbon flux across the reference depth.

Death and egested particles are allocated conservatively to log-size detritus
bins, with no elemental mass creation by rebinning. Speeds scale with particle
ESD around8m/day at50um; dead cells use cell ESD and egesta use0.1*predatorESD.
These particle-size mappings are uncertain priors. Two bins are not presumed
converged. Lysis to DOM and remineralization to inorganic pools are separate
first-order stages. Setting both rates0.15/day is not identical to a single
historical0.15/day remineralization stage. No bacteria, aggregation, fragmentation,
ballasting, oxygen limitation or refractory DOM are resolved.

## Transport, boundaries and light

Finite-volume cell averages, no ghost cells; layer thickness H/M.
Face flux of each live tracer is `-D_face*(Y[k+1]-Y[k])/dz + w*Y[k]`.
The same flux is subtracted/added at the two neighboring cells. D in m2/s is
converted by86400 into m2/day. Carbon, N reserves and abundance have identical
cell sinking speeds. Detrital C,N share the same particle speed. At the bottom,
only explicitly configured sinking exports mass; diffusion has zero boundary
flux. When export_bottom=False this version disables all sinking. This is an
explicit isolated/reflected-transport assumption, not internal sinking plus
reflection.

Volumetric DIN source, bottom areal DIN flux, dilution and optional DIC
relaxation have independent external input/output counters. DIC relaxation is a
boundary hypothesis, not a carbonate equilibrium/air-sea CO2 model. The external
reservoir acts volumetrically in every cell as `k*(DIC_external-DIC)`,
not as a physical sea-surface gas flux. The 0D
community counterpart uses it to avoid an unintended closed-carbon resource
experiment during decades of N addition. Closed cases have no exchange.
External inputs need not balance outputs: a nearly stationary phytoplankton
stock can coexist with accumulating DIN. Whole-system stationarity must be
assessed from element retention and inorganic inventories as well as biomass.

PAR uses Beer-Lambert attenuation by water and diagnostic chlorophyll, with
three local fixed-point iterations. A homogeneous0D box receives prescribed
PAR directly and has no fictitious optical attenuation across its inventory
depth. 1D calculations attenuate across cell centers/faces. Historical
`I0*(.5+.5*sin(2*pi*t/365)^2)` has period182.5d and mean.75I0. Annual and
constant controls have the same annual PAR integral. Diel square waves align
all numerical stages to half-day boundaries and allow independent phase.
DePasquale's `sin^2(pi*t/365)` nutrient source has period365d, despite the
thesis text calling it two years. Both formula and interpretation are recorded.

## Numerical method and conservation tests

SSPRK3 with embedded RK2 local-error estimate. A step is rejected if a forward
Euler substage leaves the positive orthant or normalized local error exceeds1.
No post-integration clipping and no extinction cutoff. Defaults rtol2e-4,
atol1e-8 in chemical units, abundance tolerance scaled by structural quotas.
Diagnostics are integrated with the same Runge-Kutta stages as live states.
Mass check is inventory + bottom export + HTL export + dilution output - all
external inputs, separately for C and N. This identity is tested from an
independent sum over biological inventories, random states, open and closed
boundaries, and an analytic passive-diffusion eigenmode.

Positivity alone was insufficient for a rapidly recycled Monod R* experiment:
a positive oscillatory numerical equilibrium was found during development.
The dynamic solver now controls local error; the historical didactic fixed-C:N
solver additionally caps the step by **gross DIN turnover**, not net dDIN.
The analytic R* oracle and time/tolerance refinement are required evidence,
not replaceable by mass closure.

## Inference and validation limits

Culture trajectories constrain nutrient starvation in4 taxa. Ot/Tw estimate
shared parameters; Tp/Ms are never inputs to the optimizer. Initial C,N,
cell abundance,DIN and day0 replete quotas are observed even for withheld
species: these are conditional predictions, not predictions from size alone.
The early development iteration inspected Tp/Ms scores before the kinetic
revision, so the withheld-species comparison is exploratory, not a pristine
prospective holdout. An independent dataset is still necessary.

The dataset was published in2026 with a PANGAEA notice of missing final formal
author approval. The16 replete/starved quota means cross-check within0.4percent
against2019 Table1. Two known bacterial-contaminated late stages are excluded.
At Ot day2 quotas repeat day0 exactly; uncertain time/date labels and DIN unit
label are documented in extraction provenance. Laboratory85umol PAR12:12,
18C except Micromonas6C, cannot validate ocean grazing/export or a globalQ10.

No thesis graphs are fitted as truth. Replays compare mechanistic responses,
report failed patterns and separate numerical from structural uncertainty.
Stock Z:P is reported as stock ratio. Flux TTE=Z carbon production/NPP over
one shared window. e-ratio=gravitational C flux at declared bottom/NPP integrated
over the same upper column/window; it is not export to the deep ocean or
sequestration. NPP=GPP-phyto respiration including assimilation/maintenance/growth.
A negative NPP window makes efficiency undefined. Short transient efficiencies
can exceed1 by drawing on initial stocks; no automatic clipping to[0,1].
Shannon evenness refers to size classes and weight element; it is not taxonomic
biodiversity. Spectra use biomass per delta log10ESD. Initial inventories stay
constant under grid refinement. Euphotic depth is optical, whereas prescribed
mixing depth is turbulent. A frozen-state depth integral of biological net-C
growth is explicitly distinguished from a local compensation depth, modeled
bloom onset, and the assumptions of Sverdrup's well-mixed critical depth.

### Stiff nutrient-assimilation option

Long ensembles set `method='imex'`. SSPRK3 advances the non-assimilation
processes with positivity checks and a conservative transport CFL bound; at each accepted step a
backward-Euler solve advances N assimilation. The **combined method is first
order**, despite its SSPRK3 nonstiff component. The IMEX branch does not apply the embedded error controller to its split
subproblem: rtol/atol control the explicit reference only. `max_step`
controls the IMEX discretization error. Cultures retain full explicit integration as an
independent reference. Ensemble step refinement is therefore mandatory.

For each class, abundance is fixed over the implicit substep. With per-cell
old N quota b, maximum2Nref=q, `r=Vmax*DINnew/(K+DINnew)`, the implicit uptake u
solves `u=r*(1-((b+h*u)/q)^2)`, with nonnegative feedback, and is bounded by
`kappa*CRold_per_cell/[(cost_N+rho_R)*(1+kappa*h)]`. A stable quadratic formula
and bracketed Newton solve the shared DIN equation
`DINnew+h*sum(a_i*u_i)=DINold`. The final lower bracket and an ulp safety margin
prevent donor overdraw; every C/N change uses that same uptake amount.
This avoids chemical stiffness without imposing extinction thresholds or
altering a completed state by clipping. Compare against the full explicit
solver and refine h; do not interpret the implicit solver's rtol as a bound on
full splitting error.


### Further scope limits

The positive carbon charge for N assimilation is energetic bookkeeping, not a
closed ATP/electron/redox budget. DIN lumps nitrate and ammonium; the fixed
cost_N=1.5 is a prior, not an identified or universal substrate-specific cost.
No claim of thermodynamic completeness follows from atomic C/N conservation.
The fixed-C:N Monod example and fixed-ratio culture nulls do not isolate the
effect of stoichiometry while retaining identical kinetic responses. A matched
conservative fixed-ratio counterpart is still required for that causal claim.


### Trait and encounter limits

The size class is a reference structural trait, not a prognostic measurement
of instantaneous cell volume. Carbon storage can alter density/volume in real
cells, but the grazing kernel and sinking law here do not adapt to that change.
Large Cmax values are weakly constrained by starvation data; a fitted capacity
is not proof of a physically possible carbon density. Carbon-based saturation
of feeding, rather than cell-number/biovolume-based encounters, is a specific
food-web hypothesis and is not empirically validated by the culture test.

Given the community parameters, the modeled net grazer production per unit
ingested carbon is bounded by
`min(0.396/1.2, 5.5*0.67/(prey C:N))`. The transition to N-limited growth is
at prey C:N~11.17. At C:N45 the bound is~0.082, compared with0.33 at6.625.
This is a consequence of the model assumptions, not a measured trophic
efficiency; consumer adaptation, prey biochemistry and compensatory feeding
could change the relation.


### Compatibility of inherited allometries

The ingestion curve26ESD^(-0.4) is an inherited modeling prior. Hansen1997
reports group differences and no satisfactory overall cross-group regression;
a successful saturation bound does not empirically validate this rate law.
The available carbon-ingestion table needs matched ESD/temperature/trophic
identity before it can test the specific ESD rule; no guessed body conversion
is used as an empirical fit here.

Likewise, independent trait regressions are not automatically jointly viable.
With corrected MontagnesC and a fixed6.625 reference C:N, the historical
Qmin=1.36e-12V^0.77 mmolN/cell exceeds referenceN in the smallest bins.
`allometry_admissibility.json` records that algebraic check. The new model
uses fitted minimum-quota fractions and explicitly does not claim to retain
or validate that historical absolute quota curve. This change matters when
interpreting differences from the original size spectra.
