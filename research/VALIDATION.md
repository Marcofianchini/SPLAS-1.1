# Implementation checks and auxiliary physiology of SPLAS-CN0.2

The primary objective is a Pattern Oriented Modelling toy study. See [the thesis-grounded protocol](pom/STUDY_PROTOCOL.md) and [current ecological diagnosis](pom/CURRENT_ASSESSMENT.md). These numerical and culture checks support that study; site prediction and absolute ocean-export matching are not its success criterion.

Frozen model SHA256: `e7de66b3983742fc84e8114862c9a3c8fe9969fe950dc303e720395fcb6c3f85`.

117 registered simulations completed, including105 ten-year cases and12 annual interventions/refinement cases. Additional mechanism protocols are reported separately.

Twenty scientific tests pass locally and the checkpoint passed GitHub Actions. This verifies the assessed implementation; it does not certify ecological realism.

Maximum relative element-inventory error: C1.617e-11, N4.229e-09. No negative saved live state.

## Empirical assessment

Partial exploratory N-starvation assessment. No species-universal validation. Micromonas persistence is better for C,N,density and C:N than the shared model. No field/grazing/export validation.

Calibration:Ot/Tw; Tp/Micromonas excluded from optimization. Initial traits,quotas,DIN and density are observed even in excluded species. Scores had been viewed during an earlier development iteration, so the final comparison is exploratory. Raw PANGAEA source lacks final formal author approval;16 published quota means agree within0.4percent.

| Species | Role | RMSE log C:N | Redfield null | Initial-C:N null |
|---|---|---:|---:|---:|
|Ot|calibration|0.111|0.585|0.666|
|Tw|calibration|0.150|1.885|1.364|
|Tp|held_out|0.259|1.025|1.117|
|Ms|held_out|0.188|0.609|0.165|

Metric =sqrt(mean(log(predicted/observed)^2)); initial sample excluded. The fixed-ratio null is descriptive; the Monod control changes kinetics and is not a matched causal test of stoichiometric complexity.

## Numerical configuration

`research/recommended_idealized.json` uses bufferedDIC,16 particle bins,70+70 size classes,40 cells of5m and max_step0.0125day. It retains idealized historical irradiance, not measured site forcing. The IMEX scheme is first order and uses positivity/CFL checks; rtol is not a full-splitting error bound.

| Same-state test | Bulk<=5%,size<=10% |
|---|---|
|time_step|pass|
|vertical_grid|not demonstrated|
|vertical_grid80_to160|not demonstrated|
|particle_grid|pass|

These are declared heuristic numerical gates. Inspect individual metrics in acceptance_assessment.json. A failed gate is retained. No test establishes field validity.

## Pattern checks

Rank correlations below describe seven deterministic light experiments. They are descriptive diagnostics, not a complete POM assessment. Fianchini p32 reports mixed diversity responses; pp39-40 use export weighted by biomass and normalized by irradiance, which is not export/NPP. Thus an e-ratio sign alone neither confirms nor contradicts that historical efficiency. Fianchini p33 deliberately prescribes two annual light maxima.

| Light family | rho C:N | rho evenness | rho e-ratio |
|---|---:|---:|---:|
|Fianchini_constant_open|-0.964|0.964|1.000|
|Fianchini_constant_closed|0.714|-0.893|undefined|
|Fianchini_semiannual_200m|-0.964|0.643|1.000|
|Fianchini_semiannual_100m|1.000|-0.964|1.000|

## Findings that change the inference

- An openN source without matching carbon replenishment can depleteDIC. The0D carbon-bath ensemble and revised column separate this boundary confound from N effects.
- Two detrital bins conserve mass but give a substantially different export from16bins. Element conservation is insufficient evidence for resolved gravitational flux.
- Historical light has two maxima per year by construction. Matched-dose annual/constant controls are required before interpreting bloom counts.
- StockZ:P and fluxTTE are distinct. Consumer production depends on actual preyC:N and on the declared fixed consumer stoichiometry.
- C:N differences are decomposed into within-class and composition terms; an algebraic decomposition is not unique causal attribution.
- Ten years do not guarantee a stationary community. Temporal means,phase sensitivity and drift flags remain part of the result.

- A nearly stationary phytoplankton stock can coexist with accumulatingDIN. Inspect last_window_element_balance before asserting whole-system stationarity.

## Remaining validity limits

The model resolvesC:N, not fullC:N:P or redox. Reference size is not instantaneous cell volume. Storage capacity,carbon-based feeding,homeostatic grazers,famine attrition,particle mapping and remineralization are model hypotheses. Carbon cost is not a closed ATP/electron budget. No matched site forcing/observations,independent grazer data,ballast,aggregation,bacteria,mixotrophy or carbonate chemistry are resolved.

Primary literature and the thesis-family map are inLITERATURE.md,MODEL.md and the experiment registries. Generated arrays and figures are embedded in the distributed PDF; repository source remains inGit.
