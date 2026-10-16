# Primary literature and model decisions

Verification date:2026-10-16. This is a targeted critical comparison, not a
systematic review and not an assertion of complete coverage. A newer model's
existence is not empirical validation of SPLAS. Numerical tests, culture scores,
qualitative ecological consistency and field validation are distinct evidence.

| Primary source | Relevant evidence | Consequence for SPLAS-CN | Test/domain limitation |
|---|---|---|---|
| Moloney1988 doctoral thesis, University of Cape Town, https://hdl.handle.net/11427/22132, pp46-47 and AUTSUB pp239-240 | Fixed group C:N in mass, C fixation linked to N assimilation. | Distinguish archival formulation from attached MoloneyField0D proposal and new dynamic model. |1991 article not fully available; no source-level equivalence claim. |
| DePasquale2018 thesis, https://arts.units.it/bitstream/11368/2929831/2/M_DePasquale_TESI_FINAL.pdf | Mechanism hierarchy: Monod assembly, quota size spectra, nutrient and grazer control, light/transport/export. | Trace each experiment family and preserve negative pattern results. | Different dynamics and corrected units preclude pixel matching as validation. |
| Fianchini2021 thesis, user-supplied corrected PDF | Open/closed columns, irradiation, seasonal light and critical-depth interpretation. | Replay stated contrasts; separate source-code discrepancies, optical from turbulent depths. | Original exact run archives absent; absolute irradiance ambiguity remains. |
| Montagnes etal1994, https://doi.org/10.4319/lo.1994.39.5.1044, author PDF https://bpb-us-w2.wpmucdn.com/sites.uwm.edu/dist/6/103/files/2016/07/Montagnes_etal_1994-1wjozj1.pdf | Published carbon coefficient0.109, exponent0.991; large interspecies and measurement uncertainty. | Correct historical0.0109. Keep measured reference quotas in culture tests. |2-60um source species do not establish validity of1-200um extrapolation. |
| Liefer etal2019, https://doi.org/10.3389/fmicb.2019.00763 | C and N quotas diverge during N starvation; species differ in carbon accumulation. | Test changing C,N separately, not C:N alone. Exclude known contaminated late stages. |Four taxa, N-starvation and fixed lab light; no oceanic export validation. |
| Liefer etal2026 dataset, https://doi.org/10.1594/PANGAEA.900248 | Public raw2016 culture experiments associated with2018/2019 papers; CC-BY4.0; missing final formal author approval. | Archive hashes and extraction.16 published quota means agree within0.4percent. Ot/Tw fit; Tp/Ms excluded from estimation. | Source caveat persists after a table cross-check. Earlier developmental use means evaluation is exploratory. |
| Omta etal2024, https://doi.org/10.1016/j.jtbi.2024.111883, author PDF https://mmab.ca/pubs/omta-2024-jtb.pdf | Dynamic allocation and storage explain stoichiometric variation, but storage parameters differ across groups and remain uncertain. | Keep a reduced functional organic-N pool with explicit accompanying C and a C-only reserve. Treat transferable size-only physiology as testable. | Uses the same Liefer cultures: agreement with Omta is not an independent data replication. SPLAS does not reproduce its full macromolecular or P model. |
| Harcourt etal2024, https://doi.org/10.1371/journal.pone.0292337 | Strain differences and thermal adaptation modulate elemental quotas in nutrient-replete Synechococcus. | Universal size scaling/Q10 are hypotheses; an optimum at6C need not imply universally cold-suppressed physiology. | Qualitative domain challenge only; these raw strain data were not fitted or numerically scored here. |
| Mollica etal2025, https://doi.org/10.1002/lno.70058 | Size-fractionated Baltic observations show small fractions can dominate community stoichiometry. | Report C- and N-weighted size spectra and separate within-class stoichiometry from composition shifts. | Different environment and no available matched forcing dataset; no direct Baltic hindcast or statistical validation. |
| Karakus etal2025, https://doi.org/10.1029/2024GB008403 | Export from the surface and transfer to greater depth respond differently to sinking/ballast/remineralization assumptions. | Keep bottom flux, export efficiency, higher-trophic loss and sequestration distinct. Test particle-bin sensitivity. | No ballast/O2/viscosity or deep column in SPLAS; these mechanisms are unresolved. |
| Dogliotti etal2026, https://doi.org/10.5194/bg-23-1859-2026 | Bloom phenology depends jointly on light, independently inferred mixing and regional conditions. | Prescribe turbulent mixing depth independently of Zeu; grazer ablations and integrated growth-loss diagnostics. | No Argentine shelf site forcing or data assimilation; mechanism tests only. |
| Rogerson etal2026, https://doi.org/10.5194/bg-23-3179-2026 | PISCES physiological configurations change NPP and export differently; complexity does not monotonically improve all outputs. | Use explicit flux definitions and numerical/structural ensembles. Do not infer more export from higher NPP or model complexity. | A global model intercomparison is not an observational target for an idealized SPLAS box. |
| Papapostolou etal2026, https://doi.org/10.5194/gmd-19-6627-2026 | NUM combines size-based food webs, metabolic costs, trophic strategies, silicate traits and multiple physical setups. | Retain metabolic costs and trait-based structure while identifying absent mixotrophy, bacterial pathways and ontogeny. | No claim that current SPLAS has NUM's validation or represented mechanisms. |
| Bernish etal2026, https://doi.org/10.3389/feart.2026.1839921 | Cell allocation under Fe limitation responds to Fe/light/temperature interactions. | Explicitly restrict the N-only nutrient model's validity domain. | P,Fe,Si limitation cannot be established by changing a universal N half-saturation constant. |

## Decision rule

Accept a numerical implementation only after independent accounting, positivity,
analytical/independent-solver oracles and refinement. Accept a physiological
claim only in the assessed culture domain and after comparison with simple
baselines. A modeled pattern similar to a thesis figure supports a possible
mechanism; it does not establish unique causation. A failed pattern is retained
in the results, not repaired by fitting to the historical graph. Field
validation would require location-specific forcing and matched C,N,size,
production and particle-flux observations withheld from estimation.
