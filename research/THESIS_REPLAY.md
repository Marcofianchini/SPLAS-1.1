# Experiments: De Pasquale → Fianchini → SPLAS-CN

This is a traced reinterpretation of experiment families with corrected units,
independent elemental inventories and new physiology. It is not a compiled
replay of every historical figure or an assertion of source-level identity.
Historical Fortran is retained unchanged. The model/case SHA256 signature
rejects stale numerical results.

| Development step | Thesis question | Executable cases | Interpretation |
|---|---|---|---|
| Moloney–Field benchmark | Fixed group C:N and size organization | Attached MoloneyField0D v0.1; `splas2/npz_fixed.py` | The attachment is a developmental benchmark, not an authenticated archival1991 release. Moloney1988 documents mass ratios. |
| DePasquale2.6 | Smallest competitor; sequential predator assembly | `run_mechanisms.py`:0,1,2,4,8,16,32,70 predators; D11 quota counterparts0–16 | Analytic R* is an oracle; Monod vs reserve dynamics also changes kinetics. |
| DePasquale3.3–3.4 | Nutrient and grazing control; closure | D00,D01,D02,D03 | Vary N supply and grazer mortality/closure, preserving initial total inventory. |
| DePasquale3.5.1–3.5.4 | Selectivity, phytoplankton losses, resolution, yield | D04–D09 | New C-based ingestion and element-specific consumer constraints differ from the historical fixed-yield interpretation. |
| DePasquale3.5.5 | Periodic nutrient source | D10, two source phases and matched-mean constant input | Algebraic period is365days; phase is tested rather than treating one endpoint as equilibrium. |
| DePasquale4 | Light, vertical structure, self-shading, mixing | F00,F04,F05; R05 | Conservative transport in real cells; variable face diffusivity. R05 is an annual attenuation intervention. |
| DePasquale5 | Sinking, remineralization and export | F08,F10; R04 | R04 varies N, two-stage recycling, cell sinking and consumer yield from one common state. No new decadal equilibria are claimed for these annual interventions. |
| Fianchini3 | Open and closed reference columns | F00_open/closed | Source and bottom flux in declared units; closed case has isolated layers and no sinking. |
| Fianchini4.2 | Absolute light0.7–1.3× | F01 open/closed plus F00 at1× | Seven light levels; C:N and evenness, export and e-ratio checked separately. |
| Fianchini4.3 | Seasonally varying light | F02 | Historical two-maxima/year formula, annual formula and constant PAR share the same annual light dose. |
| Fianchini5 | Shallower column and critical-depth hypothesis | F03;18 independent mixing-depth mechanisms | Column depth100m is a boundary. Mixing is prescribed independently of optical Zeu; integrated frozen-state critical-depth diagnostic has a limited interpretation. |
| Numerical discrimination | Time, physical grid, size and particle bins | D07,D12,F07,F08,F09,R01,R02 | Same integrated initial inventories. R01/R02 share the physical ten-year endpoint and compare365-day continuations. |
| Carbon-boundary discrimination | N enrichment versus changing C limitation | D14 plus R00/F06 | A declared C reservoir separates resource-boundary effects. This is not carbonate/air-sea chemistry. |
| Physiological discrimination | N starvation and recovery | Four Liefer taxa;25 light×DIN cases;16 light phases; three recovery/chemostat protocols | Culture scores are empirical but exploratory; idealized protocols verify mechanisms without supplying independent observations. |

The historical200→100m contrast also changes total areal N input from100 to60
mmolN/m2/day (`0.4*H+20`) and truncates the physical initial profile. It is not
a single-factor test of turbulent mixing depth. The18 mechanism experiments
retain one domain and prescribe mixing independently of optical depth.

## Order and evidence

Run the commands in the main README. The finalized study contains117 registered
long/annual community cases,67 additional assembly/light/mixing cases and three
recovery protocols; four culture trajectories are evaluated separately.
`validation_summary.json` records numerical summaries, deterministic pattern
checks, temporal drift, matched-dose contrasts and C:N decomposition.
`acceptance_assessment.json` declares heuristic refinement gates and retains
failed gates. `VALIDATION.md` states what the evidence supports.

A deterministic monotone response or high linear R2 demonstrates a pattern
within these model assumptions. It does not establish causation in the ocean.
A missing historical pattern is retained; parameters are not fitted to thesis
curves to make that pattern reappear. No field/export data were used as a
validation target. Alternative particle, grazing and allocation mechanisms
remain scientifically open.

The failed40-to80 vertical gate triggered one additional160-cell annual run.
The80-to160 comparison retains the same5percent bulk/10percent size criteria;
the failed40-cell result remains in the assessment.
