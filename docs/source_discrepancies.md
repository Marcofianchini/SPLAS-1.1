# Source discrepancies and implemented decisions

The thesis is authoritative for scientific formulation; the historical code is
used to reconstruct unstated implementation details. No discrepancy below is
silently resolved.

| Topic | Thesis | Historical code | Implemented decision | Expected impact / alternative |
|---|---|---|---|---|
| Open supply S | Table 3 (PDF p.16) gives closed/open `0/1`; section 3.2 (PDF p.17) explicitly says `S=4` was used | `input_default.json` gives 4 | `thesis_consistent` uses 4, following the experiment narrative | `legacy_current` uses 1 to retain the table/code sensitivity alternative |
| Interior supply profile | Eq.3 (PDF p.8): 0 at surface, `0.1 S` for every interior layer, `S` at bottom | `solver.f90` multiplies by layer index `k` | Fixed 0.1 S in every interior layer | `legacy_reconstructed` retains `k * 0.1 S` for diagnosis |
| Time step | Chapter 3 (PDF p.16) says two hours | JSON uses 0.012 d (17.28 min); text in section 2.3 says updates are usually 30 min | 1/12 d for thesis suites | Legacy profile uses 0.012 d; neither is relabeled as two hours |
| Duration/output cadence | Ten years; no unambiguous output cadence | 3652 d, `print_step=300` counts steps | 3652 d and 4-day output, sufficient for seasonal blooms without impractical uncompressed NetCDF volume | Legacy uses 3650.4 d and 3.6 d output so dt/output/duration are integral |
| Closed remineralization | Table 3 and section 3.2: 0.5 d-1 | Alternate solver uses 0.40 | 0.5 d-1 | Legacy closed config retains 0.40 |
| Closed zoo recycling γ | Table 3 and text: 1 | Alternate solver uses 0.8 | 1 | Legacy closed config retains 0.8 and is not claimed mass-conserving |
| Closed system selection | No transport, no supply, isolated layers | Requires swapping `solver.f90`; JSON name has no effect | One core with `system.mode=closed` setting supply/diffusion/P sinking/OM sinking to zero | Eliminates source swapping |
| Dissolved-N loss and quota | Eq.6 uses cellular uptake; Eq.4 evolves P biomass with realized Droop growth; text says closed mass is constant | Literal code follows both equations but does not expose a budget | Recommended profile couples dissolved-N loss to realized P production so N+P+Z+OM closes; quota still controls realized growth | Legacy profile uses the literal historical uptake term. This is the most consequential unresolved formulation choice |
| Initial cellular quota | Not reported | `initial_p_quota=0.1` is multiplied by 1e-11, often below Qmin | Start at 1.1 Qmin | Multiplier is configurable and recorded |
| Nutrient pools | Equations describe one prognostic N pool | Code allocates 3 and assigns the same scalar tendency to N/P/C | Only N is prognostic; P and C coordinates are diagnostic Redfield equivalents | Prevents false claims of explicit P/C cycling |
| Self-shading quadrature | Eq.9 integrates P from surface to z | Historical function treats the first interior layer specially and sums discrete P inconsistently | Cumulative layer-center quadrature `sum(P)*dz` | Changes near-surface attenuation slightly |
| Seasonal forcing | Eq.16 and text: minimum I0/2, two equal annual maxima | `0.5 + 0.5 * supply_season`, with `sin^2` | `I0 * (0.5 + 0.5 sin^2(2πt/365))` | Matches both sources |
| Domain depth | Eqs.1-2: center `z_k=(k-0.5)dz`, nominal depth `h=(M-1)dz` | 41 centers extend to 202.5 m while h is 200 m | Coordinate stores layer centers; `h=200 m`/`100 m` is stored and discussed as the nominal MLD/domain boundary | Avoids conflating MLD with the last cell center |
| Diversity naming | Thesis discusses Shannon index | Historical routine computes normalized evenness | Save entropy, evenness, present-class count, and validity mask separately | Removes semantic ambiguity |
| C:N | P uses Montagnes size-carbon; detritus uses Redfield | Same, but units are incompletely documented | Both are retained with explicit variable units; Redfield detritus is labeled a model limitation | Exact C:N remains sensitive to the unreported initial quota |
| Detrital C mass conversion | Redfield specifies a molar C:N of 106:16 | Historical `det_c_export` multiplies mmol N by 14 after the molar ratio | Multiply mmol C by 12 mg mmol-1 | Legacy 14-factor result is documented but not used as a physical carbon conversion |
| Critical-depth rate | Eq.18 typography is ambiguous about multiplication/division by Δt | No pipeline exists | Use log cumulative-biomass change divided by elapsed days, and linearly interpolate the r=0 crossing | Operational interpretation is tested on synthetic profiles |
| MLD* | Eq.17 uses the 1% euphotic depth as MLD* | No post-processing exists | Store euphotic depth separately from fixed nominal MLD | Prevents calling the optical diagnostic a physical MLD |
| Export efficiency | Text says export "weighted by biomass" and divided by light; Fig.4.22 units are mgC:W | No unique implementation | Save export/biomass, export/light, and export/(biomass*light); the figure uses export/light to match units | Final identification requires original numerical arrays or digitization |
| Random forcing | Normal light option has no seed | NumPy global RNG | Configurable seed, deterministic Fortran Park-Miller/Box-Muller generator, seed in metadata | Bitwise identity is compiler-sensitive but sequence is explicit |
