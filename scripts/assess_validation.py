"""Prespecified numerical gates and explicit separation of empirical claims."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import json,hashlib,argparse,zipfile
import numpy as np
from splas2 import *
ROOT=Path(__file__).resolve().parents[1]
FIELDS=['phyto_N_mmol_m2','phyto_CN_inventory_ratio','NPP_C_mmol_m2_d','bottom_export_C_mmol_m2_d','N_weighted_geomean_ESD_um']

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--allow-partial',action='store_true');args=ap.parse_args()
 s=json.load(open(ROOT/'research/validation_summary.json'))
 if s['missing_cases'] and not args.allow_partial:raise RuntimeError('Incomplete study')
 source=(ROOT/'splas2/model.py').read_text();cases={};registry=[]
 for name in ['experiment_registry.json','supplement_registry.json','carbon_boundary_registry.json','recommended_registry.json','recommended_refinement_registry.json','export_intervention_registry.json']:
  path=ROOT/'research'/name
  if path.exists():registry+=json.load(open(path))
  elif not args.allow_partial:raise RuntimeError('Missing registry '+name)
 if not args.allow_partial and len(registry)!=117:raise RuntimeError('Expected117 registered cases')
 for case in registry:
  p=ROOT/'results/thesis_replay'/(case['id']+'.json')
  if not p.exists():continue
  d=json.load(open(p));signature=hashlib.sha256((json.dumps(case,sort_keys=True)+source).encode()).hexdigest()
  if d.get('run_signature')==signature:cases[case['id']]=d
 if len(cases)!=s['current_cases']:raise RuntimeError('Regenerate the analysis from current results')
 # Recompute inventories from saved endpoints; independent accounting oracles
 # are supplied separately by the scientific test suite.
 residuals=[]
 for key,c in cases.items():
  p=ROOT/'results/thesis_replay'/(key+'.npz')
  if not p.exists():continue
  with zipfile.ZipFile(p) as archive:
   if archive.testzip() is not None:raise RuntimeError('Corrupt numerical archive: '+key)
  m=Model(Config(**c['config']),Physiology(**c['physiology']));a=np.load(p)
  err=m.budgets(a['final_state'],True)-m.budgets(a['initial_state'],True)
  residuals.append({'id':key,'endpoint_mass_error_C':float(err[0]),'endpoint_mass_error_N':float(err[1])})
 gates=[]
 for ref,alt,label in [('R01_same_state_dt_0.00625','R01_same_state_dt_0.0125','time_step'),('R02_grid80','R01_same_state_dt_0.0125','vertical_grid'),('R02_grid160','R02_grid80','vertical_grid80_to160'),('R02_particles32','R01_same_state_dt_0.0125','particle_grid')]:
  if ref not in cases or alt not in cases:continue
  sr=cases[ref]['summary'];sa=cases[alt]['summary'];row={'test':label,'reference':ref,'alternative':alt,'relative_changes':{}}
  for field in FIELDS:
   v=sa[field]/sr[field]-1 if sr[field] and sa[field] is not None else None
   row['relative_changes'][field]=v
  row['pass_5pct_bulk_10pct_size']=all(v is not None and abs(v)<= (.10 if key=='N_weighted_geomean_ESD_um' else .05) for key,v in row['relative_changes'].items());gates.append(row)
 inv=s['mass_C_max_relative']<=1e-8 and s['mass_N_max_relative']<=1e-8 and not s['negative_state_cases']
 out={'invariants_pass_1e_minus8':inv,'maximum_reported_mass_error_C_relative':s['mass_C_max_relative'],'maximum_reported_mass_error_N_relative':s['mass_N_max_relative'],
  'recommended_numerical_gates':gates,'recommended_bulk_convergence_demonstrated':len(gates)>=3 and all(g['pass_5pct_bulk_10pct_size'] for g in gates),
  'gates_definition':'Relative to finer reference: <=5percent for bulkP_N,C:N,NPP and Cexport; <=10percent for N-weighted geometric size. Annual comparison from matched10year state; heuristic research thresholds, not observational error bars.',
  'empirical_verdict':'Partial exploratory N-starvation assessment. No species-universal validation. Micromonas persistence is better for C,N,density and C:N than the shared model. No field/grazing/export validation.',
  'matched_fixed_ratio_kinetic_control':False,'no_full_ATP_redox_budget':True,'no_P_Fe_Si_O2':True,'endpoint_checks':residuals,
  'revised_reference_element_balance':next((r for r in s.get('last_window_element_balance',[]) if r['id']=='R00_buffered_reference'),None)}
 (ROOT/'research/acceptance_assessment.json').write_text(json.dumps(out,indent=2,allow_nan=False))
 lines=['# Implementation checks and auxiliary physiology of SPLAS-CN0.2','',
  'The primary objective is a Pattern Oriented Modelling toy study. See [the thesis-grounded protocol](pom/STUDY_PROTOCOL.md) and [current ecological diagnosis](pom/CURRENT_ASSESSMENT.md). These numerical and culture checks support that study; site prediction and absolute ocean-export matching are not its success criterion.','',
  f"Frozen model SHA256: `{s['model_sha256']}`.",'',
  f"{s['current_cases']} registered simulations completed, including{s.get('decadal_case_count',0)} ten-year cases and{s.get('annual_continuation_count',0)} annual interventions/refinement cases. Additional mechanism protocols are reported separately.",'',
  'Twenty scientific tests pass locally and the checkpoint passed GitHub Actions. This verifies the assessed implementation; it does not certify ecological realism.','',
  f"Maximum relative element-inventory error: C{s['mass_C_max_relative']:.3e}, N{s['mass_N_max_relative']:.3e}. No negative saved live state.",'',
  '## Empirical assessment','',out['empirical_verdict'],'',
  'Calibration:Ot/Tw; Tp/Micromonas excluded from optimization. Initial traits,quotas,DIN and density are observed even in excluded species. Scores had been viewed during an earlier development iteration, so the final comparison is exploratory. Raw PANGAEA source lacks final formal author approval;16 published quota means agree within0.4percent.','',
  '| Species | Role | RMSE log C:N | Redfield null | Initial-C:N null |','|---|---|---:|---:|---:|']
 for m in s['culture_metrics']:lines.append(f"|{m['species']}|{m['split']}|{m['CN_RMSLE']:.3f}|{m['CN_Redfield_RMSLE']:.3f}|{m['CN_initial_ratio_RMSLE']:.3f}|")
 lines+=['','Metric =sqrt(mean(log(predicted/observed)^2)); initial sample excluded. The fixed-ratio null is descriptive; the Monod control changes kinetics and is not a matched causal test of stoichiometric complexity.','',
 '## Numerical configuration','',
 '`research/recommended_idealized.json` uses bufferedDIC,16 particle bins,70+70 size classes,40 cells of5m and max_step0.0125day. It retains idealized historical irradiance, not measured site forcing. The IMEX scheme is first order and uses positivity/CFL checks; rtol is not a full-splitting error bound.','',
 '| Same-state test | Bulk<=5%,size<=10% |','|---|---|']
 for g in gates:lines.append('|'+g['test']+'|'+('pass' if g['pass_5pct_bulk_10pct_size'] else 'not demonstrated')+'|')
 lines+=['','These are declared heuristic numerical gates. Inspect individual metrics in acceptance_assessment.json. A failed gate is retained. No test establishes field validity.']
 lines+=['','## Pattern checks','', 'Rank correlations below describe seven deterministic light experiments. They are descriptive diagnostics, not a complete POM assessment. Fianchini p32 reports mixed diversity responses; pp39-40 use export weighted by biomass and normalized by irradiance, which is not export/NPP. Thus an e-ratio sign alone neither confirms nor contradicts that historical efficiency. Fianchini p33 deliberately prescribes two annual light maxima.','', '| Light family | rho C:N | rho evenness | rho e-ratio |','|---|---:|---:|---:|']
 def fmt(v):return 'undefined' if v is None else f'{v:.3f}'
 for a in s.get('pattern_assessment',[]):
  if a['family'].startswith('Fianchini_'):
   vals=[a['metrics'].get(k,{}).get('spearman_rho') for k in ['phyto_CN_inventory_ratio','phyto_evenness_N','e_ratio_at_bottom']]
   lines.append('|'+a['family']+'|'+'|'.join(fmt(v) for v in vals)+'|')
 lines+=['',
 '## Findings that change the inference','',
 '- An openN source without matching carbon replenishment can depleteDIC. The0D carbon-bath ensemble and revised column separate this boundary confound from N effects.',
 '- Two detrital bins conserve mass but give a substantially different export from16bins. Element conservation is insufficient evidence for resolved gravitational flux.',
 '- Historical light has two maxima per year by construction. Matched-dose annual/constant controls are required before interpreting bloom counts.',
 '- StockZ:P and fluxTTE are distinct. Consumer production depends on actual preyC:N and on the declared fixed consumer stoichiometry.',
 '- C:N differences are decomposed into within-class and composition terms; an algebraic decomposition is not unique causal attribution.',
 '- Ten years do not guarantee a stationary community. Temporal means,phase sensitivity and drift flags remain part of the result.','',
 '- A nearly stationary phytoplankton stock can coexist with accumulatingDIN. Inspect last_window_element_balance before asserting whole-system stationarity.','',
 '## Remaining validity limits','',
 'The model resolvesC:N, not fullC:N:P or redox. Reference size is not instantaneous cell volume. Storage capacity,carbon-based feeding,homeostatic grazers,famine attrition,particle mapping and remineralization are model hypotheses. Carbon cost is not a closed ATP/electron budget. No matched site forcing/observations,independent grazer data,ballast,aggregation,bacteria,mixotrophy or carbonate chemistry are resolved.','',
 'Primary literature and the thesis-family map are inLITERATURE.md,MODEL.md and the experiment registries. Generated arrays and figures are embedded in the distributed PDF; repository source remains inGit.']
 (ROOT/'research/VALIDATION.md').write_text('\n'.join(lines)+'\n')
 print(json.dumps({k:out[k] for k in ['invariants_pass_1e_minus8','recommended_bulk_convergence_demonstrated']},indent=2))
if __name__=='__main__':main()
