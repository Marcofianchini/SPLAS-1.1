"""Diagnose the C-boundary confound in a nominally N-only 0D comparison."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import copy,json
from concurrent.futures import ProcessPoolExecutor,as_completed
from scripts.run_experiments import registry,execute
ROOT=Path(__file__).resolve().parents[1]

def main():
 reference=next(c for c in registry() if c['id']=='D00_reference');cases=[]
 for supply in [.001,.01,.1,.5,1,2,5,7,10]:
  c=copy.deepcopy(reference);c['config'].update(dic_exchange=1.,nutrient_supply=supply)
  c.update(id='D14_Cbath_supply_'+str(supply),family='0D carbon-replete boundary counterpart',
     interpretation='Declared artificial DIC reservoir relaxation1/d, target2200. Check actual minimumDIC. Not air-sea carbonate chemistry.')
  cases.append(c)
 c=copy.deepcopy(reference);c['config']['dic_exchange']=.1;c.update(id='D14_Cbath_tau10',family='C-boundary relaxation sensitivity');cases.append(c)
 (ROOT/'research/carbon_boundary_registry.json').write_text(json.dumps(cases,indent=2))
 with ProcessPoolExecutor(max_workers=3) as pool:
  fs={pool.submit(execute,c):c['id'] for c in cases}
  for f in as_completed(fs):
   r=f.result();s=r['summary'];print(r['id'],{k:s[k] for k in ['runtime_seconds','phyto_CN_inventory_ratio','N_weighted_geomean_ESD_um','last_year_min_DIC','NPP_C_mmol_m2_d']},flush=True)
if __name__=='__main__':main()
