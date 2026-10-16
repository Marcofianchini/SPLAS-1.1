"""Additional refinement and paired physiological hypotheses; same frozen core."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import copy,json
from concurrent.futures import ProcessPoolExecutor,as_completed
from scripts.run_experiments import registry,execute,OUT
ROOT=Path(__file__).resolve().parents[1]

def main():
 base=registry();cases=[]
 def find(name):return copy.deepcopy(next(c for c in base if c['id']==name))
 for h in [.025,.0125]:
  c=find('D00_reference');c.update(id='D12_dt_'+str(h),family='Numerical0D timestep refinement',max_step=h);cases.append(c)
 c=find('F00_open_reference');c['config']['n_detritus']=16;c.update(id='F08_detritus_bins_16',family='Particle-bin refinement');cases.append(c)
 # Keep allometry and chemical state fixed, vary lumped energetic/organic-N priors.
 for cost in [.5,3.]:
  c=find('D00_reference');c['physiology']['cost_N']=cost;c.update(id='D13_Ncost_'+str(cost),family='N assimilation energy prior');cases.append(c)
 (ROOT/'research/supplement_registry.json').write_text(json.dumps(cases,indent=2))
 with ProcessPoolExecutor(max_workers=3) as pool:
  fs={pool.submit(execute,c):c['id'] for c in cases}
  for f in as_completed(fs):
   r=f.result();print(r['id'],r['summary'],flush=True)
if __name__=='__main__':main()
