"""DePasquale chapter5 / light-profile interventions from one common state.

Annual perturbations after the revised10year spin-up; no claim of new
asymptotic equilibria at every setting.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import copy,json,hashlib,time,sys
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scripts.run_experiments import execute
ROOT=Path(__file__).resolve().parents[1]

def main():
 source=ROOT/'results/thesis_replay/R00_buffered_reference.npz';base=None
 for _ in range(360):
  if source.exists() and source.with_suffix('.json').exists():
   try:base=copy.deepcopy(json.load(open(ROOT/'research/recommended_registry.json'))[0]);break
   except json.JSONDecodeError:pass
  time.sleep(20)
 if base is None:raise RuntimeError('Revised reference did not finish')
 digest=hashlib.sha256(np.load(source)['final_state'].tobytes()).hexdigest();cases=[]
 def add(name,cfg=None,phys=None,family='DePasquale5 controlled export sensitivity'):
  c=copy.deepcopy(base);c['config'].update(cfg or {});c['physiology'].update(phys or {})
  c.update(id=name,days=365.,family=family,initial_state_file=source.relative_to(ROOT).as_posix(),initial_state_sha256=digest,
    interpretation='Annual intervention from the same revised10year state; compare R01_same_state_dt_0.0125. Not a separately equilibrated historical10year run.')
  cases.append(c)
 for mult in [.25,2.]:add('R04_Ninput_'+str(mult),{'nutrient_supply':.4*mult,'bottom_N_flux':20*mult})
 add('R04_remin_twice',phys={'det_lysis_C':.3,'det_lysis_N':.3,'remin_DOC':.3,'remin_DON':.3})
 add('R04_phyto_sink_twice',{'phyto_sinking_coefficient':.0082})
 add('R04_lower_C_yield',phys={'assimilation_C':.264})
 for k in [.015,.06]:add('R05_water_attenuation_'+str(k),phys={'water_attenuation':k},family='DePasquale4.5 light-profile sensitivity')
 (ROOT/'research/export_intervention_registry.json').write_text(json.dumps(cases,indent=2))
 with ProcessPoolExecutor(max_workers=3) as pool:
  fs=[pool.submit(execute,c) for c in cases]
  for f in as_completed(fs):
   r=f.result();print(r['id'],r['summary'],flush=True)
if __name__=='__main__':main()
