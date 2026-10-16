"""Conservative spatial/particle refinement from the revised reference state."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import copy,json,time,hashlib
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from splas2 import Config,Physiology,Model
from splas2.model import _particle_weights
from scripts.run_experiments import execute
ROOT=Path(__file__).resolve().parents[1]

def main():
 path=ROOT/'research/recommended_registry.json';base=None;source=ROOT/'results/thesis_replay/R00_buffered_reference.npz'
 for _ in range(720):
  if source.exists() and source.with_suffix('.json').exists():
   try:base=copy.deepcopy(json.load(open(path))[0]);break
   except json.JSONDecodeError:pass
  time.sleep(20)
 if base is None:raise RuntimeError('Revised reference did not finish')
 cfg=Config(**base['config']);p=Physiology(**base['physiology']);old=Model(cfg,p);y=np.load(source)['final_state'].copy();y[old.ntr:]=0
 idir=ROOT/'results/recommended_initials';idir.mkdir(parents=True,exist_ok=True);cases=[]
 for kind in ['grid80','grid160','particles32']:
  vertical=kind.startswith('grid')
  conf=copy.deepcopy(base['config']);conf['n_layers']=int(kind[4:]) if vertical else 40;conf['n_detritus']=32 if kind=='particles32' else 16
  model=Model(Config(**conf),p)
  if vertical:new=np.repeat(y,conf['n_layers']//40,axis=1)
  else:
   new=np.zeros((model.nvar,40));new[:old.dc.start]=y[:old.dc.start]
   weights=_particle_weights(old.particle_sizes,model.particle_sizes)
   new[model.dc]=weights.T@y[old.dc];new[model.dn]=weights.T@y[old.dn]
  new[model.ntr:]=0
  np.testing.assert_allclose(model.budgets(new),old.budgets(y),rtol=2e-13,atol=1e-9)
  ipath=idir/(kind+'.npz');np.savez_compressed(ipath,final_state=new)
  c=copy.deepcopy(base);c.update(id='R02_'+kind,config=conf,days=365.,family='Same-state conservative grid/particle refinement',
    initial_state_file=ipath.relative_to(ROOT).as_posix(),initial_state_sha256=hashlib.sha256(new.tobytes()).hexdigest(),
    interpretation='Same10year physical inventories; piecewise-constant conservative vertical prolongation or log-particle rebinning. Compare with R01_same_state_dt_0.0125, same365day window.')
  if kind=='grid160':c['interpretation']+=' Additional nested comparison with R02_grid80 after the40-to80 test failed.'
  cases.append(c)
 (ROOT/'research/recommended_refinement_registry.json').write_text(json.dumps(cases,indent=2))
 with ProcessPoolExecutor(max_workers=2) as pool:
  fs=[pool.submit(execute,c) for c in cases]
  for f in as_completed(fs):
   r=f.result();print(r['id'],r['summary'],flush=True)
if __name__=='__main__':main()
