"""Run the revised idealized reference and same-state step refinement."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import copy,json,hashlib
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scripts.run_experiments import registry,execute
ROOT=Path(__file__).resolve().parents[1]

def main():
 case=copy.deepcopy(next(c for c in registry() if c['id']=='F00_open_reference'))
 case['config'].update(n_detritus=16,dic_exchange=.01)
 case.update(id='R00_buffered_reference',family='Revised idealized reference',max_step=.0125,
   interpretation='Constant historical idealized PAR; C boundary buffered with explicit input/output ledger;16 particle bins;70+70 plankton classes;40 cells5m. Research configuration, not site forcing.')
 cases=[case];path=ROOT/'research/recommended_registry.json';path.write_text(json.dumps(cases,indent=2))
 (ROOT/'research/recommended_idealized.json').write_text(json.dumps(case,indent=2))
 r=execute(case);print(r['id'],r['summary'],flush=True)
 statefile='results/thesis_replay/R00_buffered_reference.npz';y=np.load(ROOT/statefile)['final_state'];digest=hashlib.sha256(y.tobytes()).hexdigest()
 for h in [.0125,.00625]:
  c=copy.deepcopy(case);c.update(id='R01_same_state_dt_'+str(h),family='Same-state1year step refinement',days=365.,max_step=h,
    initial_state_file=statefile,initial_state_sha256=digest,
    interpretation='Continuation from the same physical10year state; counters reset; constant forcing. Compare these two runs, not the original10year averaging window.')
  cases.append(c)
 path.write_text(json.dumps(cases,indent=2))
 with ProcessPoolExecutor(max_workers=2) as pool:
  fs=[pool.submit(execute,c) for c in cases[1:]]
  for f in as_completed(fs):
   r=f.result();print(r['id'],r['summary'],flush=True)
if __name__=='__main__':main()
