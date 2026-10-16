"""Paired controls prompted by first-order step sensitivity, separately registered."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import argparse,hashlib,json
from concurrent.futures import ProcessPoolExecutor,as_completed
from scripts.run_stoichiometry import ROOT,registry as initial_registry,execute


def registry():
    base={c['id']:c for c in initial_registry()};cases=[]
    for l,n in [(45,1),(405,1),(405,16)]:
        for name,h,bins in [('h0125',.0125,35),('bins70',.05,70)]:
            source=base[f'C_L{l}_N{n}_a0.335']
            cases.append({**source,'id':f'B_L{l}_N{n}_a0335_{name}','family':'paired_eta_robustness',
                'max_step':h,'config':{**source['config'],'n_phyto':bins,'n_zoo':bins},
                'note':'Paired lower-N-assimilation case for the corresponding refined reference.',
                'design_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    for name in ['N_then_light','light_then_N']:
        source=base['T_'+name]
        cases.append({**source,'id':'B_T_'+name,'family':'sequence_step_control','max_step':.0125,
            'design_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    return cases


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=3);a=parser.parse_args()
    cases=registry();(ROOT/'research/stoichiometry/robustness_registry.json').write_text(json.dumps(dict(
        timing='Second exploratory control block, after observing step-sensitive community composition.',
        reason='Test the sign of the actual paired ecological contrast, not just stability of an unpaired reference.',cases=cases),indent=2)+'\n')
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(execute,c) for c in cases]):
            r=f.result();print(json.dumps(dict(id=r['case']['id'],seconds=round(r['runtime_seconds'],2),residual=r['maximum_relative_residual'])),flush=True)

if __name__=='__main__':main()
