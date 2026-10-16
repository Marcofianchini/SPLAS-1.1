"""Explicitly labelled follow-up controls, designed after the first plane results."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import argparse,hashlib,json
from dataclasses import asdict,replace
from concurrent.futures import ProcessPoolExecutor,as_completed
from splas2 import Physiology
from scripts.run_stoichiometry import ROOT,registry as initial_registry,execute,physiology


def registry():
    original={c['id']:c for c in initial_registry()};cases=[]
    def add(id,source,family,**changes):
        c={**original[source],**changes,'id':id,'family':family,
            'control_design':'Follow-up: not part of the initial 142-case registration',
            'design_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        cases.append(c)
    p=physiology();prior=replace(Physiology(),mortality=p.mortality,exud_C=p.exud_C,exud_N=p.exud_N,
        assimilation_C=p.assimilation_C,mortality_recycled=p.mortality_recycled,
        det_lysis_C=p.det_lysis_C,det_lysis_N=p.det_lysis_N,remin_DOC=p.remin_DOC,remin_DON=p.remin_DON)
    for l,n in [(45,1),(405,1),(405,16)]:
        source=f'C_L{l}_N{n}_a0.67'
        add(f'Q_L{l}_N{n}_h0125',source,'additional_step_control',max_step=.0125)
        add(f'Q_L{l}_N{n}_flat_storage',source,'storage_allometry_control',physiology=asdict(replace(p,cmax_exponent=0)),
            note='Only cmax exponent is flattened; same capacity at5µm, no size-growing potential C storage.')
        for eta in [.335,.67]:
            add(f'Q_L{l}_N{n}_prior_a{eta:g}',source,'prior_physiology_control',
                physiology=asdict(replace(prior,assimilation_N=eta)),
                note='Unfitted shared physiology; other community rates unchanged. Several physiological rules change together.')
    for l in [15,405]:
        for h in [.025,.0125]:
            add(f'Q_M_L{l}_h{h:g}',f'M_L{l}_N1_a0.67','minimal_step_control',max_step=h,
                note='First-order IMEX splitting can separate saved-state diet C:N from stage-integrated ingestion C:N.')
    return cases


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=2);a=parser.parse_args()
    cases=registry();(ROOT/'research/stoichiometry/control_registry.json').write_text(json.dumps(dict(
        reason='Separate prey-size storage hypothesis, transferred physiological fit, and first-order split-step effects.',
        timing='Designed after the initial minimal plane and first community results; exploratory follow-up, not preregistered confirmation.',
        cases=cases),indent=2)+'\n')
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futures=[pool.submit(execute,c) for c in cases]
        for f in as_completed(futures):
            r=f.result();print(json.dumps(dict(id=r['case']['id'],seconds=round(r['runtime_seconds'],2),residual=r['maximum_relative_residual'])),flush=True)

if __name__=='__main__':main()
