"""Whole-plane refinement and a biological rare-reintroduction control."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib,json,time
from datetime import datetime,timezone
import numpy as np
from splas2 import simulate
from scripts.run_pom_validation import (registry,identifier,execute,record,build,exact_state,reusable,
    atomic_json,DOC,OUT,source_hashes,BASELINE)


def cases():
    by={c['id']:c for c in registry()};cc=[]
    protocol=hashlib.sha256((DOC/'REFINEMENT_PROTOCOL.md').read_bytes()).hexdigest()
    for light in [3.,30.,300.]:
        for n in [.5,2.,8.]:
            c=json.loads(json.dumps(by[identifier('full',light,n)]))
            c.update(id=c['id']+'_h00625',family='fine_plane',max_step=.00625,refinement_protocol_sha256=protocol);cc.append(c)
    for light,n in [(3.,.5),(300.,.5),(300.,8.)]:
        c=json.loads(json.dumps(by[identifier('full',light,n)]))
        c.update(id=c['id']+'_h003125',family='fine_corner',max_step=.003125,refinement_protocol_sha256=protocol);cc.append(c)
    for parent in [identifier('full',300,.5),identifier('full',300,.5)+'_h0125']:
        meta=json.loads((OUT/(parent+'.json')).read_text())
        for fraction in [0.,1e-5]:
            c=json.loads(json.dumps(by[identifier('full',300,.5)]))
            c.update(id='B_'+parent[2:]+('_seeded' if fraction else '_resident'),family='rare_reintroduction',days=730.,max_step=.00625,
                seed=dict(parent_id=parent,archive_sha256=meta['archive_sha256'],time=3650.),
                redistribute_PN_fraction=fraction,refinement_protocol_sha256=protocol,
                refinement_runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());cc.append(c)
    return cc


def execute_extra(c):
    if c['family']!='rare_reintroduction':return execute(c)
    old=reusable(c,OUT/(c['id']+'.npz'))
    if old:return old
    started=time.perf_counter();m=build(c);seed=c['seed']
    y=exact_state(OUT/(seed['parent_id']+'.npz'),seed['archive_sha256'],at=seed['time'])
    original=y.copy();before=m.budgets(y);fraction=c['redistribute_PN_fraction']
    if fraction:
        pn=float(m.biomass(y)[1].sum());new=m.initial(phyto_N=pn*fraction,zoo_N=0.)
        for pool in [m.a,m.cr,m.nr]:y[pool]=(1-fraction)*y[pool]+new[pool]
        deltaC=m.budgets(y)[0]-before[0];y[1,0]-=deltaC
    np.testing.assert_allclose(m.budgets(y),before,rtol=1e-14,atol=1e-13)
    r=simulate(m,c['days'],y0=y,max_step=c['max_step'],save_interval=1.,method='imex')
    return record(c,m,r.time,r.state,r.steps,r.rejected_steps,started,
        initial_parent_live_CN=before.tolist(),DIC_redistribution=float(y[1,0]-original[1,0]),
        N_redistributed=float(fraction*m.biomass(original)[1].sum()))


def main():
    cc=cases();path=DOC/'refinement_registry.json'
    reg=dict(baseline_commit=BASELINE,registered_before_integration=True,source_hashes=source_hashes(),cases=cc,
        refinement_protocol_sha256=hashlib.sha256((DOC/'REFINEMENT_PROTOCOL.md').read_bytes()).hexdigest())
    if path.exists():
        old=json.loads(path.read_text());old.pop('registered_at_UTC',None)
        if old!=reg:raise ValueError('Changed refinement registry')
    else:
        reg['registered_at_UTC']=datetime.now(timezone.utc).isoformat();atomic_json(path,reg)
    print('Registered 16 additional controls before their integration',flush=True)
    with ProcessPoolExecutor(max_workers=4) as pool:
        ff={pool.submit(execute_extra,c):c['id'] for c in cc}
        for f in as_completed(ff):
            r=f.result();print(json.dumps(dict(id=r['case']['id'],seconds=r['runtime_seconds'],residual=r['maximum_relative_residual'])),flush=True)

if __name__=='__main__':main()
