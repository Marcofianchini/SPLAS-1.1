"""Final conservative redistribution contrast at both trophic levels."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime,timezone
import hashlib,json,time
import numpy as np
from splas2 import simulate
from scripts.run_pom_validation import (registry,identifier,record,build,exact_state,reusable,
    atomic_json,DOC,OUT,source_hashes,BASELINE)


def cases():
    base=next(c for c in registry() if c['id']==identifier('full',300,.5));cc=[]
    for parent in [base['id'],base['id']+'_h0125']:
        meta=json.loads((OUT/(parent+'.json')).read_text())
        for pf,zf,label in [(0.,1e-5,'Zseeded'),(1e-5,1e-5,'PZseeded'),(1e-7,1e-7,'PZseeded_small')]:
            c=json.loads(json.dumps(base));c.update(id='B_'+parent[2:]+'_'+label,family='rare_reintroduction',days=730.,max_step=.00625,
                seed=dict(parent_id=parent,archive_sha256=meta['archive_sha256'],time=3650.),
                redistribute_PN_fraction=pf,redistribute_ZN_fraction=zf,
                rare_pool_protocol_sha256=hashlib.sha256((DOC/'RARE_POOL_PROTOCOL.md').read_bytes()).hexdigest(),
                rare_pool_runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());cc.append(c)
    return cc


def execute_extra(c):
    old=reusable(c,OUT/(c['id']+'.npz'))
    if old:return old
    start=time.perf_counter();m=build(c);seed=c['seed']
    y=exact_state(OUT/(seed['parent_id']+'.npz'),seed['archive_sha256'],at=seed['time']);original=y.copy()
    before=m.budgets(y);fraction=c['redistribute_PN_fraction'];pn=float(m.biomass(y)[1].sum())
    if fraction:
        new=m.initial(phyto_N=pn*fraction,zoo_N=0.)
        for pool in [m.a,m.cr,m.nr]:y[pool]=(1-fraction)*y[pool]+new[pool]
        y[1,0]-=m.budgets(y)[0]-before[0]
    zf=c['redistribute_ZN_fraction'];zC=float(y[m.z,0].sum())
    y[m.z,0]=(1-zf)*y[m.z,0]+zf*zC/m.nZ
    np.testing.assert_allclose(m.budgets(y),before,rtol=1e-14,atol=1e-13)
    r=simulate(m,730.,y0=y,max_step=c['max_step'],save_interval=1.,method='imex')
    return record(c,m,r.time,r.state,r.steps,r.rejected_steps,start,
        initial_parent_live_CN=before.tolist(),DIC_redistribution=float(y[1,0]-original[1,0]),
        N_redistributed=pn*fraction,N_redistributed_Z=zf*zC/m.physiology.zoo_CN)


def main():
    cc=cases();path=DOC/'rare_pool_registry.json'
    reg=dict(baseline_commit=BASELINE,registered_before_integration=True,source_hashes=source_hashes(),cases=cc,
        rare_pool_protocol_sha256=hashlib.sha256((DOC/'RARE_POOL_PROTOCOL.md').read_bytes()).hexdigest())
    if path.exists():
        old=json.loads(path.read_text());old.pop('registered_at_UTC',None)
        if old!=reg:raise ValueError('Changed rare-pool registry')
    else:reg['registered_at_UTC']=datetime.now(timezone.utc).isoformat();atomic_json(path,reg)
    print('Registered 6 final controls before their integration',flush=True)
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures={pool.submit(execute_extra,c):c['id'] for c in cc}
        for f in as_completed(futures):
            r=f.result();print(json.dumps(dict(id=r['case']['id'],seconds=r['runtime_seconds'],residual=r['maximum_relative_residual'])),flush=True)

if __name__=='__main__':main()
