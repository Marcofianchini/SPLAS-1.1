"""Prespecified ecological POM challenges, retaining every result.

State archives are exact float64; process observations can be reconstructed.
Atomic, fsynced writes plus SHA/ZIP checks prevent reuse of incomplete files.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dataclasses import asdict,replace
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse,hashlib,json,os,platform,time,zipfile
from datetime import datetime,timezone
import numpy as np
from splas2 import Model,Config,Physiology,simulate
from splas2.model import ACCELERATED
from scripts.run_stoichiometry import physiology

ROOT=Path(__file__).resolve().parents[1]
DOC=ROOT/'research/pom_validation_0d'
OUT=ROOT/'results/pom_validation_0d'
SOURCES=['splas2/model.py','splas2/stoichiometry.py','splas2/pom_observation.py','scripts/run_stoichiometry.py','scripts/run_pom_validation.py']
BASELINE='8b14dda650cde76a9559e0d8ae73ae7bc3beb0a1'


def source_hashes():return {s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in SOURCES}


def identifier(structure,light,n):return f'V_{structure}_L{light:g}_N{n:g}'


def registry():
    p=physiology();cases=[]
    base=Config(n_phyto=35,n_zoo=35,n_detritus=2,depth=1,self_shading=False,
        dilution=.02,inflow_DIN=2,surface_PAR=30)
    pri=replace(Physiology(),mortality=p.mortality,exud_C=p.exud_C,exud_N=p.exud_N,
        assimilation_C=p.assimilation_C,assimilation_N=p.assimilation_N,
        mortality_recycled=p.mortality_recycled,det_lysis_C=p.det_lysis_C,
        det_lysis_N=p.det_lysis_N,remin_DOC=p.remin_DOC,remin_DON=p.remin_DON)
    alternatives={'full':p,'own':replace(p,zoo_closure=1.),'flat':replace(p,cmax_exponent=0.),
                  'uniform':p,'prior':pri,'noZ':p}
    for structure,pp in alternatives.items():
        for light in [3.,30.,300.]:
            for n in [.5,2.,8.]:
                c=replace(base,surface_PAR=light,inflow_DIN=n)
                init=dict(DIN_value=1.,DIC_value=2200.,phyto_N=.7,zoo_N=.07)
                if structure=='noZ':
                    c=replace(c,n_zoo=0);init.update(DIN_value=1.07,DIC_value=2200+.07*pp.zoo_CN,zoo_N=0.)
                cases.append(dict(id=identifier(structure,light,n),family='holdout_plane',structure=structure,
                    config=asdict(c),physiology=asdict(pp),initial=init,days=3650.,max_step=.025,
                    method='imex',selectivity='uniform' if structure=='uniform' else 'selective'))
    by={c['id']:c for c in cases}
    for structure in ['full','own']:
        for n in [.5,2.,8.]:
            for mortality in [.05,.4]:
                cc=json.loads(json.dumps(by[identifier(structure,30,n)]))
                cc.update(id=cc['id']+f'_m{mortality:g}',family='closure_gradient')
                cc['physiology']['zoo_mortality']=mortality;cases.append(cc)
    for n in [.5,8.]:
        for eta in [.335,.67]:
            for rate in [.05,.15,.45]:
                if eta==.67 and rate==.15:continue
                cc=json.loads(json.dumps(by[identifier('full',30,n)]))
                cc.update(id=cc['id']+f'_a{eta:g}_r{rate:g}',family='N_retention_recycling')
                cc['physiology'].update(assimilation_N=eta,det_lysis_N=rate,remin_DON=rate);cases.append(cc)
    for structure in ['full','own']:
        for light,n in [(3.,.5),(300.,.5),(300.,8.)]:
            for resolution in ['h0125','bins70']:
                cc=json.loads(json.dumps(by[identifier(structure,light,n)]))
                cc.update(id=cc['id']+'_'+resolution,family='numerical_control')
                if resolution=='h0125':cc['max_step']=.0125
                else:cc['config'].update(n_phyto=70,n_zoo=70)
                cases.append(cc)
    assert len(cases)==88 and len({c['id'] for c in cases})==88
    return cases


def build(case,light=None):
    c=Config(**case['config'])
    if light is not None:c=replace(c,surface_PAR=float(light))
    m=Model(c,Physiology(**case['physiology']))
    if case['selectivity']=='uniform':m.phi[:]=1.
    return m


def signature(case):return hashlib.sha256(json.dumps(dict(case=case,sources=source_hashes()),sort_keys=True).encode()).hexdigest()


def atomic_json(path,data):
    tmp=path.with_suffix(path.suffix+'.partial')
    with tmp.open('w') as f:
        json.dump(data,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    tmp.replace(path)


def atomic_archive(path,**arrays):
    tmp=path.with_suffix(path.suffix+'.partial')
    with tmp.open('wb') as f:
        np.savez_compressed(f,**arrays);f.flush();os.fsync(f.fileno())
    with zipfile.ZipFile(tmp) as z:
        if z.testzip() is not None:raise IOError('Failed ZIP integrity '+str(tmp))
    tmp.replace(path)


def reusable(case,path):
    jp=path.with_suffix('.json')
    if jp.exists() and path.exists():
        old=json.loads(jp.read_text())
        if old.get('run_signature')==signature(case) and old.get('archive_sha256')==hashlib.sha256(path.read_bytes()).hexdigest():
            with zipfile.ZipFile(path) as z:
                if z.testzip() is None:return old
    return None


def record(case,m,times,states,steps,rejected,started,**extra):
    path=OUT/(case['id']+'.npz');OUT.mkdir(parents=True,exist_ok=True)
    atomic_archive(path,time=times,state=states,sizes=m.sizes,zsizes=m.zsizes,traits=m.traits)
    initial=m.budgets(states[0],True);adjusted=np.array([m.budgets(y,True) for y in states])
    import scipy,numba,llvmlite
    rec=dict(case=case,run_signature=signature(case),source_hashes=source_hashes(),
        initial_total_CN=m.budgets(states[0]).tolist(),initial_adjusted_CN=initial.tolist(),steps=steps,rejected_steps=rejected,
        runtime_seconds=time.perf_counter()-started,archive_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        maximum_relative_residual=(abs(adjusted-initial).max(axis=0)/np.maximum(abs(initial),1.)).tolist(),
        minimum_live_state=float(states[:,:m.ntr].min()),minimum_DIC=float(states[:,1].min()),
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
            numba=numba.__version__,llvmlite=llvmlite.__version__),**extra)
    atomic_json(path.with_suffix('.json'),rec)
    return rec


def execute(case):
    if not ACCELERATED:raise RuntimeError('Numba acceleration required')
    path=OUT/(case['id']+'.npz');old=reusable(case,path)
    if old:return old
    started=time.perf_counter();m=build(case)
    if case['family'] in ['sequence_memory','reserve_common_garden']:return execute_history(case,started)
    times=np.unique(np.r_[np.arange(0.,case['days'],2.),np.arange(1825.,3651.,365.),2737.5,case['days']])
    r=simulate(m,case['days'],y0=m.initial(**case['initial']),max_step=case['max_step'],times=times,method=case['method'])
    return record(case,m,r.time,r.state,r.steps,r.rejected_steps,started)


def memory_registry():
    by={c['id']:c for c in registry()};cases=[];arms={
        'control':[(3,3,0),(3,3,0),(30,3,0)],
        'Ne':[(3,3,.5),(3,3,0),(30,3,0)],
        'Nl':[(3,3,0),(3,3,.5),(30,3,0)],
        'Le':[(3,30,0),(3,3,0),(30,3,0)],
        'Ll':[(3,3,0),(3,30,0),(30,3,0)],
        'NeLl':[(3,3,.5),(3,30,0),(30,3,0)],
        'LeNl':[(3,30,0),(3,3,.5),(30,3,0)]}
    for n in [.5,2.]:
        parent=by[identifier('full',3,n)];path=OUT/(parent['id']+'.npz')
        meta=json.loads(path.with_suffix('.json').read_text())
        if not reusable(parent,path):raise ValueError('Invalid common-state parent')
        with np.load(path) as z:
            times=z['time'];state=z['state'];m=build(parent)
            bn=m.traits[:,4][None,:]*state[:,m.a,0]+state[:,m.nr,0]
            eligible=np.flatnonzero(times>=3285.)
            indices={'min':int(eligible[np.argmin(bn[eligible].sum(axis=1))]),
                     'max':int(eligible[np.argmax(bn[eligible].sum(axis=1))])}
        for phase,index in indices.items():
            for step in [.025,.0125,.00625]:
                for arm,segments in arms.items():
                    cc=json.loads(json.dumps(parent));cc.update(id=f'T_N{n:g}_{phase}_{arm}_h{step:g}',
                        family='sequence_memory',days=36.,max_step=step,segments=segments,
                        seed=dict(parent_id=parent['id'],archive_sha256=meta['archive_sha256'],index=index,time=float(times[index])),
                        arm=arm,phase=phase)
                    cases.append(cc)
    return cases


def garden_registry():
    memory=memory_registry();cases=[]
    for donor in memory:
        if donor['arm'] not in ['NeLl','LeNl']:continue
        control=next(c for c in memory if c['phase']==donor['phase'] and c['config']['inflow_DIN']==donor['config']['inflow_DIN']
            and c['max_step']==donor['max_step'] and c['arm']=='control')
        metas=[]
        for case in [donor,control]:
            path=OUT/(case['id']+'.npz')
            meta=reusable(case,path)
            if not meta:raise ValueError('Missing valid donor/control '+case['id'])
            metas.append(meta)
        cc=json.loads(json.dumps(donor));cc.update(id='G_'+donor['id'][2:],family='reserve_common_garden',days=1.,
            transplant=dict(donor=donor['id'],donor_sha256=metas[0]['archive_sha256'],
                control=control['id'],control_sha256=metas[1]['archive_sha256'],time=6.,tiny_fraction=1e-10))
        cases.append(cc)
    return cases


def exact_state(path,sha,index=None,at=None):
    if hashlib.sha256(path.read_bytes()).hexdigest()!=sha:raise ValueError('Seed hash mismatch')
    with np.load(path) as z:
        if index is None:
            hit=np.flatnonzero(abs(z['time']-at)<1e-9)
            if len(hit)!=1:raise ValueError('Ambiguous seed time')
            index=int(hit[0])
        return z['state'][index].copy()


def execute_history(case,started):
    m=build(case);extra={}
    if case['family']=='sequence_memory':
        seed=case['seed'];y=exact_state(OUT/(seed['parent_id']+'.npz'),seed['archive_sha256'],index=seed['index'])
        segments=case['segments']
    else:
        tr=case['transplant']
        donor=exact_state(OUT/(tr['donor']+'.npz'),tr['donor_sha256'],at=6.)
        y=exact_state(OUT/(tr['control']+'.npz'),tr['control_sha256'],at=6.)
        before=m.budgets(y);_,bn=m.biomass(donor)
        good=(bn[:,0]>=tr['tiny_fraction']*bn.sum())&(donor[m.a,0]>0)
        receiver_N=m.biomass(y)[1][:,0];fallback=float(receiver_N[~good].sum()/receiver_N.sum())
        for pool in [m.cr,m.nr]:
            per=np.divide(donor[pool,0],donor[m.a,0],out=np.zeros(m.nP),where=donor[m.a,0]>0)
            yy=y[pool,0].copy();yy[good]=per[good]*y[m.a,0][good];y[pool,0]=yy
        # Diagnostic counters remain those of the receiving control; new initial
        # live inventories are intentionally different and explicitly measured.
        extra=dict(transplanted_CN=(m.budgets(y)-before).tolist(),fallback_receiver_N_fraction=fallback)
        segments=[(1.,3.,0.)]
    states=[y.copy()];times=[0.];clock=0.;steps=rejected=0;events=[]
    initial_exact_sha=hashlib.sha256(y.tobytes()).hexdigest()
    for duration,light,pulse in segments:
        if pulse:
            y[0,0]+=pulse;y[m.ntr+11,0]+=pulse
            times.append(clock);states.append(y.copy());events.append(dict(time=clock,N=pulse))
        mm=build(case,light)
        grid=np.unique(np.r_[np.arange(0.,duration,.1),duration])
        r=simulate(mm,duration,y0=y,max_step=case['max_step'],times=grid,method=case['method'])
        times.extend(clock+r.time[1:]);states.extend(r.state[1:]);y=r.state[-1].copy();clock+=duration
        steps+=r.steps;rejected+=r.rejected_steps
    return record(case,m,np.asarray(times),np.asarray(states),steps,rejected,started,
        events=events,common_initial_state_sha256=initial_exact_sha,**extra)


def register(cases,name):
    DOC.mkdir(parents=True,exist_ok=True);path=DOC/(name+'.json')
    payload=dict(baseline_commit=BASELINE,registered_before_integration=True,
        protocol_sha256=hashlib.sha256((DOC/'PROTOCOL.md').read_bytes()).hexdigest(),
        source_hashes=source_hashes(),cases=cases)
    if path.exists():
        previous=json.loads(path.read_text());previous.pop('registered_at_UTC',None)
        if previous!=payload:raise ValueError('Registry changed: use an explicit new registry')
    else:
        payload['registered_at_UTC']=datetime.now(timezone.utc).isoformat();atomic_json(path,payload)
    return payload


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--stage',choices=['community','memory','garden'],default='community')
    parser.add_argument('--register-only',action='store_true');parser.add_argument('--ids',nargs='*');a=parser.parse_args()
    cases={'community':registry,'memory':memory_registry,'garden':garden_registry}[a.stage]()
    register(cases,a.stage+'_registry')
    if a.ids:
        selected=[c for c in cases if c['id'] in a.ids]
        if len(selected)!=len(a.ids):raise ValueError('Unknown or duplicate ID')
    else:selected=cases
    print(json.dumps(dict(stage=a.stage,registered=len(cases),selected=len(selected),accelerated=ACCELERATED)),flush=True)
    if a.register_only:return
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futures={pool.submit(execute,c):c['id'] for c in selected}
        for f in as_completed(futures):
            r=f.result();print(json.dumps(dict(id=r['case']['id'],seconds=round(r['runtime_seconds'],2),
                residual=r['maximum_relative_residual'],minimum=r['minimum_live_state'])),flush=True)

if __name__=='__main__':main()
