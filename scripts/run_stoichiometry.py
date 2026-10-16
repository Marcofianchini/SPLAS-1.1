"""Registered, matched-inventory experiments for emergent C/N patterns.

PYTHONPATH=<numba environment> python scripts/run_stoichiometry.py --workers 4
Registry is written before any integration; no old results or solver are changed.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dataclasses import asdict,replace
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse,hashlib,json,platform,time,zipfile
import numpy as np
from splas2 import Model,Config,Physiology,simulate
from splas2.model import ACCELERATED
from splas2.stoichiometry import observe

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/stoichiometry'
SOURCES=['splas2/model.py','splas2/stoichiometry.py','scripts/run_stoichiometry.py']

def physiology():
    return replace(Physiology(**json.loads((ROOT/'research/calibrated_physiology.json').read_text())),
        mortality=.1,exud_C=.02,exud_N=.01,assimilation_C=.396,assimilation_N=.67,
        mortality_recycled=1.,det_lysis_C=.15,det_lysis_N=.15,remin_DOC=.15,remin_DON=.15)

def registry():
    p=physiology(); cases=[]
    def add(id,family,c,pp=p,days=730.,step=.05,selectivity='selective',initial=None,segments=None,note=''):
        cases.append(dict(id=id,family=family,config=asdict(c),physiology=asdict(pp),days=days,
            max_step=step,method='imex',selectivity=selectivity,
            initial=initial or dict(DIN_value=1.,phyto_N=.7,zoo_N=.07),segments=segments,note=note))
    # Keep both total C and N equal when consumers are removed.
    mono=Config(n_phyto=1,n_zoo=1,n_detritus=2,depth=1,phyto_min=5,phyto_max=5,
        zoo_min=42,zoo_max=42,self_shading=False,dilution=.05,inflow_DIN=1)
    for light in [15,45,135,405,1215]:
        for nutrient in [.25,1,4,16]:
            c=replace(mono,surface_PAR=light,inflow_DIN=nutrient)
            for eta in [.335,.67,1.]:
                add(f'M_L{light}_N{nutrient:g}_a{eta:g}','minimal_resource_plane',c,replace(p,assimilation_N=eta))
            add(f'M_L{light}_N{nutrient:g}_noZ','minimal_no_grazers',replace(c,n_zoo=0),
                initial=dict(DIN_value=1.07,DIC_value=2200+.07*p.zoo_CN,phyto_N=.7,zoo_N=0))
    community=replace(mono,n_phyto=35,n_zoo=35,phyto_min=1,phyto_max=200,zoo_min=2,zoo_max=1000,dilution=.02)
    for light in [45,135,405]:
        for nutrient in [1,4,16]:
            c=replace(community,surface_PAR=light,inflow_DIN=nutrient)
            for eta in [.335,.67]:
                add(f'C_L{light}_N{nutrient}_a{eta:g}','community_resource_plane',c,replace(p,assimilation_N=eta),days=3650.)
            add(f'C_L{light}_N{nutrient}_uniform','community_uniform_encounter',c,days=3650.,selectivity='uniform',
                note='phi=1 exactly; equal per-prey preference, same maximum ingestion per consumer, encounter availability also changes.')
            add(f'C_L{light}_N{nutrient}_noZ','community_no_grazers',replace(c,n_zoo=0),days=3650.,
                initial=dict(DIN_value=1.07,DIC_value=2200+.07*p.zoo_CN,phyto_N=.7,zoo_N=0))
    # Thesis-style ordered assembly, corrected to exactly equal total inventories.
    grid=np.geomspace(2,1000,70)
    for light in [45,405]:
        for nz in [0,1,2,4,8,16]:
            c=replace(community,n_phyto=70,n_zoo=nz,zoo_max=grid[max(nz-1,0)],dilution=0,surface_PAR=light)
            init=dict(DIN_value=1-.01*nz,DIC_value=2200-.01*nz*p.zoo_CN,phyto_N=7.,zoo_N=.01*nz)
            add(f'A_L{light}_Z{nz}','matched_closed_assembly',c,days=3650.,initial=init)
    # Preselected resolution/step controls, not selected by their eventual sign.
    for light,nutrient in [(45,1),(405,1),(405,16)]:
        c=replace(community,surface_PAR=light,inflow_DIN=nutrient)
        add(f'R_L{light}_N{nutrient}_h025','time_refinement',c,days=3650.,step=.025)
        for n in [70,140]:
            add(f'R_L{light}_N{nutrient}_bins{n}','size_refinement',replace(c,n_phyto=n,n_zoo=n),days=3650.)
    # Sequence reversal has equal light dose and N pulse, same final environment.
    pulsebase=replace(community,surface_PAR=135,inflow_DIN=4)
    histories={
        'N_then_light':[(730,135,0),(3,135,1),(3,405,0),(30,135,0)],
        'light_then_N':[(730,135,0),(3,405,0),(3,135,1),(30,135,0)],
        'N_only':[(730,135,0),(6,135,1),(30,135,0)],
        'light_only':[(730,135,0),(3,405,0),(3,135,0),(30,135,0)],
        'control':[(730,135,0),(6,135,0),(30,135,0)]}
    for name,segs in histories.items():
        add('T_'+name,'sequence_memory',pulsebase,days=766.,segments=segs,
            note='Pulse at segment start is recorded in external N ledger; light dose identical only for the two reversed sequences.')
    return cases

def source_hashes():
    return {s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in SOURCES}

def signature(case):
    return hashlib.sha256(json.dumps(dict(case=case,sources=source_hashes()),sort_keys=True).encode()).hexdigest()

def build(case):
    m=Model(Config(**case['config']),Physiology(**case['physiology']))
    if case['selectivity']=='uniform': m.phi[:]=1.
    elif case['selectivity']!='selective':raise ValueError('Unknown selection rule')
    return m

def execute(case):
    if not ACCELERATED:raise RuntimeError('Numba is required for this registered ensemble')
    OUT.mkdir(parents=True,exist_ok=True); path=OUT/(case['id']+'.npz'); jp=path.with_suffix('.json')
    sig=signature(case)
    if jp.exists() and path.exists():
        old=json.loads(jp.read_text())
        if old.get('run_signature')==sig and old.get('archive_sha256')==hashlib.sha256(path.read_bytes()).hexdigest():
            with zipfile.ZipFile(path) as z:
                if z.testzip() is None:return old
    m=build(case); y=m.initial(**case['initial']); baseline=m.budgets(y,True)
    start=time.perf_counter(); steps=rejected=0
    # Segment boundaries are preserved, including pre/post pulse snapshots.
    if not case['segments']:
        days=case['days']; window=days-min(365.,days/2)
        times=np.unique(np.r_[np.arange(0,days,2.),window,days])
        result=simulate(m,days,y0=y,max_step=case['max_step'],times=times,method=case['method'])
        states=result.state; steps=result.steps;rejected=result.rejected_steps; events=[]
    else:
        times=[0.];states=[y.copy()];t=0.;events=[]
        for duration,light,pulse in case['segments']:
            mm=Model(replace(m.config,surface_PAR=light),m.physiology)
            if pulse:
                y[0]+=pulse;y[m.ntr+11]+=pulse
                # At a pulse the same time has left/right states; never divide by dt here.
                times.append(t);states.append(y.copy());events.append(dict(time=t,N=pulse))
            dt=2. if t<730 else .1
            rr=simulate(mm,duration,y0=y,max_step=case['max_step'],save_interval=dt,method=case['method'])
            times.extend(t+rr.time[1:]);states.extend(rr.state[1:]);y=rr.state[-1].copy();t+=duration
            steps+=rr.steps;rejected+=rr.rejected_steps
        times=np.asarray(times);states=np.asarray(states)
    observations={}
    for j,(t,y) in enumerate(zip(times,states)):
        # For temporal trials preserve the actual environmental history.
        mm=m
        if case['segments']:
            bound=0.
            for duration,light,pulse in case['segments']:
                if t<bound+duration-1e-9 or bound+duration>=case['days']:
                    mm=Model(replace(m.config,surface_PAR=light),m.physiology);break
                bound+=duration
        obs=observe(mm,y,float(t))
        for key,value in obs.items():observations.setdefault(key,[]).append(value)
    arrays={k:np.asarray(v) for k,v in observations.items()}
    np.savez_compressed(path,time=times,state=states,sizes=m.sizes,zsizes=m.zsizes,traits=m.traits,**arrays)
    adj=np.array([m.budgets(z,True) for z in states]); residual=abs(adj-baseline).max(axis=0)
    import scipy,numba,llvmlite
    record=dict(case=case,run_signature=sig,source_hashes=source_hashes(),steps=steps,rejected_steps=rejected,
        runtime_seconds=time.perf_counter()-start,initial_total_CN=baseline.tolist(),
        minimum_live_state=float(states[:,:m.ntr].min()),maximum_relative_residual=(residual/np.maximum(abs(baseline),1.)).tolist(),
        minimum_DIC=float(states[:,1].min()),events=events,
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,numba=numba.__version__,llvmlite=llvmlite.__version__),
        archive_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    jp.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    return record

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--ids',nargs='*');parser.add_argument('--register-only',action='store_true');a=parser.parse_args()
    cases=registry(); destination=ROOT/'research/stoichiometry/registry.json'
    destination.write_text(json.dumps(dict(design='Resource quality -> trophic feedback -> size structure -> temporal memory',
        registered_before_integration=True,baseline_commit='a8192dc360d432b325fb181b4da70dd58b4bd440',cases=cases),indent=2)+'\n')
    if a.register_only:print(f'Registered {len(cases)} cases');return
    selected=[c for c in cases if not a.ids or c['id'] in a.ids]
    if a.ids and len(selected)!=len(a.ids):raise ValueError('Unknown or duplicate experiment ID')
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futures={pool.submit(execute,c):c['id'] for c in selected}
        for f in as_completed(futures):
            r=f.result();print(json.dumps(dict(id=r['case']['id'],seconds=round(r['runtime_seconds'],2),residual=r['maximum_relative_residual'])),flush=True)

if __name__=='__main__':main()
