"""Explicit thesis experiment registry and conservative reinterpretations.
Run from repository root: python scripts/run_experiments.py --workers 4.
Every case records actual configuration and assumptions; no archival bitwise replay.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dataclasses import replace,asdict
import json,time,argparse,traceback,hashlib,zipfile,os
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from splas2 import *
from splas2.metrics import summarize
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/thesis_replay'

def registry():
 # Physiology inferred only from Ot/Tw (exploratory transfer to community).
 raw=json.loads((ROOT/'research/calibrated_physiology.json').read_text())
 p=replace(Physiology(**raw),mortality=.1,exud_C=.02,exud_N=.01,
     assimilation_C=.396,assimilation_N=.67,mortality_recycled=0.,
     det_lysis_C=.15,det_lysis_N=.15,remin_DOC=.15,remin_DON=.15)
 base=Config(nutrient_supply=1,export_bottom=True,phyto_sinking_coefficient=.01,
     detritus_sinking_reference=8,depth=100,dic_exchange=.01,surface_PAR=300)
 cases=[]
 def add(name,family,c=base,pp=p,days=3652.,step=.05,init=None,note=''):
  cases.append({'id':name,'family':family,'config':asdict(c),'physiology':asdict(pp),
    'days':days,'max_step':step,'method':'imex','initial':init or {'DIN_value':1,'phyto_N':7,'zoo_N':.7},'interpretation':note})
 add('D00_reference','DePasquale 3.3-3.4',note='Dynamic C:N counterpart. New detrital delay and respiration; HTL losses explicit; DIC exchange declared.')
 for s in [.001,.01,.1,.5,2,5,7,10]:add('D01_supply_'+str(s),'DePasquale 3.4.1-2',replace(base,nutrient_supply=s))
 for xi in [.05,.1,.15,.25]:add('D02_zmort_'+str(xi),'DePasquale 3.4.3',pp=replace(p,zoo_mortality=xi))
 add('D03_own_closure','DePasquale 3.3.2',pp=replace(p,zoo_closure=1))
 for b in [.177,.2,.354,.5,1.]:add('D04_kernel_'+str(b),'DePasquale 3.5.1',replace(base,kernel_width=b))
 for m in [.05,.2,.4]:add('D05_pmort_'+str(m),'DePasquale 3.5.2',pp=replace(p,mortality=m))
 for s in [0.,.001,.01]:add('D06_sink_'+str(s),'DePasquale 3.5.2',replace(base,phyto_sinking_coefficient=s*base.depth))
 for n in [35,140]:add('D07_size_bins_'+str(n),'DePasquale 3.5.3',replace(base,n_phyto=n,n_zoo=n),note='Same initial integrated P_N=7,Z_N=0.7 for all resolutions.')
 for eps in [.16,.5,.66]:add('D08_yield_'+str(eps),'DePasquale 3.5.4',replace(base,nutrient_supply=5),replace(p,assimilation_C=eps*(1+p.cost_growth)),note='Net C-limited gross-growth yield = eps; N assimilation held fixed.')
 for k in [1,2,4,5]:add('D09_Kz_'+str(k),'DePasquale 3.5.4',replace(base,nutrient_supply=5,grazing_half_C=k*6.625))
 add('D10_constant_mean','DePasquale 3.5.5',replace(base,nutrient_supply=.5))
 for phase in [0,91.25]:add('D10_periodic_phase_'+str(phase),'DePasquale 3.5.5',replace(base,periodic_supply=True,supply_phase=phase),note='sin^2(pi*t/365) has period365 days and mean0.5; text calling it2years is algebraically inconsistent.')
 for nz in [0,1,2,4,8,16]:
  zmax=np.geomspace(2,1000,70)[max(0,nz-1)]
  add('D11_quota_assembly_'+str(nz),'DePasquale 2.6',replace(base,n_zoo=nz,zoo_max=zmax,nutrient_supply=0,export_bottom=False,dic_exchange=0),
      replace(p,mortality_recycled=1),init={'DIN_value':1,'phyto_N':7,'zoo_N':.01*nz},note='Quota counterpart of separately executed fixed-C:N Monod control; same initial totals as historical assembly.')
 # Fianchini 200m reference: 40 REAL cells of5m, not41 with ghost endpoints.
 f=replace(base,n_layers=40,depth=200,phyto_sinking_coefficient=.0041,diffusivity_m2_s=1e-4,
    nutrient_supply=.4,bottom_N_flux=20,surface_PAR=959.83*4.57,dic_exchange=0)
 fp=replace(p,mortality=.05,mortality_recycled=.5)
 ini={'DIN_value':1,'phyto_N':.7,'zoo_N':.07,'gradient':True}
 closed=replace(f,diffusivity_m2_s=0,export_bottom=False,nutrient_supply=0,bottom_N_flux=0)
 cp=replace(fp,mortality_recycled=1,det_lysis_C=.5,det_lysis_N=.5,remin_DOC=.5,remin_DON=.5)
 add('F00_open_reference','Fianchini 3',f,fp,init=ini,note='Source0.1S per layer + bottom4 mmolN/m3/d over5m =20 mmolN/m2/d. PAR=4.57*959.83 assuming statedW are PAR. Idealized high light; no site attribution.')
 add('F00_closed_reference','Fianchini 3',closed,cp,init=ini,note='Isolated layers, no external inputs/outputs. Dynamic stoichiometry and explicit C costs.')
 for q in [.7,.8,.9,1.1,1.2,1.3]:
  add('F01_open_light_'+str(q),'Fianchini 4.2',replace(f,surface_PAR=f.surface_PAR*q),fp,init=ini)
  add('F01_closed_light_'+str(q),'Fianchini 4.2',replace(closed,surface_PAR=f.surface_PAR*q),cp,init=ini)
 for force,label in [(1,'semiannual'),(2,'annual'),(3,'constant_matched')]:
  add('F02_'+label,'Fianchini 4.3',replace(f,forcing=force),fp,init=ini,note='All three annual PAR integrals are equal; historical sin^2(2*pi*t/365) has two light maxima/year.')
 for q in [.7,.8,.9,1.1,1.2,1.3]:add('F02_seasonal_light_'+str(q),'Fianchini 4.3',replace(f,forcing=1,surface_PAR=f.surface_PAR*q),fp,init=ini)
 for q in [.7,.8,.9,1.,1.1,1.2,1.3]:add('F03_100m_light_'+str(q),'Fianchini 5',replace(f,depth=100,n_layers=20,forcing=1,surface_PAR=f.surface_PAR*q),fp,init=ini,note='Column depth100m is a boundary, not diagnosed MLD. Integrated critical-depth test is separate.')
 for mult in [.1,10]:add('F04_mixing_'+str(mult),'DePasquale 4.6 / Fianchini mechanism',replace(f,forcing=1,diffusivity_m2_s=1e-4*mult),fp,init=ini)
 add('F04_variable_mixing','DePasquale 4.6',replace(f,forcing=1,deep_diffusivity_m2_s=1e-5),fp,init=ini,note='Conservative face diffusivity, surface1e-4 to deep1e-5m2/s.')
 add('F05_no_grazers','Fianchini 5 mechanism',replace(f,n_zoo=0,forcing=1),fp,init=ini)
 add('F05_no_shading','DePasquale 4.5.1',replace(f,self_shading=False,forcing=1),fp,init=ini)
 # Distinguish historical irradiance ambiguity from process response.
 add('F06_PAR_fraction_043','Fianchini irradiance units sensitivity',replace(f,surface_PAR=f.surface_PAR*.43),fp,init=ini,note='Alternative: reportedW were total shortwave;43percent assignedPAR. Conversion assumption, not measured.')
 add('F06_DIC_exchange','C boundary sensitivity',replace(f,dic_exchange=.01),fp,init=ini)
 # Numerical/particle sensitivities on full10year reference, same physical initial N(z).
 for n in [20,80]:
  dz=f.depth/n;ii={**ini,'gradient':False}
  add('F07_layers_'+str(n),'Numerical vertical refinement',replace(f,n_layers=n),fp,init=ii,note='Initial DIN(z)=0.5+z/5; fixed physical profile, not index-dependent across grids.')
 # Use same center-based physical DIN on reference too (equals k+1 atdz5).
 for d in [4,8]:add('F08_detritus_bins_'+str(d),'New particle C:N resolution',replace(f,n_detritus=d),fp,init=ini)
 for dt in [.025,.0125]:add('F09_dt_'+str(dt),'Numerical temporal refinement',f,fp,step=dt,init=ini)
 for formulation in [1]:add('F10_product_colimitation','Structural sensitivity',f,replace(fp,colimitation=formulation),init=ini)
 add('F10_prior_physiology','Physiological transfer sensitivity',f,replace(Physiology(),mortality=.05,assimilation_C=.396,assimilation_N=.67,mortality_recycled=.5,det_lysis_C=.15,det_lysis_N=.15,remin_DOC=.15,remin_DON=.15),init=ini)
 add('F10_selective_remin','Detrital C:N hypothesis',f,replace(fp,det_lysis_N=.3,remin_DON=.3),init=ini,note='N lysis/remin twice C: deliberately alternative, not an empirically calibrated rate.')
 return cases

def execute(case):
 OUT.mkdir(parents=True,exist_ok=True)
 path=OUT/(case['id']+'.npz');jp=OUT/(case['id']+'.json')
 signature=hashlib.sha256((json.dumps(case,sort_keys=True)+(ROOT/'splas2/model.py').read_text()).encode()).hexdigest()
 if jp.exists() and path.exists():
  try:
   cached=json.loads(jp.read_text())
   if cached.get('run_signature')==signature:
    with zipfile.ZipFile(path) as archive:
     if archive.testzip() is None:return cached
  except (OSError,ValueError,zipfile.BadZipFile):pass
 m=Model(Config(**case['config']),Physiology(**case['physiology']))
 if case.get('initial_state_file'):
  y=np.load(ROOT/case['initial_state_file'])['final_state'].copy()
  digest=hashlib.sha256(y.tobytes()).hexdigest()
  if digest!=case['initial_state_sha256']:raise ValueError('Initial-state checksum mismatch')
  y[m.ntr:]=0.
 else:y=m.initial(**case['initial'])
 if case['id'].startswith('F07_'):y[0]=.5+m.depth/5
 # Monthly early record, every2days in final2years resolves seasonal means/phase.
 days=case['days'];times=np.unique(np.r_[np.arange(0,max(0,days-730),30),np.arange(max(0,days-730),days,2),days])
 r=simulate(m,days,y0=y,max_step=case['max_step'],times=times,method=case['method'])
 summary,arr=summarize(r)
 rec={**case,'run_signature':signature,'summary':summary}
 temporary=OUT/(case['id']+'.'+str(os.getpid())+'.tmp.npz')
 np.savez_compressed(temporary,**arr);temporary.replace(path)
 tempjson=OUT/(case['id']+'.'+str(os.getpid())+'.tmp.json')
 tempjson.write_text(json.dumps(rec,indent=2,allow_nan=False));tempjson.replace(jp)
 return rec

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=4);ap.add_argument('--filter',default='');args=ap.parse_args()
 OUT.mkdir(parents=True,exist_ok=True);cases=registry();(ROOT/'research/experiment_registry.json').write_text(json.dumps(cases,indent=2))
 filters=args.filter.split(',')
 todo=[c for c in cases if any(part in c['id'] for part in filters)];start=time.perf_counter();errors=[]
 with ProcessPoolExecutor(max_workers=args.workers) as pool:
  fs={pool.submit(execute,c):c['id'] for c in todo}
  for f in as_completed(fs):
   try:
    r=f.result();s=r['summary'];print(json.dumps({'id':r['id'],'runtime':s['runtime_seconds'],'N':s['phyto_N_mmol_m2'],'CN':s['phyto_CN_inventory_ratio'],'C_error':s['max_relative_mass_error_C'],'N_error':s['max_relative_mass_error_N']}),flush=True)
   except Exception as exc:
    err={'id':fs[f],'error':repr(exc)};errors.append(err);print(json.dumps(err),flush=True);traceback.print_exc()
 status='run_status'+('_'+''.join(ch if ch.isalnum() else '_' for ch in args.filter) if args.filter else '')+'.json'
 (OUT/status).write_text(json.dumps({'cases_requested':len(todo),'errors':errors,'wall_seconds':time.perf_counter()-start},indent=2))
 if errors:sys.exit(1)
if __name__=='__main__':main()
