"""R* oracle, predator assembly, culture phase uncertainty, and optical/mixing separation."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dataclasses import replace,asdict
import json,time
import numpy as np
from scipy.optimize import brentq
from splas2 import *
from splas2.model import INPUT_N,UPTAKE_N,GPP,RESP_P
from splas2.npz_fixed import run_fixed
from splas2.metrics import summarize
from scripts.calibrate_cultures import model_prediction,DATA,PRIOR
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/mechanisms'

def main():
 OUT.mkdir(parents=True,exist_ok=True)
 p=Physiology(**json.loads((ROOT/'research/calibrated_physiology.json').read_text()))
 out={'fixed_CN_assembly':[],'culture_phase':[],'physiology_grid':[],'mixing_experiments':[]}
 for nz in [0,1,2,4,8,16,32,70]:
  z=run_fixed(nz);y=z['final'];y0=z['initial'];n=70
  np.savez_compressed(OUT/('fixed_assembly_'+str(nz)),**{k:v for k,v in z.items() if isinstance(v,np.ndarray)})
  d={k:v for k,v in z.items() if not isinstance(v,np.ndarray)}
  d.update(n_zoo=nz,final_DIN=float(y[0]),winning_phyto_ESD=float(z['sizes'][np.argmax(y[2:72])]),final_P_N=float(y[2:72].sum()),final_Z_N=float(y[72:].sum()))
  if nz==0:
   nt=y0[0]+y0[2:].sum();ct=y0[1]+6.625*y0[2:].sum()
   def f(din):
    dic=ct-6.625*(nt-din);return 2.6*din/(.1+din)*dic/(10+dic)-.1001
   root=brentq(f,1e-8,1)
   d['analytic_Rstar']=root;d['Rstar_absolute_error']=float(abs(y[0]-root))
  out['fixed_CN_assembly'].append(d);print('fixed',nz,d['final_DIN'],flush=True)
 for key in DATA:
  for phase in [.0,.25,.5,.75]:
   pr,r=model_prediction(key,p,max_step=.01,phase=phase)
   actual=np.array([o['CN_molar'] for o in DATA[key]['observations']]);cn=pr[:,0]/pr[:,1]*14/12
   out['culture_phase'].append({'species':key,'light_phase_days':phase,'CN_RMSLE':float(np.sqrt(np.mean(np.log(cn[1:]/actual[1:])**2))),'final_CN':float(cn[-1]),'max_relative_mass_error_N':r.diagnostics()['max_relative_mass_error_N']})
 for light in [0,10,85,300,1000]:
  for din in [0,.02,.25,2,20]:
   m=Model(Config(n_phyto=1,n_zoo=0,phyto_min=10,phyto_max=10,surface_PAR=light,depth=1),p)
   r=simulate(m,30,y0=m.initial(DIN_value=din,phyto_N=.1,zoo_N=0),save_interval=.5)
   s,a=summarize(r,window=10);out['physiology_grid'].append({'PAR':light,'initial_DIN':din,**s})
   if din in [0,2] and light in [0,85,1000]:np.savez_compressed(OUT/f'physiology_L{light}_N{din}',**a)
 # Distinct prescribed turbulent mixing depths. Short mechanistic tests, not
 # asymptotic thesis replacements. Uniform DIN removes index/grid confounding.
 fp=replace(p,mortality=.05,assimilation_C=.396,assimilation_N=.67,mortality_recycled=.5,
    det_lysis_C=.15,det_lysis_N=.15,remin_DOC=.15,remin_DON=.15)
 for mix_depth in [25.,50.,100.]:
  for light in [40.,120.,360.]:
   for grazers in [False,True]:
    cfg=Config(n_phyto=35,n_zoo=35 if grazers else 0,n_layers=20,depth=100,surface_PAR=light,
      forcing=2,export_bottom=True,nutrient_supply=0,bottom_N_flux=2,dic_exchange=0)
    m=Model(cfg,fp);faces=np.arange(1,20)*m.dz
    m.Dface=np.where(faces<=mix_depth,1e-2,1e-6)*86400
    r=simulate(m,365,y0=m.initial(DIN_value=5,phyto_N=.07,zoo_N=.007),max_step=.05,save_interval=2,method='imex')
    ss,ar=summarize(r,window=180);first=r.state[0];last=r.state[-1]
    net=np.diff(np.log(np.maximum(ar['phyto_C'],1e-200)))/np.diff(r.time)
    # Freeze instantaneous state to diagnose BIOLOGICAL net-C balance.
    # Zero sinking/diffusion: transport is not falsely interpreted as growth.
    mb=Model(replace(cfg,export_bottom=False,diffusivity_m2_s=0,bottom_N_flux=0),fp)
    thresholds=[];euphotic=[];gamma_at_M=[]
    for t,y in zip(r.time,r.state):
     rhs=mb.rhs(t,y);netc=(mb.traits[:,3,None]*rhs[mb.a]+rhs[mb.cr]+fp.reserve_N_CN*rhs[mb.nr]).sum(axis=0)
     integral=np.cumsum(netc)*m.dz
     # First + to - crossing; if absent, explicit status instead of false scalar.
     root=None
     for j in range(1,len(integral)):
      if integral[j-1]>0 and integral[j]<=0:
       root=float(np.interp(0,[-integral[j-1],-integral[j]],m.depth[j-1:j+1]));break
     thresholds.append(root)
     irradiance,chl=mb.light(t,y);surface=cfg.surface_PAR # ratio threshold uses true instantaneous surface
     from splas2.model import surface_light
     surface=surface_light(t,cfg.surface_PAR,cfg.forcing)
     ids=np.where(irradiance>=.01*surface)[0] if surface>0 else []
     euphotic.append(float(m.depth[ids[-1]]) if len(ids) else None)
     n=int(mix_depth/m.dz);gamma_at_M.append(float(netc[:n].sum()*m.dz))
    case=f'mixing_H{int(mix_depth)}_L{int(light)}_Z{int(grazers)}'
    np.savez_compressed(OUT/case,**ar,net_growth_C=net,critical_depth=np.array([v if v is not None else np.nan for v in thresholds]),
      Zeu=np.array([v if v is not None else np.nan for v in euphotic]),integrated_bio_growth_at_MLD=np.array(gamma_at_M),Dface_m2_d=m.Dface)
    out['mixing_experiments'].append({'id':case,'mixing_depth_m_prescribed':mix_depth,'surface_PAR':light,'grazers':grazers,'days':365,
       'D_mixed_m2_s':1e-2,'D_below_m2_s':1e-6,'peak_log_C_growth_d':float(net.max()),
       'fraction_days_positive_integrated_bio_growth_at_MLD':float((np.array(gamma_at_M)>0).mean()),
       'fraction_defined_critical_depth':float(np.mean([v is not None for v in thresholds])),**ss})
    print(case,'done',flush=True)
 (OUT/'mechanisms.json').write_text(json.dumps(out,indent=2,allow_nan=False))
if __name__=='__main__':main()
