"""Shared parameters on Ot/Tw; predict Tp/Ms without refitting; no field validation."""
from pathlib import Path
import json,time
from dataclasses import replace,asdict
import numpy as np
from scipy.optimize import least_squares
from splas2 import *
ROOT=Path(__file__).resolve().parents[1]
DATA=json.loads((ROOT/'data/processed/cultures.json').read_text())['species']
PARAMS=['mu_scale','alpha_Chl','nmin_fraction','nmin_exponent','cmax_ratio','cmax_exponent']
BOUNDS=([.3,.1,.08,-.6,1.2,0.],[4.,4.,.40,0.,5.,.8])
# Costs are declared priors; lab death rate deliberately separate from grazing community.
PRIOR=Physiology(mortality=.01,exud_C=.01,exud_N=.005,det_lysis_C=0,det_lysis_N=0,remin_DOC=0,remin_DON=0)

def model_prediction(key,p,max_step=.02,phase=0):
 d=DATA[key];obs=d['observations'];o=obs[0]
 c=Config(n_phyto=1,n_zoo=0,n_detritus=1,depth=1,surface_PAR=85,forcing=4,self_shading=False,
          dic_exchange=.1,inflow_DIC=4000,light_phase=phase)
 m=Model(c,replace(p,temperature=d['temperature_C']),volumes=[d['volume_um3']],
         reference_C=[d['reference_C_pg']/120],reference_N=[d['reference_N_pg']/140])
 y=m.initial(phyto_N=0,zoo_N=0,DIN_value=o['DIN_uM'],DIC_value=4000)
 a=o['cells_ml']*.01  # cells/mL * 1e6 mL/m3 / 1e8 cells
 y[m.a]=a;y[m.nr]=a*(o['N_pg']/140-m.traits[0,4]);y[m.cr]=a*(o['C_pg']/120-m.traits[0,3])-p.reserve_N_CN*y[m.nr]
 if y[m.cr,0]<0 or y[m.nr,0]<0:raise ValueError('Initial reserves below structural requirement')
 times=np.array([z['t'] for z in obs]);r=simulate(m,float(times[-1]),y0=y,times=times,max_step=max_step)
 aa=r.state[:,m.a,0].ravel();bc,bn=[],[]
 for s in r.state:
  cc,nn=m.biomass(s);bc.append(cc[0,0]);bn.append(nn[0,0])
 pred=np.column_stack((np.array(bc)/aa*120,np.array(bn)/aa*140,aa*100,r.state[:,0,0]))
 return pred,r

def residual(v):
 p=replace(PRIOR,**dict(zip(PARAMS,v)));res=[]
 for key in ['Ot','Tw']:
  try:pr,r=model_prediction(key,p)
  except (ValueError,RuntimeError,FloatingPointError):return np.full(31,1e3)
  obs=DATA[key]['observations']
  for j,name in enumerate(['C_pg','N_pg','cells_ml']):
   measured=np.array([o[name] for o in obs]);res.extend(np.log(np.maximum(pr[1:,j],1e-30)/measured[1:]))
  # DIN uses asinh with 1uM scale; zeros are valid measured zeros.
  for i,o in enumerate(obs[1:],1):
   if o['DIN_uM'] is not None:res.append(.3*(np.arcsinh(pr[i,3])-np.arcsinh(o['DIN_uM'])))
 return np.array(res)

def main():
 t=time.perf_counter();outdir=ROOT/'results/cultures';outdir.mkdir(parents=True,exist_ok=True)
 starts=[[1,.5,.25,-.15,2.2,.2],[2,1,.3,-.25,3,.35],[.7,.3,.18,-.1,1.8,.15]]
 length=len(residual(np.array(starts[0])));print('residuals',length,flush=True)
 solutions=[]
 for seed,x in enumerate(starts):
  fit=least_squares(residual,x,bounds=BOUNDS,loss='soft_l1',f_scale=.2,max_nfev=180,diff_step=2e-3,ftol=1e-7,xtol=1e-6,gtol=1e-6)
  rec={'start':seed,'parameters':dict(zip(PARAMS,fit.x.tolist())),'cost':float(fit.cost),'nfev':fit.nfev,'success':bool(fit.success),'message':fit.message}
  solutions.append((fit,rec));print(json.dumps(rec),flush=True)
 best=min(solutions,key=lambda x:x[0].cost)[0];p=replace(PRIOR,**dict(zip(PARAMS,best.x)))
 sv=np.linalg.svd(best.jac,compute_uv=False); metrics=[];predictions={}
 for key,d in DATA.items():
  pr,r=model_prediction(key,p,max_step=.005);obs=d['observations'];r.save(outdir/key)
  actual=np.array([[o[k] if o[k] is not None else np.nan for k in ['C_pg','N_pg','cells_ml','DIN_uM']] for o in obs])
  item={'species':key,'split':d['split'],**r.diagnostics()}
  for j,name in enumerate(['C_pg','N_pg','cells_ml']):
   item[name+'_RMSLE']=float(np.sqrt(np.mean(np.log(pr[1:,j]/actual[1:,j])**2)))
   item[name+'_persistence_RMSLE']=float(np.sqrt(np.mean(np.log(actual[0,j]/actual[1:,j])**2)))
  cn=pr[:,0]/pr[:,1]*14/12;cnobs=np.array([o['CN_molar'] for o in obs])
  # Report ratios of mean quotas and mean replicate ratios separately.
  item['CN_RMSLE']=float(np.sqrt(np.mean(np.log(cn[1:]/cnobs[1:])**2)))
  item['CN_Redfield_RMSLE']=float(np.sqrt(np.mean(np.log(6.625/cnobs[1:])**2)))
  item['CN_initial_ratio_RMSLE']=float(np.sqrt(np.mean(np.log(cn[0]/cnobs[1:])**2)))
  item['final_CN_predicted']=float(cn[-1]);item['final_CN_observed']=float(cnobs[-1]);metrics.append(item)
  predictions[key]={'times':[o['day'] for o in obs],'observed':actual.tolist(),'predicted':pr.tolist(),'CN_predicted':cn.tolist(),'CN_observed':cnobs.tolist()}
  print(json.dumps(item),flush=True)
 report={'version':'0.2.0','split':'Ot and Tw calibration; Tp and Ms held-out; conditional on observed initial values and day0 reference quotas',
    'development_disclosure':'Initial kinetics showed weak identifiability and an incorrect N factor in the low-light slope. Corrected algebraically and reduced fitted parameters using the training failure. Initial Tp/Ms scores had already been inspected; current withheld-species assessment is exploratory, not a pristine prospective holdout. No Tp/Ms measurements enter the optimizer.',
    'objective':'shared parameters, robust soft_l1 log C,N,density + 0.3 asinh DIN; all future held-out observations excluded',
    'phase_assumption':'first sample at lights-on; phase sensitivity not fitted',
    'prior_parameters':asdict(PRIOR),'fitted_parameters':asdict(p),'bounds':dict(zip(PARAMS,zip(*BOUNDS))),
    'solutions':[z[1] for z in solutions],'jacobian_singular_values':sv.tolist(),
    'jacobian_condition':float(sv[0]/sv[-1]),'metrics':metrics,'predictions':predictions,'runtime_seconds':time.perf_counter()-t}
 # JSON has no nonstandard NaN tokens.
 def clean(x):
  if isinstance(x,float) and not np.isfinite(x):return None
  if isinstance(x,list):return [clean(z) for z in x]
  if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
  return x
 (outdir/'validation.json').write_text(json.dumps(clean(report),indent=2,allow_nan=False))
 (ROOT/'research/calibrated_physiology.json').write_text(json.dumps(asdict(p),indent=2))
if __name__=='__main__':main()
