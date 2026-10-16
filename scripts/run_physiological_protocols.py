"""Small carbon/energy challenge protocols with explicit pulse/dilution ledgers."""
from pathlib import Path
from dataclasses import replace
import json
import numpy as np
from splas2 import *
from splas2.model import INPUT_N,UPTAKE_N,GPP,RESP_P
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/protocols'

def main():
 OUT.mkdir(parents=True,exist_ok=True);p=Physiology(**json.load(open(ROOT/'research/calibrated_physiology.json')));records=[]
 for name,segments in [('dark_recovery',[(5,85,0),(5,0,0),(10,85,0)]),
                       ('N_pulse_recovery',[(8,85,0),(12,85,2)])]:
  cfg=Config(n_phyto=1,n_zoo=0,n_detritus=1,phyto_min=5,phyto_max=5,depth=1,surface_PAR=85)
  model=Model(cfg,p);y=model.initial(DIN_value=.1,phyto_N=.1,zoo_N=0)
  initial_budget=model.budgets(y,True);tt=[0.];yy=[y.copy()];t=0.;steps=0
  for duration,light,pulse in segments:
   m=Model(replace(cfg,surface_PAR=light),p)
   if pulse:
    y[0]+=pulse;y[m.ntr+INPUT_N]+=pulse;tt.append(t);yy.append(y.copy())
   r=simulate(m,duration,y0=y,max_step=.02,save_interval=.1,method='explicit');steps+=r.steps
   tt.extend((t+r.time[1:]).tolist());yy.extend(list(r.state[1:]));t+=duration;y=r.state[-1].copy()
  state=np.array(yy);times=np.array(tt);bud=np.array([model.budgets(z,True) for z in state])-initial_budget
  bc=np.array([model.biomass(z)[0][0,0] for z in state]);bn=np.array([model.biomass(z)[1][0,0] for z in state])
  np.savez_compressed(OUT/name,time=times,state=state,phyto_C=bc,phyto_N=bn,CN=bc/bn)
  records.append({'id':name,'segments_days_PAR_Npulse':segments,'steps':steps,'max_abs_mass_error_C':float(abs(bud[:,0]).max()),'max_abs_mass_error_N':float(abs(bud[:,1]).max()),'minimum_state':float(state[:,:model.ntr].min()),'final_CN':float(bc[-1]/bn[-1]),'interpretation':'Mechanistic energy/storage test, not an empirical trajectory.'})
 cfg=Config(n_phyto=1,n_zoo=0,n_detritus=1,phyto_min=5,phyto_max=5,depth=1,surface_PAR=85,dilution=.1,inflow_DIN=5)
 m=Model(cfg,p);r=simulate(m,100,y0=m.initial(phyto_N=.1,zoo_N=0),max_step=.02,save_interval=.5);r.save(OUT/'chemostat.npz')
 records.append({'id':'chemostat','inflow_DIN_uM':5,'dilution_d':.1,**r.diagnostics()})
 (OUT/'protocols.json').write_text(json.dumps(records,indent=2,allow_nan=False));print(json.dumps(records,indent=2))
if __name__=='__main__':main()
