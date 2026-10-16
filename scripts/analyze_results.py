"""Audit run identities, aggregate exact metrics, and make shareable scientific figures."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import json,hashlib,argparse,csv,platform,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from splas2 import Config,Physiology,Model
from splas2.metrics import normalized_spectrum
ROOT=Path(__file__).resolve().parents[1];RDIR=ROOT/'results/thesis_replay';FIG=ROOT/'results/figures'
COLORS=['#166b8c','#d16f30','#6b5294','#3a916f','#b54c5a','#94733b','#4c749f','#8b729e']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,
   'axes.prop_cycle':plt.cycler(color=COLORS),'savefig.dpi':180,'figure.constrained_layout.use':True})

def cn_decomposition(a,b):
 ca=a['spectrum_C'];na=a['spectrum_N'];cb=b['spectrum_C'];nb=b['spectrum_N']
 ra=np.divide(ca,na,out=np.zeros_like(ca),where=na>0);rb=np.divide(cb,nb,out=np.zeros_like(cb),where=nb>0)
 ra=np.where(na==0,rb,ra);rb=np.where(nb==0,ra,rb)
 wa=na/na.sum();wb=nb/nb.sum()
 intra=float(np.sum((wa+wb)/2*(rb-ra)));composition=float(np.sum((ra+rb)/2*(wb-wa)))
 delta=float(cb.sum()/nb.sum()-ca.sum()/na.sum())
 return {'total_CN_change':delta,'within_class_term':intra,'composition_term':composition,'identity_error':delta-intra-composition,
    'note':'Symmetric N-weighted decomposition of ratios of time-mean stocks. Absent-class quota is copied from the other case; not a unique mechanistic attribution.'}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--allow-partial',action='store_true');args=ap.parse_args()
 FIG.mkdir(parents=True,exist_ok=True);registry=json.load(open(ROOT/'research/experiment_registry.json'));source=(ROOT/'splas2/model.py').read_text()
 for extra in ['supplement_registry.json','carbon_boundary_registry.json','recommended_registry.json','recommended_refinement_registry.json','export_intervention_registry.json']:
  supplement=ROOT/'research'/extra
  if supplement.exists():registry+=json.load(open(supplement))
  elif not args.allow_partial:raise RuntimeError('Unfinished study: missing registry '+extra)
 if not args.allow_partial and len(registry)!=117:raise RuntimeError('Unfinished study: expected117 registered cases')
 cases={};missing=[]
 for c in registry:
  p=RDIR/(c['id']+'.json')
  if not p.exists():missing.append(c['id']);continue
  try:r=json.load(open(p))
  except json.JSONDecodeError:missing.append(c['id']);continue
  sig=hashlib.sha256((json.dumps(c,sort_keys=True)+source).encode()).hexdigest()
  if r.get('run_signature')!=sig:missing.append(c['id']);continue
  cases[c['id']]=r
 if missing and not args.allow_partial:raise RuntimeError('Missing current runs: '+','.join(missing))
 rows=[{'id':k,'family':v['family'],'days':v['days'],**v['summary']} for k,v in cases.items()]
 if rows:
  keys=list(rows[0]);f=(ROOT/'research/experiment_summary.csv').open('w',newline='');w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows);f.close()
 out={'current_cases':len(cases),'requested_cases':len(registry),'missing_cases':missing,
  'decadal_case_count':sum(c['days']>=3650 for c in cases.values()),'annual_continuation_count':sum(c['days']<3650 for c in cases.values()),
  'model_sha256':hashlib.sha256(source.encode()).hexdigest(),'mass_C_max_relative':max([r['summary']['max_relative_mass_error_C'] for r in cases.values()],default=0),
  'mass_N_max_relative':max([r['summary']['max_relative_mass_error_N'] for r in cases.values()],default=0),
  'minimum_state':min([r['summary']['minimum_live_state'] for r in cases.values()],default=0),
  'negative_state_cases':[k for k,v in cases.items() if v['summary']['minimum_live_state']<0],
  'drifting_cases':[k for k,v in cases.items() if abs(v['summary']['last_year_N_drift_fraction'] or 0)>.1],
  'elapsed_compute_seconds':sum(v['summary']['runtime_seconds'] for v in cases.values()),'numerical_comparisons':[], 'stoichiometric_decompositions':[]}
 def arr(key):return np.load(RDIR/(key+'.npz'))
 def figsave(fig,name):fig.savefig(FIG/(name+'.png'),facecolor='white');fig.savefig(FIG/(name+'.pdf'));plt.close(fig)
 # Culture trajectories: measurements, shared fitted dynamics, Redfield null.
 v=json.load(open(ROOT/'results/cultures/validation.json'));D=json.load(open(ROOT/'data/processed/cultures.json'))['species']
 fig,ax=plt.subplots(2,4,figsize=(11.8,6.1))
 for j,key in enumerate(['Ot','Tw','Tp','Ms']):
  pr=v['predictions'][key];dd=D[key];t=np.array(pr['times']);obs=np.array([[np.nan if x is None else x for x in z] for z in pr['observed']]);pred=np.array(pr['predicted'])
  sd=np.array([o['C_pg_sd'] for o in dd['observations']]);ax[0,j].errorbar(t,obs[:,0],sd,fmt='o',color='#273b46',label='Dati ± SD')
  ax[0,j].plot(t,pred[:,0],label='Modello C:N');ax[0,j].set_title(key+' | '+('stima' if j<2 else 'valutazione'))
  ax[0,j].set_ylabel('C (pg/cellula)' if j==0 else '')
  ax[1,j].errorbar(t,pr['CN_observed'],[o['CN_molar_sd'] for o in dd['observations']],fmt='o',color='#273b46')
  ax[1,j].plot(t,pr['CN_predicted']);ax[1,j].axhline(6.625,color='#a9a9a9',ls='--',label='Redfield')
  ax[1,j].set_ylabel('C:N molare' if j==0 else '');ax[1,j].set_xlabel('Giorni dalla diluizione')
 ax[0,0].legend(frameon=False,fontsize=8);ax[1,0].legend(frameon=False,fontsize=8);figsave(fig,'01_cultures')
 fig,ax=plt.subplots(figsize=(8,3.5));x=np.arange(4);vals=v['metrics']
 # Keep original data order Ot,Tw,Tp,Ms, same as extraction dictionary.
 for j,(key,label) in enumerate([('CN_RMSLE','Dinamico'),('CN_Redfield_RMSLE','Redfield'),('CN_initial_ratio_RMSLE','C:N iniziale costante')]):
  ax.bar(x+(j-1)*.24,[r[key] for r in vals],.23,label=label)
 ax.set_xticks(x,[r['species'] for r in vals]);ax.set_ylabel('RMSE log(C:N predetto/osservato)');ax.legend(frameon=False,fontsize=8);figsave(fig,'02_culture_scores')
 supply=['D01_supply_0.001','D01_supply_0.01','D01_supply_0.1','D01_supply_0.5','D00_reference','D01_supply_2','D01_supply_5','D01_supply_7','D01_supply_10']
 have=[k for k in supply if k in cases]
 if have:
  fig,ax=plt.subplots(2,2,figsize=(10,6.6));S=np.array([cases[k]['config']['nutrient_supply'] for k in have]);ss=[cases[k]['summary'] for k in have]
  for key in have:
   a=arr(key);ax[0,0].plot(a['sizes'],normalized_spectrum(a['sizes'],a['spectrum_N']),label=str(cases[key]['config']['nutrient_supply']))
  ax[0,0].set(xscale='log',yscale='log',xlabel='ESD (µm)',ylabel='N / Δlog10 ESD (mmol/m²)');ax[0,0].legend(title='S (mmolN/m³/d)',fontsize=7,ncol=3,frameon=False)
  ax[0,1].plot(S,[s['N_weighted_geomean_ESD_um'] for s in ss],'o-');ax[0,1].set(xscale='log',xlabel='S (mmolN/m³/d)',ylabel='ESD geometrico pesato N (µm)')
  ax[1,0].plot([s['phyto_N_mmol_m2'] for s in ss],[s['phyto_evenness_N'] for s in ss],'o-');ax[1,0].set(xscale='log',xlabel='Biomassa P (mmolN/m²)',ylabel='Evenness per taglia, peso N')
  ax[1,1].plot(S,[s['TTE_flux_C'] for s in ss],'o-',label='TTE di flusso');ax[1,1].plot(S,[s['Z_to_P_stock_C'] for s in ss],'s--',label='Stock Z:P')
  ax[1,1].set(xscale='log',xlabel='S (mmolN/m³/d)',ylabel='Rapporto C');ax[1,1].legend(frameon=False,fontsize=8);figsave(fig,'03_bottom_up')
 carbon=[k for k in cases if k.startswith('D14_Cbath_supply_')]
 if carbon:
  carbon.sort(key=lambda k:cases[k]['config']['nutrient_supply'])
  fig,ax=plt.subplots(1,3,figsize=(11,3.8))
  for keys,label in [(have,'Rilassamento DIC0.01/d'),(carbon,'Reservoir DIC1/d')]:
   ss=[cases[k]['summary'] for k in keys];x=[cases[k]['config']['nutrient_supply'] for k in keys]
   for a,field,ylabel in zip(ax,['N_weighted_geomean_ESD_um','phyto_CN_inventory_ratio','last_year_min_DIC'],['ESD geometrico peso N (µm)','C:N molare','Minimo DIC (mmolC/m³)']):
    a.plot(x,[s[field] for s in ss],'o-',label=label);a.set(xscale='log',xlabel='S (mmolN/m³/d)',ylabel=ylabel)
  ax[0].legend(frameon=False,fontsize=7);figsave(fig,'10_carbon_boundary')
 if 'D00_reference' in cases:
  reference=cases['D00_reference']['summary']
  for k in ['D07_size_bins_35','D07_size_bins_140','D12_dt_0.025','D12_dt_0.0125']:
   if k in cases:
    alt=cases[k]['summary'];out['numerical_comparisons'].append({'reference':'D00_reference','alternative':k,**{x+'_relative_change':(alt[x]/reference[x]-1) if reference[x] else None for x in ['phyto_N_mmol_m2','phyto_CN_inventory_ratio','N_weighted_geomean_ESD_um','NPP_C_mmol_m2_d']}})
 if 'F00_open_reference' in cases:
  ref=cases['F00_open_reference']['summary']
  for key in ['F07_layers_20','F07_layers_80','F08_detritus_bins_4','F08_detritus_bins_8','F08_detritus_bins_16','F09_dt_0.025','F09_dt_0.0125']:
   if key in cases:
    alt=cases[key]['summary'];out['numerical_comparisons'].append({'reference':'F00_open_reference','alternative':key,**{x+'_relative_change':(alt[x]/ref[x]-1) if ref[x] else None for x in ['phyto_N_mmol_m2','phyto_CN_inventory_ratio','N_weighted_geomean_ESD_um','NPP_C_mmol_m2_d','bottom_export_C_mmol_m2_d']}})
  fig,ax=plt.subplots(1,4,figsize=(11,4.8))
  for key,label in [('F00_open_reference','Aperto'),('F00_closed_reference','Chiuso'),('F08_detritus_bins_8','8 classi detrito')]:
   if key not in cases:continue
   a=arr(key);z=a['depth'];ax[0].plot(a['profile_N'],z,label=label);ax[1].plot(a['profile_DIN'],z);ax[2].plot(a['profile_CN'],z)
   data=a['profile_P_size_N'];total=data.sum(axis=0);mean=np.exp((np.log(a['sizes'])[:,None]*data).sum(axis=0)/np.maximum(total,1e-300));ax[3].plot(mean,z)
  for x in ax:x.invert_yaxis();x.set_ylim(200,0)
  ax[0].set(xlabel='P (mmolN/m³)',ylabel='Profondità (m)');ax[1].set_xlabel('DIN (mmolN/m³)');ax[2].set_xlabel('C:N molare');ax[3].set_xlabel('ESD pesato N (µm)');ax[0].legend(frameon=False,fontsize=8);figsave(fig,'04_vertical_profiles')
  # Independent within-class vs composition decomposition.
  for a,b in [('F01_open_light_0.7','F01_open_light_1.3'),('F01_closed_light_0.7','F01_closed_light_1.3')]:
   if a in cases and b in cases:out['stoichiometric_decompositions'].append({'reference':a,'alternative':b,**cn_decomposition(arr(a),arr(b))})
 season=[k for k in ['F02_semiannual','F02_annual','F02_constant_matched'] if k in cases]
 if season:
  fig,ax=plt.subplots(2,2,figsize=(10,6))
  for k in season:
   a=arr(k);t=a['time'];ix=t>=t[-1]-730;xx=t[ix]-(t[-1]-730);label={'F02_semiannual':'Semestrale storica','F02_annual':'Annuale','F02_constant_matched':'Costante, stessa dose'}[k]
   for sub,key,unit in [(ax[0,0],'phyto_C','P (mmolC/m²)'),(ax[0,1],'phyto_N','P (mmolN/m²)'),(ax[1,0],'CN','C:N molare'),(ax[1,1],'zoo_C','Z (mmolC/m²)')]:
    sub.plot(xx,a[key][ix],label=label);sub.set_ylabel(unit);sub.set_xlabel('Giorni negli ultimi 2 anni')
  ax[0,0].legend(frameon=False,fontsize=8);figsave(fig,'05_equal_light_dose')
 # Deterministic sensitivity plot - variability, not confidence intervals.
 sensitivity=[k for k in ['F00_open_reference','F04_mixing_0.1','F04_mixing_10','F04_variable_mixing','F05_no_grazers','F05_no_shading','F06_PAR_fraction_043','F06_DIC_exchange','F08_detritus_bins_4','F08_detritus_bins_8','F08_detritus_bins_16','F10_product_colimitation','F10_prior_physiology','F10_selective_remin'] if k in cases]
 if sensitivity:
  fig,ax=plt.subplots(1,2,figsize=(11,5.3));yy=np.arange(len(sensitivity))
  ax[0].barh(yy,[cases[k]['summary']['NPP_C_mmol_m2_d'] for k in sensitivity]);ax[1].barh(yy,[cases[k]['summary']['bottom_export_C_mg_m2_d'] for k in sensitivity])
  ax[0].set_yticks(yy,[k.replace('F00_','').replace('F04_','').replace('F05_','').replace('F06_','').replace('F08_','').replace('F10_','') for k in sensitivity],fontsize=8);ax[1].set_yticks(yy,[])
  ax[0].set_xlabel('NPP (mmolC/m²/d)');ax[1].set_xlabel('Export al fondo (mgC/m²/d)');figsave(fig,'06_structural_sensitivity')
 if out['numerical_comparisons']:
  comparisons=out['numerical_comparisons'];fig,ax=plt.subplots(figsize=(10,4.5));xx=np.arange(len(comparisons))
  for i,field in enumerate(['phyto_N_mmol_m2','phyto_CN_inventory_ratio','NPP_C_mmol_m2_d','bottom_export_C_mmol_m2_d']):
   ax.plot(xx,[100*c.get(field+'_relative_change',np.nan) if c.get(field+'_relative_change') is not None else np.nan for c in comparisons],'o-',label=field.split('_mmol')[0])
  ax.axhline(0,color='gray',lw=.8);ax.set_xticks(xx,[c['alternative'] for c in comparisons],rotation=45,ha='right',fontsize=8);ax.set_ylabel('Variazione rispetto al riferimento (%)');ax.legend(frameon=False,fontsize=8);figsave(fig,'07_refinement')
 long_stats=[]
 for key in ['D00_reference','D12_dt_0.025','D12_dt_0.0125','F00_open_reference','F09_dt_0.025','F09_dt_0.0125']:
  if key not in cases:continue
  a=arr(key);c=cases[key];ntr=4+3*c['config']['n_phyto']+c['config']['n_zoo']+2*c['config']['n_detritus'];dz=c['config']['depth']/c['config']['n_layers']
  ledger=(a['final_state'][ntr:]-a['initial_state'][ntr:]).sum(axis=1)*dz/c['days'];t=a['time'];ix=t>=t[-1]-730;tt=t[ix]
  meanC=np.trapezoid(a['phyto_C'][ix],tt)/(tt[-1]-tt[0]);meanN=np.trapezoid(a['phyto_N'][ix],tt)/(tt[-1]-tt[0])
  long_stats.append({'id':key,'two_year_CN_ratio':float(meanC/meanN),'two_year_phyto_N':float(meanN),'whole_run_mean_NPP_C':float(ledger[0]-ledger[1]),'whole_run_mean_export_C':float(ledger[6])})
 out['time_window_diagnostic']=long_stats
 if all(k in cases for k in ['R01_same_state_dt_0.0125','R02_grid80','R02_grid160']):
  keys=['R01_same_state_dt_0.0125','R02_grid80','R02_grid160'];cells=[40,80,160]
  fig,axes=plt.subplots(2,2,figsize=(9.5,6.1))
  for ax,field,label in zip(axes.ravel(),['phyto_N_mmol_m2','N_weighted_geomean_ESD_um','NPP_C_mmol_m2_d','bottom_export_C_mg_m2_d'],['P_N (mmolN/m²)','ESD geometrico pesato N (µm)','NPP C (mmolC/m²/d)','Export C a200m (mgC/m²/d)']):
   ax.plot(cells,[cases[k]['summary'][field] for k in keys],'o-');ax.set_xscale('log',base=2);ax.set_xticks(cells,[str(n) for n in cells]);ax.set(xlabel='Numero di celle, stesso stato iniziale',ylabel=label)
  figsave(fig,'13_nested_vertical_refinement')
 # A stable phytoplankton stock does not imply a stationary inorganic inventory.
 out['last_window_element_balance']=[]
 for key,c in cases.items():
  a=arr(key);d=a['diag_final_year'];dz=c['config']['depth']/c['config']['n_layers']
  out['last_window_element_balance'].append({'id':key,'window_days':c['summary']['window_end_day']-c['summary']['window_start_day'],
   'C_external_net_mmol_m2_d':float(d[10]-d[12]),'N_external_input_mmol_m2_d':float(d[11]),
   'N_external_output_mmol_m2_d':float(d[7]+d[9]+d[13]),
   'N_inventory_accumulation_mmol_m2_d':float(d[11]-d[7]-d[9]-d[13]),
   'DIN_whole_run_change_mean_mmol_m2_d':float((a['final_state'][0]-a['initial_state'][0]).sum()*dz/c['days']),
   'interpretation':'Net element retention from integrated external ledgers; consistent with conservation, distinct from phytoplankton-stock stationarity.'})
 mechanisms=ROOT/'results/mechanisms/mechanisms.json'
 if mechanisms.exists():
  mech=json.load(open(mechanisms));out['Rstar_oracle']=mech['fixed_CN_assembly'][0]
  out['mechanism_case_counts']={k:len(mech[k]) for k in ['fixed_CN_assembly','culture_phase','physiology_grid','mixing_experiments']}
  mix=mech['mixing_experiments'];fig,ax=plt.subplots(1,2,figsize=(9.8,4.3))
  for grazing,marker in [(False,'o'),(True,'s')]:
   for light in [40,120,360]:
    ss=[m for m in mix if m['grazers']==grazing and m['surface_PAR']==light]
    if ss:
     ax[0].plot([m['mixing_depth_m_prescribed'] for m in ss],[m['phyto_N_mmol_m2'] for m in ss],marker+'-',label=f'PAR {light}, Z{int(grazing)}')
     ax[1].plot([m['mixing_depth_m_prescribed'] for m in ss],[m['fraction_defined_critical_depth'] for m in ss],marker+'-')
  ax[0].set(xlabel='Profondità mixing prescritta (m)',ylabel='P medio ultimi180d (mmolN/m²)');ax[1].set(xlabel='Profondità mixing prescritta (m)',ylabel='Frazione campioni con Dcr definita');ax[0].legend(frameon=False,fontsize=8,ncol=2);figsave(fig,'08_mixing_mechanism')
  out['culture_phase_uncertainty']={k:{'min_final_CN':min(s['final_CN'] for s in mech['culture_phase'] if s['species']==k),'max_final_CN':max(s['final_CN'] for s in mech['culture_phase'] if s['species']==k)} for k in D}
  fig,ax=plt.subplots(1,2,figsize=(10,4.2))
  for nz in [0,1,2,8,32,70]:
   a=np.load(ROOT/f'results/mechanisms/fixed_assembly_{nz}.npz');ax[0].plot(a['sizes'],normalized_spectrum(a['sizes'],a['final'][2:72]),label=str(nz))
   key=f'D11_quota_assembly_{nz}'
   if key in cases:
    b=arr(key);ax[1].plot(b['sizes'],normalized_spectrum(b['sizes'],b['spectrum_N'])/100,label=str(nz))
  for a in ax:a.set(xscale='log',yscale='log',xlabel='ESD (µm)',ylabel='P (mmolN/m³ / Δlog10 ESD)')
  ax[0].set_title('Controllo Monod C:N fisso');ax[1].set_title('Riserve C:N dinamiche');ax[0].legend(title='Numero predatori',frameon=False,fontsize=8,ncol=2);figsave(fig,'09_predator_assembly')
 # Analytical consequence of the declared consumer stoichiometry.
 cn=np.linspace(4,80,300);yieldC=np.minimum(.396/1.2,5.5*.67/cn)
 fig,ax=plt.subplots(figsize=(7.8,3.6));ax.plot(cn,yieldC);ax.axvline(5.5*.67/(.396/1.2),color='gray',ls='--')
 ax.set(xlabel='C:N molare della preda',ylabel='Produzione Z_C / C ingerito',title='Vincolo stechiometrico del modello, non dato empirico')
 figsave(fig,'11_food_quality_hypothesis')
 # Revised reference and controlled annual interventions share an initial state.
 out['annual_intervention_comparisons']=[]
 if 'R01_same_state_dt_0.0125' in cases:
  ref=cases['R01_same_state_dt_0.0125']['summary']
  keys=[k for k in cases if k.startswith(('R04_','R05_'))]
  fields=['phyto_N_mmol_m2','phyto_CN_inventory_ratio','NPP_C_mmol_m2_d','bottom_export_C_mmol_m2_d','N_weighted_geomean_ESD_um']
  for key in keys:
   alt=cases[key]['summary']
   out['annual_intervention_comparisons'].append({'reference':'R01_same_state_dt_0.0125','alternative':key,
    'relative_changes':{field:alt[field]/ref[field]-1 if ref[field] and alt[field] is not None else None for field in fields}})
  if keys:
   fig,ax=plt.subplots(1,2,figsize=(10.6,4.8));yy=np.arange(len(keys))
   for a,field,title in zip(ax,['NPP_C_mmol_m2_d','bottom_export_C_mmol_m2_d'],['NPP C','Export C al fondo200m']):
    vals=[100*r['relative_changes'][field] for r in out['annual_intervention_comparisons']]
    a.barh(yy,vals);a.axvline(0,color='gray',lw=.8);a.set_xlabel('Variazione annuale rispetto al controllo (%)');a.set_title(title)
   ax[0].set_yticks(yy,[k.replace('R04_','').replace('R05_','') for k in keys],fontsize=8);ax[1].set_yticks(yy,[])
   figsave(fig,'12_annual_interventions')
 # Deterministic pattern checks, not significance tests or observations.
 out['pattern_assessment']=[]
 groups={
  'Fianchini_constant_open':['F01_open_light_0.7','F01_open_light_0.8','F01_open_light_0.9','F00_open_reference','F01_open_light_1.1','F01_open_light_1.2','F01_open_light_1.3'],
  'Fianchini_constant_closed':['F01_closed_light_0.7','F01_closed_light_0.8','F01_closed_light_0.9','F00_closed_reference','F01_closed_light_1.1','F01_closed_light_1.2','F01_closed_light_1.3'],
  'Fianchini_semiannual_200m':['F02_seasonal_light_0.7','F02_seasonal_light_0.8','F02_seasonal_light_0.9','F02_semiannual','F02_seasonal_light_1.1','F02_seasonal_light_1.2','F02_seasonal_light_1.3'],
  'Fianchini_semiannual_100m':['F03_100m_light_'+str(q) for q in [.7,.8,.9,1.,1.1,1.2,1.3]],
  'DePasquale_N_supply':supply,
  'N_supply_carbon_bath':['D14_Cbath_supply_'+str(q) for q in [.001,.01,.1,.5,1,2,5,7,10]]}
 for name,keys in groups.items():
  if any(k not in cases for k in keys):continue
  x=np.array([cases[k]['config']['nutrient_supply' if 'supply' in name else 'surface_PAR'] for k in keys]);rec={'family':name,'cases':keys,'metrics':{},'interpretation':'Deterministic rank/monotonicity and descriptive linearity only; not a statistical validation against data.'}
  for field in ['phyto_N_mmol_m2','phyto_CN_inventory_ratio','phyto_evenness_N','N_weighted_geomean_ESD_um','N_weighted_ESD95_um','bottom_export_C_mmol_m2_d','e_ratio_at_bottom']:
   values=[cases[k]['summary'][field] for k in keys]
   if any(v is None for v in values):continue
   y=np.array(values);constant=np.ptp(y)<=1e-12
   coef=np.polyfit(x,y,1);fit=np.polyval(coef,x);ss=float(np.sum((y-y.mean())**2))
   rec['metrics'][field]={'spearman_rho':None if constant else float(spearmanr(x,y).statistic),
    'nondecreasing':bool(np.all(np.diff(y)>=-1e-12)),'nonincreasing':bool(np.all(np.diff(y)<=1e-12)),
    'linear_slope':float(coef[0]),'linear_R2':None if ss<=1e-24 else float(1-np.sum((y-fit)**2)/ss),'values':values}
  out['pattern_assessment'].append(rec)
 # Store allometry checks and environment in human-readable, versionable files.
 out['culture_metrics']=v['metrics'];out['culture_fit_bounds']=v['bounds'];out['culture_jacobian_condition']=v['jacobian_condition']
 coverage=[]
 for key,d in D.items():
  pr=np.array(v['predictions'][key]['predicted']);obs=d['observations'];rec={'species':key}
  for j,field in enumerate(['C_pg','N_pg','cells_ml']):
   residual=np.array([(pr[i,j]-o[field])/o[field+'_sd'] for i,o in enumerate(obs) if i>0 and o[field+'_sd'] and o[field+'_sd']>0])
   rec[field+'_fraction_within_2replicateSD']=float(np.mean(abs(residual)<=2));rec[field+'_max_abs_SD_residual']=float(max(abs(residual)))
  coverage.append(rec)
 out['culture_descriptive_replicate_SD_check']=coverage
 out['assessed_domain']='Atom-conservative C:N model; N-starvation cultures and idealized thesis counterpart experiments. No site-level field/export validation.'
 (ROOT/'research/validation_summary.json').write_text(json.dumps(out,indent=2,allow_nan=False))
 env={'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'scipy':__import__('scipy').__version__,'matplotlib':matplotlib.__version__,'numba':__import__('numba').__version__,'model_sha256':out['model_sha256']}
 (ROOT/'research/reference_environment.json').write_text(json.dumps(env,indent=2))
 print(json.dumps({k:out[k] for k in ['current_cases','requested_cases','mass_C_max_relative','mass_N_max_relative','negative_state_cases','missing_cases']},indent=2))
if __name__=='__main__':main()
