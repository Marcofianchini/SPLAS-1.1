"""Read the CC-BY data without changing the Excel workbook; deterministic extraction."""
from pathlib import Path
import json, hashlib, warnings
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
RAW=ROOT/'data/raw/Liefer2019/Final_MASS_Database.xlsx'
SPECIES={'Ot':('Ostreococcus tauri',1.8,18.,2.,'calibration'),
         'Tw':('Thalassiosira weissflogii',1630.,18.,3.,'calibration'),
         'Tp':('Thalassiosira pseudonana',158.,18.,1.,'held_out'),
         'Ms':('Micromonas sp.',1.8,6.,6.,'held_out')}
# Published Table 1, rounded numerical facts, DOI 10.3389/fmicb.2019.00763.
TABLE={'Ot':[.277,.0536,.362,.0269], 'Tw':[139.,26.6,368.,8.30],
       'Tp':[15.1,2.70,26.7,1.29], 'Ms':[.582,.105,.634,.0604]}

def main():
 with warnings.catch_warnings():
  warnings.simplefilter('ignore',UserWarning)
  df=pd.read_excel(RAW,sheet_name='Full Data Summary')
 df=df[df.Species.notna()].copy()
 cols={'Days in N-free Media':'day','Growth Stage':'stage','Sample ID':'sample_id',
       'Cell Density':'cells_ml','Growth Rate':'growth_d','C':'C_pg','N':'N_pg','DIN':'DIN_uM','C:N':'CN_reported'}
 out=df[['Species',*cols]].rename(columns=cols)
 out['CN_molar']=out.C_pg/out.N_pg*14/12
 out['excluded_contamination']=((out.Species.isin(['Ostreococcus tauri','Thalassiosira pseudonana']))&(out.stage=='Late Stationary'))
 out.to_csv(ROOT/'data/processed/culture_replicates.csv',index=False)
 data={};checks=[]
 for key,(name,vol,temp,start,split) in SPECIES.items():
  g=out[out.Species==name];valid=g[~g.excluded_contamination]
  ref=g[g.day==0]
  means=valid.groupby(['day','stage'],sort=True)[['C_pg','N_pg','cells_ml','DIN_uM','CN_molar']].mean().reset_index()
  sd=valid.groupby(['day','stage'],sort=True)[['C_pg','N_pg','cells_ml','DIN_uM','CN_molar']].std().reset_index()
  obs=[]
  for i,r in means.iterrows():
   if r.day<start:continue
   d={'day':float(r.day),'t':float(r.day-start),'stage':r.stage}
   for k in ['C_pg','N_pg','cells_ml','DIN_uM','CN_molar']:
    d[k]=float(r[k]) if pd.notna(r[k]) else None
    d[k+'_sd']=float(sd.iloc[i][k]) if pd.notna(sd.iloc[i][k]) else None
   obs.append(d)
  data[key]={'name':name,'volume_um3':vol,'temperature_C':temp,'split':split,'start_day':start,
             'reference_C_pg':float(ref.C_pg.mean()),'reference_N_pg':float(ref.N_pg.mean()),'observations':obs}
  mid=valid[valid.stage=='Mid Stationary'];midday=mid.day.max();mid=mid[mid.day==midday]
  actual=[ref.C_pg.mean(),ref.N_pg.mean(),mid.C_pg.mean(),mid.N_pg.mean()]
  discrepancy=max(abs(np.array(actual)/TABLE[key]-1))
  checks.append({'species':key,'published':[float(x) for x in TABLE[key]],'raw_means':[float(x) for x in actual],
                 'max_relative_difference':float(discrepancy),'agrees_to_rounding_1pct':bool(discrepancy<.01)})
 payload={'dataset_doi':'10.1594/PANGAEA.900248','paper_doi':'10.3389/fmicb.2019.00763',
   'license':'CC-BY-4.0; Liefer et al. (2026), original experiments 2016, article 2019',
   'curation_caveat':'No final formal author approval per PANGAEA. Raw units label micromolar (mM) is contradictory; article establishes micromolar (= mmol m^-3). O. tauri day2 C,N exactly repeat day0. Dates not used to infer time.',
   'protocol':'Conditional predictions after unquantified dilution; initial cell C,N,density,DIN observed. Day0 replete C,N and volumes are reference traits, including held-out species. No subsequent held-out data used for estimation. 85 umol photons m^-2 s^-1, 12:12 light:dark. Initial phase assumed lights-on; evaluate phase sensitivity.',
   'replicate_rows':len(out),'excluded_contamination_rows':int(out.excluded_contamination.sum()),
   'raw_sha256':hashlib.sha256(RAW.read_bytes()).hexdigest(),'table1_crosscheck':checks,'species':data}
 (ROOT/'data/processed/cultures.json').write_text(json.dumps(payload,indent=2,allow_nan=False))
 print(json.dumps({'rows':len(out),'excluded':payload['excluded_contamination_rows'],'table1':checks},indent=2))
 if not all(c['agrees_to_rounding_1pct'] for c in checks):raise RuntimeError('Published table mismatch')
if __name__=='__main__':main()
