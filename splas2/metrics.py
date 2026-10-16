"""Observable definitions. Stocks are never silently interpreted as efficiencies."""
import numpy as np
from .model import *

def normalized_spectrum(sizes,biomass):
    """Inventory per unit log10 ESD, with endpoints treated as half-bins."""
    x=np.log10(sizes)
    if len(x)==1:return np.asarray(biomass)
    edges=np.r_[x[0],.5*(x[1:]+x[:-1]),x[-1]]
    return np.asarray(biomass)/np.diff(edges)

def evenness(b):
    b=np.asarray(b);total=b.sum()
    if total<=1e-12 or len(b)<2:return None
    p=b[b>0]/total
    p=p[p>0] # subnormal biomass can underflow during normalization
    return float(-(p*np.log(p)).sum()/np.log(len(b)))

def summarize(result,window=365.):
    r=result;m=r.model;t=r.time;y=r.state
    k=max(0,int(np.searchsorted(t,max(0,t[-1]-window),side='left')))
    if k==len(t)-1:k=max(0,k-1)
    tt=t[k:];dt=tt[-1]-tt[0]
    bc=m.traits[:,3][None,:,None]*y[:,m.a]+y[:,m.cr]+m.physiology.reserve_N_CN*y[:,m.nr]
    bn=m.traits[:,4][None,:,None]*y[:,m.a]+y[:,m.nr]
    meanC=np.trapezoid(bc[k:],tt,axis=0)/dt
    meanN=np.trapezoid(bn[k:],tt,axis=0)/dt
    psC=meanC.sum(axis=1)*m.dz;psN=meanN.sum(axis=1)*m.dz
    zmean=np.trapezoid(y[k:,m.z],tt,axis=0)/dt
    zsC=zmean.sum(axis=1)*m.dz
    d=((y[-1,m.ntr:]-y[k,m.ntr:]).sum(axis=1))*m.dz/dt
    C=psC.sum();N=psN.sum();Z=zsC.sum()
    npp=float(d[GPP]-d[RESP_P]);export=float(d[EXPORT_C])
    weighted=float(np.exp((np.log(m.sizes)*psN).sum()/N)) if N>1e-12 else None
    cc=psN.cumsum()/N if N>1e-12 else np.zeros(m.nP)
    d95=float(m.sizes[min(m.nP-1,np.searchsorted(cc,.95))]) if N>1e-12 else None
    profile=meanN.sum(axis=0)
    dbm=float(m.depth[np.argmax(profile)]) if N>1e-12 else None
    cn_time=bc.sum(axis=(1,2))/np.maximum(bn.sum(axis=(1,2)),1e-300)
    biomassN=bn.sum(axis=(1,2))*m.dz
    # Numerical steady/periodic drift: last-year vs previous-year mean.
    prior=(t>=t[-1]-2*window)&(t<t[k]);drift=None
    if prior.sum()>1:
      v=np.trapezoid(biomassN[prior],t[prior])/(t[prior][-1]-t[prior][0])
      drift=float((N-v)/max(1e-12,v))
    summary={**r.diagnostics(),'window_start_day':float(tt[0]),'window_end_day':float(tt[-1]),
      'phyto_C_mmol_m2':float(C),'phyto_N_mmol_m2':float(N),'zoo_C_mmol_m2':float(Z),
      'phyto_CN_inventory_ratio':float(C/N) if N>1e-12 else None,
      'N_weighted_geomean_ESD_um':weighted,'N_weighted_ESD95_um':d95,
      'phyto_evenness_N':evenness(psN),'phyto_evenness_C':evenness(psC),'zoo_evenness_C':evenness(zsC),
      'NPP_C_mmol_m2_d':npp,'GPP_C_mmol_m2_d':float(d[GPP]),'Z_production_C_mmol_m2_d':float(d[Z_PROD]),
      'TTE_flux_C':float(d[Z_PROD]/npp) if npp>1e-12 else None,
      'Z_to_P_stock_C':float(Z/C) if C>1e-12 else None,
      'bottom_export_C_mmol_m2_d':export,'bottom_export_C_mg_m2_d':12*export,
      'bottom_export_N_mmol_m2_d':float(d[EXPORT_N]),
      'export_CN':float(export/d[EXPORT_N]) if d[EXPORT_N]>1e-12 else None,
      'e_ratio_at_bottom':float(export/npp) if npp>1e-12 else None,
      'DBM_depth_m':dbm,'last_year_N_drift_fraction':drift,
      'last_year_min_DIC':float(y[k:,DIC].min()),'last_year_max_DIN':float(y[k:,DIN].max())}
    arrays={'time':t,'phyto_N':biomassN,'phyto_C':bc.sum(axis=(1,2))*m.dz,
      'zoo_C':y[:,m.z].sum(axis=(1,2))*m.dz,'CN':cn_time,
      'sizes':m.sizes,'zsizes':m.zsizes,'spectrum_N':psN,'spectrum_C':psC,'spectrum_Z_C':zsC,
      'depth':m.depth,'profile_N':profile,'profile_C':meanC.sum(axis=0),'profile_DIN':np.trapezoid(y[k:,DIN],tt,axis=0)/dt,
      'profile_CN':meanC.sum(axis=0)/np.maximum(profile,1e-300),'profile_P_size_N':meanN,'diag_final_year':d,
      'final_state':y[-1],'initial_state':y[0]}
    return summary,arrays
