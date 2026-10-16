"""Read-only process observations for conservative 0D stoichiometry experiments.

The numerical model is deliberately unchanged. These expressions reconstruct
instantaneous biological transfers and are checked against its independent RHS.
Rates use mmol m^-3 d^-1; sizes are ESD, not changing cell volume.
"""
from __future__ import annotations
import numpy as np
from .model import DIN, DIC


def observe(model, state, time=0.):
    """Resolve physiology and actual ingested C/N, for one homogeneous box."""
    if model.config.n_layers != 1:
        raise ValueError("This observation operator is defined for 0D boxes")
    m=model; p=m.physiology; c=m.config; tr=m.traits
    y=np.asarray(state); a=y[m.a,0]; cr=y[m.cr,0]; nr=y[m.nr,0]
    bc,bn=(z[:,0] for z in m.biomass(y)); temp=p.Q10**((p.temperature-18)/10)
    irradiance=float(m.light(time,y)[0][0]); tiny=1e-300
    nstatus=np.minimum(1.,bn/np.maximum(tr[:,2]*a,tiny))
    theta=p.theta_min+(p.theta_max-p.theta_min)*nstatus/(1+irradiance/p.acclimation_PAR)
    ccap=np.maximum(tr[:,8]*a-tr[:,3]*a-p.reserve_N_CN*nr,tiny)
    cf=1-np.minimum(1.,cr/ccap)**2
    nf=1-np.minimum(1.,bn/np.maximum(2*tr[:,2]*a,tiny))**2
    light=-np.expm1(-p.alpha_Chl*theta*irradiance/np.maximum(tr[:,6]/tr[:,1]*nstatus,tiny))
    fix=a*tr[:,6]*temp*light*y[DIC,0]/(y[DIC,0]+p.half_DIC)*nstatus*cf
    uptake=a*tr[:,7]*temp*y[DIN,0]/(p.half_N_scale*tr[:,0]+y[DIN,0])*nf
    uptake=np.minimum(uptake,p.mobilization*temp*cr/(p.cost_N+p.reserve_N_CN))
    maintenance=p.maintenance*temp*tr[:,3]*a
    paid=np.minimum(maintenance,p.mobilization*temp*cr)
    starve=np.maximum(0.,maintenance-paid)/np.maximum(bc,tiny)
    lc=cr/np.maximum(bc,tiny); ln=nr/np.maximum(bn,tiny)
    lim=np.minimum(lc,ln) if p.colimitation==0 else lc*ln
    division=a*tr[:,5]*temp*lim
    division=np.minimum(division,p.mobilization*temp*cr/((1+p.cost_growth)*tr[:,3]-p.reserve_N_CN*tr[:,4]))
    division=np.minimum(division,p.mobilization*temp*nr/tr[:,4])
    available=m.phi.T@bc
    factor=26*m.zsizes**(-.4)*temp*y[m.z,0]/(c.grazing_half_C+available)
    grazing=m.phi@factor
    ingestion_C=available*factor; ingestion_N=(m.phi.T@bn)*factor
    carbon_potential=p.assimilation_C*ingestion_C/(1+p.cost_growth)
    nitrogen_potential=p.zoo_CN*p.assimilation_N*ingestion_N
    production_Z=np.minimum(carbon_potential,nitrogen_potential)
    respiration_Z=p.assimilation_C*ingestion_C-production_Z+p.zoo_respiration*temp*y[m.z,0]
    excretion_N=p.assimilation_N*ingestion_N-production_Z/p.zoo_CN+p.zoo_respiration*temp*y[m.z,0]/p.zoo_CN
    losses=p.mortality*temp+starve+grazing+c.dilution
    respiration_P=p.cost_N*uptake+paid+p.cost_growth*tr[:,3]*division+starve*bc
    net_PN=uptake-p.exud_N*nr-losses*bn
    net_PC=(1-p.pdoc_fraction)*fix-p.cost_N*uptake-paid-p.cost_growth*tr[:,3]*division-p.exud_C*cr-p.reserve_N_CN*p.exud_N*nr-losses*bc
    return dict(phyto_C=bc,phyto_N=bn,carbon_reserve=cr,nitrogen_reserve=nr,
        division=division,division_specific=division/np.maximum(a,tiny),
        grazing_specific=grazing,starvation_specific=starve,
        growth_net_N=net_PN,growth_net_C=net_PC,fixation_C=fix,
        respiration_P_C=respiration_P,uptake_N=uptake,
        ingestion_C=ingestion_C,ingestion_N=ingestion_N,
        production_Z_C=production_Z,carbon_potential_Z=carbon_potential,
        nitrogen_potential_Z=nitrogen_potential,respiration_Z_C=respiration_Z,
        excretion_Z_N=excretion_N,N_limited_Z=(nitrogen_potential<carbon_potential),
        prey_flux_C=bc*grazing,prey_flux_N=bn*grazing)


def time_mean(time, values, start):
    """Trapezoidal time mean; callers must explicitly save the start endpoint."""
    selected=np.asarray(time)>=start-1e-9
    t=np.asarray(time)[selected]; v=np.asarray(values)[selected]
    if len(t)<2 or abs(t[0]-start)>1e-7:
        raise ValueError("Averaging window must start on a saved endpoint")
    return np.trapezoid(v,t,axis=0)/(t[-1]-t[0])


def size_metrics(biomass, sizes, bins=14):
    """Geometric size, ESD95 and Shannon evenness on common physical bins."""
    b=np.asarray(biomass); total=float(b.sum())
    if total<=1e-30:
        return dict(geomean=None,ESD95=None,evenness=None,pico=None,nano=None,micro=None)
    w=b/total; edges=np.geomspace(1.,200.,bins+1)
    # Treat each class as a density over its nearest log-size control volume.
    x=np.log(np.asarray(sizes))
    if len(x)==1:
        rebinned=np.histogram(sizes,bins=edges,weights=b)[0]
    else:
        old=np.r_[np.log(1.),.5*(x[:-1]+x[1:]),np.log(200.)]
        overlap=np.maximum(0.,np.minimum(old[1:,None],np.log(edges[1:]))-np.maximum(old[:-1,None],np.log(edges[:-1])))
        rebinned=(b/np.diff(old))@overlap
    q=rebinned/total; pos=q>0
    return dict(geomean=float(np.exp(w@x)),ESD95=float(sizes[min(np.searchsorted(np.cumsum(w),.95),len(w)-1)]),
        evenness=float(-np.sum(q[pos]*np.log(q[pos]))/np.log(bins)),
        pico=float(w[np.asarray(sizes)<2].sum()),nano=float(w[(np.asarray(sizes)>=2)&(np.asarray(sizes)<20)].sum()),
        micro=float(w[np.asarray(sizes)>=20].sum()))
