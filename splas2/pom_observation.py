"""Ecological observations for the independent 0D POM challenge.

No feedback to the integrator. N flows are mmol N m^-3 d^-1. Size spectra are
densities on fixed physical log-ESD intervals, shared across numerical grids.
"""
from __future__ import annotations
import numpy as np
from .stoichiometry import observe

N_POOLS=('DIN','P','Z','detritus','DON')
N_EDGES=(('uptake',0,1),('starvation',1,0),('grazing',1,2),
         ('P_exudation',1,4),('P_mort_detritus',1,3),('P_mort_DON',1,4),
         ('Z_egestion',2,3),('Z_excretion',2,0),('Z_mort_detritus',2,3),
         ('Z_mort_external',2,-1),('lysis',3,4),('remineralization',4,0))


def nitrogen_network(m,y,t=0.):
    """Material routing, with ingestion represented before assimilation.

    The Z node includes instantaneous passage of ingested food. Its net change
    therefore subtracts egestion/excretion as well as mortality and dilution.
    """
    o=observe(m,y,t);p=m.physiology;c=m.config
    temp=p.Q10**((p.temperature-18.)/10.)
    pn=o['phyto_N'];zn=y[m.z,0]/p.zoo_CN
    closure=np.sum(zn) if p.zoo_closure==0 else zn
    zm=float(np.sum(p.zoo_mortality*temp*closure*zn))
    mort=float(p.mortality*temp*pn.sum())
    v=np.array([o['uptake_N'].sum(),np.dot(o['starvation_specific'],pn),
        o['ingestion_N'].sum(),p.exud_N*y[m.nr,0].sum(),
        p.mortality_to_detritus*mort,(1-p.mortality_to_detritus)*mort,
        (1-p.assimilation_N)*o['ingestion_N'].sum(),o['excretion_Z_N'].sum(),
        p.mortality_recycled*zm,(1-p.mortality_recycled)*zm,
        p.det_lysis_N*temp*y[m.dn,0].sum(),p.remin_DON*temp*y[3,0]])
    inventory=np.array([y[0,0],pn.sum(),zn.sum(),y[m.dn,0].sum(),y[3,0]])
    outflow=c.dilution*inventory
    source=c.nutrient_supply
    if c.periodic_supply:source*=np.sin(np.pi*(t+c.supply_phase)/365.)**2
    source+=c.bottom_N_flux/m.dz+c.dilution*c.inflow_DIN
    derivative=-outflow;derivative[0]+=source
    for value,(_,donor,receiver) in zip(v,N_EDGES):
        derivative[donor]-=value
        if receiver>=0:derivative[receiver]+=value
    return dict(inventory=inventory,flows=v,external_input=source,
                dilution=outflow,derivative=derivative,production_Z_N=o['production_Z_C'].sum()/p.zoo_CN)


def rebin(biomass,sizes,lo,hi,bins=14):
    """Conservative overlap of log-size control volumes; no class counting."""
    b=np.asarray(biomass);x=np.log(np.asarray(sizes));target=np.linspace(np.log(lo),np.log(hi),bins+1)
    if len(x)==1:
        q=np.histogram(x,bins=target,weights=b)[0]
    else:
        old=np.r_[np.log(lo),.5*(x[:-1]+x[1:]),np.log(hi)]
        overlap=np.maximum(0.,np.minimum(old[1:,None],target[1:])-np.maximum(old[:-1,None],target[:-1]))
        q=(b/np.diff(old))@overlap
    return np.exp(.5*(target[:-1]+target[1:])),q


def spectrum(biomass,sizes,lo,hi,bins=14):
    from scipy.signal import find_peaks
    b=np.asarray(biomass);total=float(b.sum());x,q=rebin(b,sizes,lo,hi,bins)
    if total<=1e-12:
        return dict(total=total,geomean=None,ESD95=None,evenness=None,modes=0,
                    bin_centres=x.tolist(),bin_fraction=np.zeros(bins).tolist(),small_fraction=None)
    w=b/total;fractions=q/total;positive=fractions>0
    # Padding permits a true endpoint maximum, without declaring two adjacent
    # numerical bins on the same peak to be two ecological groups.
    peaks,_=find_peaks(np.r_[0.,fractions,0.],prominence=.1*fractions.max())
    peaks=peaks-1;separate=[]
    for k in peaks:
        if not separate:separate.append(int(k));continue
        j=separate[-1]
        if k-j>=2 and fractions[j:k+1].min()<=.5*min(fractions[j],fractions[k]):
            separate.append(int(k))
        elif fractions[k]>fractions[j]:separate[-1]=int(k)
    return dict(total=total,geomean=float(np.exp(w@np.log(sizes))),
        ESD95=float(sizes[min(np.searchsorted(np.cumsum(w),.95),len(w)-1)]),
        evenness=float(-np.sum(fractions[positive]*np.log(fractions[positive]))/np.log(bins)),
        modes=len(separate),mode_sizes=x[separate].tolist(),bin_centres=x.tolist(),bin_fraction=fractions.tolist(),
        small_fraction=float(w[np.asarray(sizes)<5.85].sum()))
