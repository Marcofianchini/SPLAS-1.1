"""Conservative C:N plankton physiology, trophic transfers and finite volumes.

Units: mmol element m^-3, days, metres, ESD micrometres, PAR micromol
photons m^-2 s^-1. Abundance is in 1e8 cells m^-3. All C:N ratios are molar.
Diagnostics are integrated fluxes, in mmol element m^-3, not live pools.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict, replace
from pathlib import Path
import json
import time
import numpy as np

try:
    from numba import njit
    ACCELERATED = True
except ImportError:
    ACCELERATED = False
    def njit(*args, **kwargs):
        def decorate(f): return f
        return decorate

DIN,DIC,DOC,DON=range(4)
GPP,RESP_P,UPTAKE_N,Z_PROD,INGEST_C,INGEST_N,EXPORT_C,EXPORT_N,HTL_C,HTL_N,INPUT_C,INPUT_N,OUT_C,OUT_N,STARVE_C,RESP_Z=range(16)
DIAGNOSTICS=('gpp_C','respiration_P_C','uptake_N','production_Z_C','ingestion_C','ingestion_N',
             'bottom_export_C','bottom_export_N','htl_export_C','htl_export_N','input_C','input_N',
             'external_out_C','external_out_N','starvation_respiration_C','respiration_Z_C')

@dataclass(frozen=True)
class Physiology:
    # Priors/model hypotheses unless explicitly estimated in a calibration file.
    mu_scale: float=1.0
    photo_scale: float=2.4
    nmin_fraction: float=.35
    nmin_exponent: float=-.12
    cmax_ratio: float=2.2
    cmax_exponent: float=.20
    structural_CN: float=5.0
    balanced_CN: float=6.625
    uptake_scale: float=1.0
    cost_N: float=1.5
    cost_growth: float=.20
    maintenance: float=.06
    mobilization: float=3.0
    mortality: float=.05
    exud_C: float=.02
    exud_N: float=.01
    theta_min: float=.005
    theta_max: float=.05
    acclimation_PAR: float=150.
    alpha_Chl: float=.5
    half_DIC: float=10.
    half_N_scale: float=.1
    zoo_CN: float=5.5
    assimilation_C: float=.67
    assimilation_N: float=.67
    zoo_respiration: float=.03
    zoo_mortality: float=.2
    zoo_closure: float=0.  # 0 community-dependent, 1 own quadratic
    mortality_recycled: float=1.
    det_lysis_C: float=.08
    det_lysis_N: float=.08
    remin_DOC: float=.10
    remin_DON: float=.10
    water_attenuation: float=.03
    chl_attenuation: float=.02
    Q10: float=1.6
    temperature: float=18.
    colimitation: float=0. # 0 Liebig, 1 product
    mortality_to_detritus: float=.8
    pdoc_fraction: float=.05
    reserve_N_CN: float=4. # organic N reserve carries C; 0 = inorganic storage hypothesis

    def validate(self):
        v=asdict(self)
        if not all(np.isfinite(x) for x in v.values()): raise ValueError('Nonfinite parameter')
        signed={'nmin_exponent','cmax_exponent','temperature'}
        if any(x<0 for k,x in v.items() if k not in signed): raise ValueError('Negative rate or quota')
        for k in ['assimilation_C','assimilation_N','mortality_recycled','mortality_to_detritus','pdoc_fraction']:
            if not 0<=v[k]<=1: raise ValueError(k)
        if self.theta_min>self.theta_max or not .01<=self.nmin_fraction<.9: raise ValueError('Invalid physiological bounds')
        if self.cmax_ratio<=1 or min(self.half_DIC,self.half_N_scale,self.zoo_CN,self.balanced_CN,self.Q10,self.mobilization,self.structural_CN,self.acclimation_PAR)<=0:
            raise ValueError('Invalid positive parameter')
        if self.structural_CN>=self.balanced_CN: raise ValueError('Structural C:N must leave C reserves at reference state')
        if self.reserve_N_CN>=self.structural_CN: raise ValueError('Reserve N C:N must be below structural C:N')
        if self.colimitation not in (0,1) or self.zoo_closure not in (0,1): raise ValueError('Invalid discrete formulation')

    def array(self):
        self.validate()
        return np.array(list(asdict(self).values()),dtype=float)

@dataclass(frozen=True)
class Config:
    n_phyto: int=70
    n_zoo: int=70
    n_detritus: int=2
    n_layers: int=1
    depth: float=100.
    phyto_min: float=1.
    phyto_max: float=200.
    zoo_min: float=2.
    zoo_max: float=1000.
    kernel_width: float=.25
    grazing_half_C: float=19.875
    surface_PAR: float=300.
    forcing: int=0 # constant, historical semiannual, matched annual, matched constant, 12:12
    diffusivity_m2_s: float=0.
    deep_diffusivity_m2_s: float=-1. # optional linear profile; -1 = constant
    export_bottom: bool=False
    phyto_sinking_coefficient: float=.0041
    detritus_sinking_reference: float=8.
    nutrient_supply: float=0. # uniform volumetric source
    bottom_N_flux: float=0. # mmol N m^-2 d^-1
    periodic_supply: bool=False
    dilution: float=0.
    inflow_DIN: float=0.
    inflow_DIC: float=2200.
    self_shading: bool=True
    dic_exchange: float=0. # relaxation to inflow_DIC, external C ledger
    light_phase: float=0. # days, phase of diel/seasonal light
    supply_phase: float=0. # days, independent phase of nutrient forcing

    def validate(self):
        if self.n_phyto<1 or self.n_zoo<0 or self.n_detritus<1 or self.n_layers<1: raise ValueError('Invalid dimensions')
        for k in ['depth','phyto_min','zoo_min','kernel_width','grazing_half_C']:
            if not np.isfinite(getattr(self,k)) or getattr(self,k)<=0: raise ValueError(k)
        if self.phyto_max<self.phyto_min or self.zoo_max<self.zoo_min: raise ValueError('Invalid size range')
        for k in ['surface_PAR','diffusivity_m2_s','phyto_sinking_coefficient','detritus_sinking_reference','nutrient_supply','bottom_N_flux','dilution','inflow_DIN','inflow_DIC','dic_exchange']:
            if not np.isfinite(getattr(self,k)) or getattr(self,k)<0: raise ValueError(k)
        if self.forcing not in range(5): raise ValueError('Unknown forcing')
        if self.deep_diffusivity_m2_s<0 and self.deep_diffusivity_m2_s!=-1: raise ValueError('Invalid deep diffusivity')
        if not np.isfinite(self.light_phase) or not np.isfinite(self.supply_phase): raise ValueError('Invalid phase')

def montagnes_carbon_pg(volume):
    """Montagnes et al. 1994, Table 3: 0.109 V^0.991 (NOT 0.0109)."""
    return .109*np.asarray(volume)**.991

def _log_grid(lo,hi,n):
    if n==0: return np.empty(0)
    return np.array([np.sqrt(lo*hi)]) if n==1 else np.geomspace(lo,hi,n)

def _particle_weights(size,grid):
    w=np.zeros((len(size),len(grid)))
    for i,s in enumerate(size):
        if len(grid)==1 or s<=grid[0]: w[i,0]=1.
        elif s>=grid[-1]: w[i,-1]=1.
        else:
            j=np.searchsorted(grid,s)-1
            z=(np.log(s)-np.log(grid[j]))/(np.log(grid[j+1])-np.log(grid[j]))
            w[i,j]=1-z;w[i,j+1]=z
    return w

@njit(cache=True)
def surface_light(t,base,forcing):
    if forcing==1: return base*(.5+.5*np.sin(2*np.pi*t/365.)**2)
    if forcing==2: return base*.75*(1.+np.cos(2*np.pi*(t-91.25)/365.)/3.)
    if forcing==3: return base*.75
    if forcing==4: return base if t%1.<.5 else 0.
    return base

@njit(cache=True)
def _light(y,t,tr,p,nP,nZ,nD,dz,base,forcing,shading):
    m=y.shape[1];out=np.zeros(m);chl=np.zeros(m)
    top=surface_light(t,base,forcing)
    for k in range(m):
        center=top if m==1 else top*np.exp(-p[33]*dz*.5)
        for _ in range(3):
            chlor=0.
            for i in range(nP):
                a=y[4+i,k];cr=y[4+nP+i,k];nr=y[4+2*nP+i,k]
                bc=tr[i,3]*a+cr+p[40]*nr;bn=tr[i,4]*a+nr
                status=min(1.,bn/(tr[i,2]*a+1e-300)) if a>0 else 0.
                theta=p[16]+(p[17]-p[16])*status/(1.+center/p[18])
                chlor+=12.*theta*bc
            att=p[33]+p[34]*chlor if shading else p[33]
            center=top if m==1 else top*np.exp(-att*dz*.5)
        out[k]=center;chl[k]=chlor
        top=top*np.exp(-att*dz)
    return out,chl

@njit(cache=True)
def _rhs(y,t,tr,xz,phi,pw,zw,p,phys,nP,nZ,nD,ntr,dz,base,forcing,shading,Dface,w,supply,bottom_N,dilution,inDIN,inDIC,periodic,dic_exchange,light_phase,supply_phase,enable_uptake):
    m=y.shape[1];dy=np.zeros_like(y)
    irr,chl=_light(y,t+light_phase,tr,p,nP,nZ,nD,dz,base,forcing,shading)
    z0=4+3*nP;dc0=z0+nZ;dn0=dc0+nD
    temp=p[35]**((p[36]-18.)/10.)
    for k in range(m):
        dic=y[DIC,k];din=y[DIN,k]
        for i in range(nP):
            a=y[4+i,k];cr=y[4+nP+i,k];nr=y[4+2*nP+i,k]
            c0=tr[i,3];n0=tr[i,4];cref=tr[i,1];nref=tr[i,2]
            bc=c0*a+cr+p[40]*nr;bn=n0*a+nr
            if a<=0.: continue
            rcap=max(tr[i,8]*a-c0*a-p[40]*nr,1e-300)
            cfeedback=max(0.,1.-(cr/rcap)**2)
            nfeedback=max(0.,1.-(bn/(2.*nref*a+1e-300))**2)
            nstatus=min(1.,bn/(nref*a+1e-300))
            theta=p[16]+(p[17]-p[16])*nstatus/(1.+irr[k]/p[18])
            pspecific=tr[i,6]/cref
            # The N-dependent maximum belongs INSIDE the saturation denominator.
            # Thus the low-light derivative is alpha*theta, rather than an
            # unintended second multiplication by nitrogen status.
            light=1.-np.exp(-p[19]*theta*irr[k]/max(pspecific*nstatus,1e-300))
            fix=a*tr[i,6]*temp*light*dic/(dic+p[20])*nstatus*cfeedback
            uptake=a*tr[i,7]*temp*din/(p[21]*tr[i,0]+din)*nfeedback
            if not enable_uptake:uptake=0.
            if p[9]+p[40]>0.: uptake=min(uptake,p[12]*temp*cr/(p[9]+p[40]))
            ra=p[9]*uptake
            maintenance=p[11]*temp*c0*a
            rm=min(maintenance,p[12]*temp*cr)
            deficit=max(0.,maintenance-rm)
            starve=deficit/max(bc,1e-300)
            lc=cr/max(bc,1e-300);ln=nr/max(bn,1e-300)
            limitation=min(lc,ln) if p[37]==0. else lc*ln
            j=a*tr[i,5]*temp*limitation
            j=min(j,p[12]*temp*cr/((1+p[10])*c0-p[40]*n0),p[12]*temp*nr/n0)
            rg=p[10]*c0*j
            ec=p[14]*cr;en=p[15]*nr;pd=p[39]*fix
            mort=p[13]*temp
            dy[4+i,k]+=j-(mort+starve)*a
            dy[4+nP+i,k]+=fix-pd-ra-p[40]*uptake-rm-rg-ec-(c0-p[40]*n0)*j-(mort+starve)*cr
            dy[4+2*nP+i,k]+=uptake-en-n0*j-(mort+starve)*nr
            # Starved cells pay missing maintenance by structural attrition:
            # all their C is respired, all N remineralized. Explicit hypothesis.
            dy[DIC,k]+=-fix+ra+rm+rg+starve*bc
            dy[DIN,k]+=-uptake+starve*bn
            dy[DOC,k]+=pd+ec+p[40]*en+(1-p[38])*mort*bc
            dy[DON,k]+=en+(1-p[38])*mort*bn
            for d in range(nD):
                dy[dc0+d,k]+=p[38]*mort*bc*pw[i,d]
                dy[dn0+d,k]+=p[38]*mort*bn*pw[i,d]
            dy[ntr+GPP,k]+=fix
            dy[ntr+RESP_P,k]+=ra+rm+rg+starve*bc
            dy[ntr+UPTAKE_N,k]+=uptake
            dy[ntr+STARVE_C,k]+=starve*bc

        # Carbon-based encounter; prey removed in its actual instantaneous C:N.
        losses=np.zeros(nP)
        totalZN=0.
        for j in range(nZ): totalZN+=y[z0+j,k]/p[22]
        for j in range(nZ):
            zc=y[z0+j,k]
            if zc<=0.: continue
            available=0.
            for i in range(nP):
                bc=tr[i,3]*y[4+i,k]+y[4+nP+i,k]+p[40]*y[4+2*nP+i,k]
                available+=phi[i,j]*bc
            factor=26.*xz[j]**(-.4)*temp*zc/(phys[0]+available)
            ic=0.;inn=0.
            for i in range(nP):
                bc=tr[i,3]*y[4+i,k]+y[4+nP+i,k]+p[40]*y[4+2*nP+i,k]
                bn=tr[i,4]*y[4+i,k]+y[4+2*nP+i,k]
                rate=factor*phi[i,j]
                losses[i]+=rate
                ic+=rate*bc;inn+=rate*bn
            ac=p[23]*ic;an=p[24]*inn
            growth=min(ac/(1.+p[10]),p[22]*an)
            resp=ac-growth
            exn=an-growth/p[22]
            rz=p[25]*temp*zc
            closure=totalZN if p[27]==0. else zc/p[22]
            zm=p[26]*temp*closure*zc
            dy[z0+j,k]+=growth-rz-zm
            dy[DIC,k]+=resp+rz
            dy[DIN,k]+=exn+rz/p[22]
            for d in range(nD):
                dy[dc0+d,k]+=((1-p[23])*ic+p[28]*zm)*zw[j,d]
                dy[dn0+d,k]+=((1-p[24])*inn+p[28]*zm/p[22])*zw[j,d]
            dy[ntr+Z_PROD,k]+=growth
            dy[ntr+INGEST_C,k]+=ic;dy[ntr+INGEST_N,k]+=inn
            dy[ntr+RESP_Z,k]+=resp+rz
            dy[ntr+HTL_C,k]+=(1-p[28])*zm
            dy[ntr+HTL_N,k]+=(1-p[28])*zm/p[22]
        for i in range(nP):
            dy[4+i,k]-=losses[i]*y[4+i,k]
            dy[4+nP+i,k]-=losses[i]*y[4+nP+i,k]
            dy[4+2*nP+i,k]-=losses[i]*y[4+2*nP+i,k]
        for d in range(nD):
            c=p[29]*temp*y[dc0+d,k];n=p[30]*temp*y[dn0+d,k]
            dy[dc0+d,k]-=c;dy[dn0+d,k]-=n
            dy[DOC,k]+=c;dy[DON,k]+=n
        c=p[31]*temp*y[DOC,k];n=p[32]*temp*y[DON,k]
        dy[DOC,k]-=c;dy[DIC,k]+=c;dy[DON,k]-=n;dy[DIN,k]+=n
        src=supply
        if periodic: src*=np.sin(np.pi*(t+supply_phase)/365.)**2
        if k==m-1: src+=bottom_N/dz
        dy[DIN,k]+=src;dy[ntr+INPUT_N,k]+=src
        carbon_exchange=dic_exchange*(inDIC-dic)
        dy[DIC,k]+=carbon_exchange
        dy[ntr+INPUT_C,k]+=max(0.,carbon_exchange)
        dy[ntr+OUT_C,k]+=max(0.,-carbon_exchange)
        if dilution>0.:
            cexport=0.;nexport=0.
            for r in range(ntr):
                loss=dilution*y[r,k];dy[r,k]-=loss
                if r==DIC or r==DOC or dc0<=r<dn0: cexport+=loss
                if r==DIN or r==DON or dn0<=r<ntr: nexport+=loss
                if 4<=r<4+nP:
                    cexport+=tr[r-4,3]*loss;nexport+=tr[r-4,4]*loss
                elif 4+nP<=r<4+2*nP: cexport+=loss
                elif 4+2*nP<=r<z0: nexport+=loss;cexport+=p[40]*loss
                elif z0<=r<dc0: cexport+=loss;nexport+=loss/p[22]
            dy[DIN,k]+=dilution*inDIN;dy[DIC,k]+=dilution*inDIC
            dy[ntr+INPUT_C,k]+=dilution*inDIC;dy[ntr+INPUT_N,k]+=dilution*inDIN
            dy[ntr+OUT_C,k]+=cexport;dy[ntr+OUT_N,k]+=nexport
    # Finite-volume face fluxes. All inventories and abundance are transported.
    for r in range(ntr):
        for k in range(m-1):
            flux=-Dface[k]*(y[r,k+1]-y[r,k])/dz+w[r]*y[r,k]
            dy[r,k]-=flux/dz;dy[r,k+1]+=flux/dz
        if w[r]>0.:
            flux=w[r]*y[r,m-1]/dz
            dy[r,m-1]-=flux
            if 4<=r<4+nP:
                dy[ntr+EXPORT_C,m-1]+=tr[r-4,3]*flux
                dy[ntr+EXPORT_N,m-1]+=tr[r-4,4]*flux
            elif 4+nP<=r<4+2*nP: dy[ntr+EXPORT_C,m-1]+=flux
            elif 4+2*nP<=r<z0:
                dy[ntr+EXPORT_N,m-1]+=flux;dy[ntr+EXPORT_C,m-1]+=p[40]*flux
            elif dc0<=r<dn0: dy[ntr+EXPORT_C,m-1]+=flux
            elif dn0<=r<ntr: dy[ntr+EXPORT_N,m-1]+=flux
    return dy

class Model:
    def __init__(self,config=None,physiology=None,*,reference_C=None,reference_N=None,volumes=None,mu_cap=None):
        self.config=Config() if config is None else config
        self.physiology=Physiology() if physiology is None else physiology
        c=self.config;p=self.physiology
        c.validate();p.validate()
        self.p=p.array();self.nP=c.n_phyto;self.nZ=c.n_zoo;self.nD=c.n_detritus
        self.a=slice(4,4+self.nP);self.cr=slice(self.a.stop,self.a.stop+self.nP)
        self.nr=slice(self.cr.stop,self.cr.stop+self.nP);self.z=slice(self.nr.stop,self.nr.stop+self.nZ)
        self.dc=slice(self.z.stop,self.z.stop+self.nD);self.dn=slice(self.dc.stop,self.dc.stop+self.nD)
        self.ntr=self.dn.stop;self.nvar=self.ntr+16;self.dz=c.depth/c.n_layers
        self.depth=(np.arange(c.n_layers)+.5)*self.dz
        self.sizes=_log_grid(c.phyto_min,c.phyto_max,self.nP)
        self.zsizes=_log_grid(c.zoo_min,c.zoo_max,self.nZ)
        vol=np.pi/6*self.sizes**3 if volumes is None else np.asarray(volumes,dtype=float)
        if vol.shape!=(self.nP,) or np.any(vol<=0): raise ValueError('Invalid cell volumes')
        self.sizes=(6*vol/np.pi)**(1/3)
        cref=montagnes_carbon_pg(vol)/120. if reference_C is None else np.asarray(reference_C,dtype=float)
        nref=cref/p.balanced_CN if reference_N is None else np.asarray(reference_N,dtype=float)
        if cref.shape!=(self.nP,) or nref.shape!=(self.nP,) or min(cref.min(),nref.min())<=0: raise ValueError('Invalid reference quotas')
        nfrac=np.clip(p.nmin_fraction*(self.sizes/5.)**p.nmin_exponent,.05,.8)
        n0=nref*nfrac
        # Bound structural carbon below measured/reference carbon, explicitly.
        c0=np.minimum(p.structural_CN*n0,.85*cref)
        if np.any(c0<=p.reserve_N_CN*n0) or np.any(cref<c0+p.reserve_N_CN*(nref-n0)):
            raise ValueError('Reference quotas incompatible with organic N reserves')
        mcap=2.6*self.sizes**(-.45)*p.mu_scale if mu_cap is None else np.asarray(mu_cap,dtype=float)
        cmax=cref*np.clip(p.cmax_ratio*(self.sizes/5.)**p.cmax_exponent,1.05,8.)
        self.traits=np.column_stack((self.sizes,cref,nref,c0,n0,mcap,p.photo_scale*mcap*cref,9.1e-4*vol**.67*p.uptake_scale,cmax))
        opt=.16*self.zsizes*np.exp(-.02*np.log(self.zsizes)**2)
        self.phi=np.exp(-((np.log10(self.sizes[:,None])-np.log10(opt[None,:]))/c.kernel_width)**2)
        self.particle_sizes=_log_grid(2.,200.,self.nD)
        self.pw=_particle_weights(self.sizes,self.particle_sizes)
        self.zw=_particle_weights(np.maximum(2.,.1*self.zsizes),self.particle_sizes)
        deep=c.diffusivity_m2_s if c.deep_diffusivity_m2_s<0 else c.deep_diffusivity_m2_s
        self.Dface=np.linspace(c.diffusivity_m2_s,deep,max(1,c.n_layers-1))*86400.
        self.w=np.zeros(self.ntr)
        # Transport can have internal sinking with reflecting bottom? This version
        # enables sinking only with an explicitly exporting bottom.
        if c.export_bottom:
            wp=c.phyto_sinking_coefficient*self.sizes**1.17
            self.w[self.a]=wp;self.w[self.cr]=wp;self.w[self.nr]=wp
            self.w[self.dc]=c.detritus_sinking_reference*self.particle_sizes/50.
            self.w[self.dn]=self.w[self.dc]
        self.wC=np.zeros(self.ntr);self.wN=np.zeros(self.ntr)
        self.wC[[DIC,DOC]]=1;self.wN[[DIN,DON]]=1
        self.wC[self.a]=c0;self.wN[self.a]=n0
        self.wC[self.cr]=1;self.wN[self.nr]=1;self.wC[self.nr]=p.reserve_N_CN
        self.wC[self.z]=1;self.wN[self.z]=1/p.zoo_CN
        self.wC[self.dc]=1;self.wN[self.dn]=1
        self.names=['DIN','DIC','DOC','DON']+[f'a_{i}' for i in range(self.nP)]+[f'Cres_{i}' for i in range(self.nP)]+[f'Nres_{i}' for i in range(self.nP)]+[f'ZC_{i}' for i in range(self.nZ)]+[f'DC_{i}' for i in range(self.nD)]+[f'DN_{i}' for i in range(self.nD)]+list(DIAGNOSTICS)

    def initial(self,DIN_value=1.,DIC_value=2200.,phyto_N=.7,zoo_N=.07,gradient=False):
        c=self.config
        y=np.zeros((self.nvar,c.n_layers))
        y[DIN]=DIN_value*(np.arange(c.n_layers)+1) if gradient else DIN_value
        y[DIC]=DIC_value
        a=phyto_N/self.nP/self.traits[:,2]
        y[self.a]=a[:,None]
        y[self.cr]=((self.traits[:,1]-self.traits[:,3]-self.physiology.reserve_N_CN*(self.traits[:,2]-self.traits[:,4]))*a)[:,None]
        y[self.nr]=((self.traits[:,2]-self.traits[:,4])*a)[:,None]
        if self.nZ: y[self.z]=zoo_N*self.physiology.zoo_CN/self.nZ
        return y

    def rhs(self,t,y):
        c=self.config
        return _rhs(np.asarray(y),float(t),self.traits,self.zsizes,self.phi,self.pw,self.zw,self.p,np.array([c.grazing_half_C]),self.nP,self.nZ,self.nD,self.ntr,self.dz,c.surface_PAR,c.forcing,c.self_shading,self.Dface,self.w,c.nutrient_supply,c.bottom_N_flux,c.dilution,c.inflow_DIN,c.inflow_DIC,c.periodic_supply,c.dic_exchange,c.light_phase,c.supply_phase,True)

    def light(self,t,y):
        c=self.config
        return _light(y,float(t)+c.light_phase,self.traits,self.p,self.nP,self.nZ,self.nD,self.dz,c.surface_PAR,c.forcing,c.self_shading)

    def budgets(self,y,adjusted=False):
        b=np.array([np.sum(self.wC@y[:self.ntr])*self.dz,np.sum(self.wN@y[:self.ntr])*self.dz])
        if adjusted:
            d=np.sum(y[self.ntr:],axis=1)*self.dz
            b+=np.array([d[EXPORT_C]+d[HTL_C]+d[OUT_C]-d[INPUT_C],d[EXPORT_N]+d[HTL_N]+d[OUT_N]-d[INPUT_N]])
        return b

    def biomass(self,y):
        return self.traits[:,3,None]*y[self.a]+y[self.cr]+self.physiology.reserve_N_CN*y[self.nr], self.traits[:,4,None]*y[self.a]+y[self.nr]

@dataclass
class Result:
    time: np.ndarray
    state: np.ndarray
    model: Model
    steps: int
    rejected_steps: int
    runtime_seconds: float
    rtol: float=2e-4
    atol: float=1e-8
    method: str='explicit'

    def diagnostics(self):
        m=self.model;y=self.state
        bc=m.traits[:,3][None,:,None]*y[:,m.a]+y[:,m.cr]+m.physiology.reserve_N_CN*y[:,m.nr]
        bn=m.traits[:,4][None,:,None]*y[:,m.a]+y[:,m.nr]
        C=bc.sum(axis=(1,2))*m.dz;N=bn.sum(axis=(1,2))*m.dz
        adj=np.array([m.budgets(s,True) for s in y])
        err=adj-adj[0]
        return {'steps':self.steps,'rejected_steps':self.rejected_steps,'runtime_seconds':self.runtime_seconds,
          'minimum_live_state':float(y[:,:m.ntr].min()),'max_abs_mass_error_C':float(abs(err[:,0]).max()),
          'max_abs_mass_error_N':float(abs(err[:,1]).max()),'max_relative_mass_error_C':float(abs(err[:,0]).max()/max(1.,abs(adj[0,0]))),
          'max_relative_mass_error_N':float(abs(err[:,1]).max()/max(1.,abs(adj[0,1]))),
          'final_live_CN':float(C[-1]/N[-1]) if N[-1]>0 else None,
          'accelerated':ACCELERATED,'rtol':self.rtol,'atol':self.atol,'method':self.method}

    def save(self,path):
        path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(path,time=self.time,state=self.state,sizes=self.model.sizes,zsizes=self.model.zsizes,
                            particle_sizes=self.model.particle_sizes,depth=self.model.depth,traits=self.model.traits)
        path.with_suffix('.json').write_text(json.dumps({'config':asdict(self.model.config),'physiology':asdict(self.model.physiology),**self.diagnostics()},indent=2))

@njit(cache=True)
def _implicit_uptake(y,h,args):
    """Backward-Euler coupled DIN/organic-N assimilation; scalar monotone root.

    All donor and receiver changes use the same uptake flux. Newton is
    bracketed; the lower bracket guarantees no DIN overdraw without clipping.
    Intracellular luxury feedback and energetic availability are implicit too.
    """
    tr=args[0];p=args[5];nP=args[7];ntr=args[10]
    temp=p[35]**((p[36]-18.)/10.);energy=p[9]+p[40]
    for k in range(y.shape[1]):
        old=y[DIN,k]
        if old<=0.:continue
        cap=np.zeros(nP);bn=np.zeros(nP);q=np.ones(nP);rmax=np.zeros(nP);ks=np.ones(nP);abund=np.zeros(nP)
        slope=0.
        for i in range(nP):
            a=y[4+i,k];cr=y[4+nP+i,k];nr=y[4+2*nP+i,k]
            if a<=0.:continue
            # Per-cell algebra avoids squaring subnormal abundances at extinction.
            abund[i]=a;bn[i]=tr[i,4]+nr/a;q[i]=2.*tr[i,2]
            rmax[i]=tr[i,7]*temp;ks[i]=p[21]*tr[i,0]
            cap[i]=p[12]*temp*(cr/a)/(energy*(1.+p[12]*temp*h)) if energy>0 else 1e300
            slope+=a*rmax[i]/ks[i]*max(0.,1.-(bn[i]/q[i])**2)
        lo=0.;hi=old;x=old/(1.+h*slope)
        for iteration in range(50):
            total=0.;derivative=1.
            for i in range(nP):
                feedback=max(0.,1.-(bn[i]/q[i])**2)
                rr=rmax[i]*x/(ks[i]+x)
                aa=rr*h*h/(q[i]*q[i]);bb=1.+2.*rr*h*bn[i]/(q[i]*q[i])
                u=2.*rr*feedback/(bb+np.sqrt(bb*bb+4.*aa*rr*feedback))
                du=(rmax[i]*ks[i]/(ks[i]+x)**2)*max(0.,1.-((bn[i]+h*u)/q[i])**2)/(1.+2.*rr*h*(bn[i]+h*u)/(q[i]*q[i]))
                if u>=cap[i]:u=cap[i];du=0.
                total+=abund[i]*u;derivative+=h*abund[i]*du
            value=x+h*total-old
            if value<=0.:lo=x
            else:hi=x
            tolerance=2e-13*max(1.,old)
            if abs(value)<=tolerance:
                # Root-residual stopping avoids repeated bisection after Newton
                # has converged in floating arithmetic. Stay below the root.
                lo=max(lo,x-2.*tolerance/derivative)
                break
            if hi-lo<=2e-13*max(1.,old):break
            trial=x-value/derivative
            if trial<=lo or trial>=hi:trial=.5*(lo+hi)
            x=trial
        # Stay just below the root; fluxes conserve elements without clipping
        # the DIN state when summation order differs by an ulp.
        lo*=1.-1e-12
        total=0.
        for i in range(nP):
            rr=rmax[i]*lo/(ks[i]+lo);feedback=max(0.,1.-(bn[i]/q[i])**2)
            aa=rr*h*h/(q[i]*q[i]);bb=1.+2.*rr*h*bn[i]/(q[i]*q[i])
            u=min(cap[i],2.*rr*feedback/(bb+np.sqrt(bb*bb+4.*aa*rr*feedback)))
            amount=h*u*abund[i];total+=amount
            y[4+nP+i,k]-=energy*amount;y[4+2*nP+i,k]+=amount
            y[DIC,k]+=p[9]*amount
            y[ntr+RESP_P,k]+=p[9]*amount;y[ntr+UPTAKE_N,k]+=amount
        y[DIN,k]-=total
    return y

@njit(cache=True)
def _advance(y,t,end,dt_max,args,rtol,atol,imex):
    """SSPRK3 reference, or positive first-order IMEX. Never state clipping."""
    steps=0;reject=0;proposal=dt_max
    transport_rate=2.*np.max(args[15])/(args[11]*args[11])+np.max(args[16])/args[11]+args[19]+args[23]
    transport_bound=.8/transport_rate if transport_rate>0 else dt_max
    while t<end-1e-12:
        h=min(proposal,end-t,transport_bound)
        # Exact alignment with half-day discontinuities of the 12:12 forcing.
        if args[13]==4:
            boundary=(np.floor(2*(t+args[24])+1e-10)+1.)/2.-args[24]
            h=min(h,boundary-t)
        f=_rhs(y,t,*args)
        for r in range(args[10]):
            for k in range(y.shape[1]):
                if f[r,k]<0.: h=min(h,.85*y[r,k]/(-f[r,k]))
        if h<1e-12: raise RuntimeError('Positivity step underflow; inspect kinetics')
        accepted=False
        for _ in range(50):
            yp=y+h*f
            # Evaluate the left limit at discontinuous forcing endpoints.
            f2=_rhs(yp,t+h-1e-12,*args)
            ee=yp+h*f2
            good=True
            for r in range(args[10]):
                for k in range(y.shape[1]):
                    if ee[r,k]<0. or not np.isfinite(ee[r,k]): good=False
            err=1e300
            if good:
                low=.5*y+.5*ee
                stage=.75*y+.25*ee
                third=stage+h*_rhs(stage,t+.5*h,*args)
                for r in range(args[10]):
                    for k in range(y.shape[1]):
                        if third[r,k]<0. or not np.isfinite(third[r,k]):good=False
                if good:
                    high=y/3.+2.*third/3.
                    err=0.
                    for r in range(args[10]):
                        scale=1.
                        if 4<=r<4+args[7]:scale=1./max(args[0][r-4,3],args[0][r-4,4])
                        for k in range(y.shape[1]):
                            e=abs(high[r,k]-low[r,k])/(atol*scale+rtol*max(abs(high[r,k]),abs(y[r,k])))
                            err=max(err,e)
                    if err<=1. or imex:
                        y=_implicit_uptake(high,h,args) if imex else high
                        accepted=True
                        proposal=dt_max if imex else min(dt_max,h*min(2.,.9*max(err,1e-12)**(-1./3.)))
                        break
            h*=max(.1,min(.5,.8*max(err,1e-12)**(-1./3.))) if good else .5
            reject+=1
        if not accepted: raise RuntimeError('Unable to take positive, accurate SSPRK3 step')
        t+=h;steps+=1
    return y,t,steps,reject

def simulate(model,tmax=3652.,*,y0=None,max_step=.05,save_interval=10.,times=None,rtol=2e-4,atol=1e-8,method='explicit'):
    if min(tmax,max_step,save_interval,rtol,atol)<=0: raise ValueError('Invalid time grid/tolerances')
    if method not in ('explicit','imex'):raise ValueError('Unknown integration method')
    c=model.config
    y=model.initial() if y0 is None else np.array(y0,dtype=float,copy=True)
    if y.shape!=(model.nvar,c.n_layers) or not np.isfinite(y).all() or np.any(y[:model.ntr]<0): raise ValueError('Invalid initial state')
    if times is None:
        times=np.unique(np.r_[np.arange(0,tmax,save_interval),tmax])
    else:
        times=np.asarray(times,dtype=float)
        if times[0]!=0 or times[-1]>tmax or np.any(np.diff(times)<=0): raise ValueError('Invalid output times')
    # _rhs positional argument tuple after y,t. Indices are tested explicitly.
    args=(model.traits,model.zsizes,model.phi,model.pw,model.zw,model.p,np.array([c.grazing_half_C]),
          model.nP,model.nZ,model.nD,model.ntr,model.dz,c.surface_PAR,c.forcing,c.self_shading,
          model.Dface,model.w,c.nutrient_supply,c.bottom_N_flux,c.dilution,c.inflow_DIN,c.inflow_DIC,c.periodic_supply,c.dic_exchange,c.light_phase,c.supply_phase,method=='explicit')
    out=np.empty((len(times),*y.shape));out[0]=y;t=0.;steps=rejects=0;start=time.perf_counter()
    for j,end in enumerate(times[1:],1):
        y,t,st,rj=_advance(y,t,float(end),max_step,args,rtol,atol,method=='imex')
        out[j]=y;steps+=st;rejects+=rj
    return Result(np.array(times),out,model,steps,rejects,time.perf_counter()-start,rtol,atol,method)
