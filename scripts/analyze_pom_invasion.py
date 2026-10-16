"""Local rare-class growth in frozen late resident environments.

This is a linearized ecological diagnostic, not another closed ecosystem run.
Intracellular quotas relax before P invasion growth is evaluated.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dataclasses import replace
import hashlib,json
import numpy as np
from scipy.integrate import solve_ivp
from splas2 import Model
from splas2.stoichiometry import observe
from scripts.run_pom_validation import OUT,DOC,build,atomic_json,source_hashes


def diagnose(id):
    path=OUT/(id+'.npz');meta=json.loads(path.with_suffix('.json').read_text());m=build(meta['case']);p=m.physiology
    if hashlib.sha256(path.read_bytes()).hexdigest()!=meta['archive_sha256']:raise ValueError('Parent hash')
    with np.load(path) as z:y=z['state'][-1].copy()
    bc,bn=m.biomass(y);o=observe(m,y);rP=[];residuals=[];quotas=[]
    for i,size in enumerate(m.sizes):
        c=replace(m.config,n_phyto=1,n_zoo=0,phyto_min=float(size),phyto_max=float(size))
        probe=Model(c,p,reference_C=np.array([m.traits[i,1]]),reference_N=np.array([m.traits[i,2]]))
        cref,nref=probe.traits[0,1:3]
        initial=np.array([(cref-probe.traits[0,3]-p.reserve_N_CN*(nref-probe.traits[0,4]))/cref,
                          (nref-probe.traits[0,4])/nref])
        def cell(x):
            yy=np.zeros((probe.nvar,1));yy[0,0]=y[0,0];yy[1,0]=y[1,0]
            yy[probe.a,0]=1.;yy[probe.cr,0]=x[0]*cref;yy[probe.nr,0]=x[1]*nref
            return yy
        def derivative(t,x):
            yy=cell(x);f=probe.rhs(0.,yy);a=f[probe.a,0][0]
            return np.array([f[probe.cr,0][0]/cref-x[0]*a,f[probe.nr,0][0]/nref-x[1]*a])
        r=solve_ivp(derivative,[0.,200.],initial,method='Radau',rtol=1e-9,atol=1e-12)
        if not r.success or np.any(r.y[:,-1]<0):raise ValueError('Quota relaxation failed '+id)
        yy=cell(r.y[:,-1]);f=probe.rhs(0.,yy)
        rP.append(float(f[probe.a,0][0]-o['grazing_specific'][i]))
        residuals.append(float(np.max(abs(derivative(200.,r.y[:,-1])))));quotas.append(r.y[:,-1].tolist())
    temp=p.Q10**((p.temperature-18.)/10.)
    ac=m.phi.T@bc[:,0];an=m.phi.T@bn[:,0]
    k=26*m.zsizes**(-.4)*temp/(m.config.grazing_half_C+ac)
    growth=np.minimum(p.assimilation_C*k*ac/(1+p.cost_growth),p.zoo_CN*p.assimilation_N*k*an)
    mortality=p.zoo_mortality*temp*y[m.z,0].sum()/p.zoo_CN if p.zoo_closure==0 else 0.
    rZ=growth-p.zoo_respiration*temp-mortality-m.config.dilution
    rareP=bn[:,0]<1e-10*bn.sum();rareZ=y[m.z,0]<1e-10*y[m.z,0].sum()
    return dict(id=id,parent_archive_sha256=meta['archive_sha256'],environment=dict(DIN=float(y[0,0]),DIC=float(y[1,0]),PAR=m.config.surface_PAR),
        interpretation='Frozen resident environment; local rare-class potential, not proof of nonlinear/global stability.',
        quota_relaxation_days=200.,maximum_dimensionless_quota_derivative_per_day=max(residuals),
        phyto=dict(sizes=m.sizes.tolist(),rare_in_parent=rareP.tolist(),invasion_rate_per_day=rP,
            dimensionless_quota_derivative_per_day=residuals,relaxed_CR_over_Cref_NR_over_Nref=quotas,
            positive_rare_sizes=[float(s) for s,r,rare in zip(m.sizes,rP,rareP) if r>1e-4 and rare]),
        zoo=dict(sizes=m.zsizes.tolist(),rare_in_parent=rareZ.tolist(),invasion_rate_per_day=rZ.tolist(),
            positive_rare_sizes=[float(s) for s,r,rare in zip(m.zsizes,rZ,rareZ) if r>1e-4 and rare]))


def main():
    ids=['V_full_L300_N0.5','V_full_L300_N0.5_h0125','V_full_L300_N0.5_h00625']
    rows=[diagnose(id) for id in ids]
    result=dict(source_hashes=source_hashes(),diagnostic_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),rows=rows)
    atomic_json(DOC/'invasion_diagnostic.json',result)
    for r in rows:print(json.dumps(dict(id=r['id'],quota_residual=r['maximum_dimensionless_quota_derivative_per_day'],
        positive_rare_P=r['phyto']['positive_rare_sizes'],positive_rare_Z=r['zoo']['positive_rare_sizes'])),flush=True)

if __name__=='__main__':main()
