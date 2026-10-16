"""Independent ecological analysis and explicit POM decisions on all cases."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import argparse,csv,hashlib,json,zipfile,inspect
import numpy as np
from scipy.stats import spearmanr,pearsonr
from scipy.signal import find_peaks
from splas2.model import GPP,RESP_P,Z_PROD,INGEST_C,INGEST_N
from splas2.stoichiometry import observe
from splas2.pom_observation import nitrogen_network,spectrum,rebin,N_POOLS,N_EDGES
from scripts.run_pom_validation import ROOT,DOC,OUT,build,identifier,signature,atomic_json


def load(case):
    path=OUT/(case['id']+'.npz');meta=json.loads(path.with_suffix('.json').read_text())
    if signature(case)!=meta['run_signature']:raise ValueError('Source signature mismatch '+case['id'])
    if hashlib.sha256(path.read_bytes()).hexdigest()!=meta['archive_sha256']:raise IOError('Archive SHA '+case['id'])
    with zipfile.ZipFile(path) as z:
        if z.testzip() is not None:raise IOError('Archive CRC '+case['id'])
    with np.load(path) as z:data={k:z[k] for k in z.files}
    return meta,data


def average(t,v,lo,hi=None):
    if hi is None:hi=t[-1]
    select=(t>lo)&(t<hi);tt=np.r_[lo,t[select],hi]
    v=np.asarray(v);shape=v.shape[1:]
    flat=v.reshape((len(t),-1));ends=np.array([[np.interp(q,t,flat[:,i]) for i in range(flat.shape[1])] for q in [lo,hi]])
    vv=np.concatenate((ends[:1],flat[select],ends[1:]))
    result=np.trapezoid(vv,tt,axis=0)/(hi-lo)
    return result.reshape(shape) if shape else float(result[0])


def at(t,v,q,side='right'):
    hit=np.flatnonzero(abs(t-q)<1e-8)
    if len(hit):return v[hit[-1] if side=='right' else hit[0]]
    v=np.asarray(v)
    if v.ndim==1:return np.interp(q,t,v)
    return np.array([np.interp(q,t,v[:,i]) for i in range(v.shape[1])])


def ratio(x,y):return float(x/y) if y>1e-12 else None


def summarize(case):
    meta,d=load(case);m=build(case);t=d['time'];y=d['state'];p=m.physiology
    pc=m.traits[:,3][None,:]*y[:,m.a,0]+y[:,m.cr,0]+p.reserve_N_CN*y[:,m.nr,0]
    pn=m.traits[:,4][None,:]*y[:,m.a,0]+y[:,m.nr,0];zc=y[:,m.z,0]
    start=1825.;sel=t>=start;xx=t[sel]
    pcmean=average(t,pc,start);pnmean=average(t,pn,start);zmean=average(t,zc,start)
    ledger=(at(t,y[:,m.ntr:,0],t[-1])-at(t,y[:,m.ntr:,0],start))/(t[-1]-start)
    oo=[observe(m,state,ti) for state,ti in zip(y[sel],xx)]
    net=[nitrogen_network(m,state,ti) for state,ti in zip(y[sel],xx)]
    av=lambda key:average(xx,np.asarray([o[key] for o in oo]),start)
    nav=lambda key:average(xx,np.asarray([o[key] for o in net]),start)
    ledger_start=at(t,y[:,m.ntr:,0],start)
    annual=[]
    for lo in np.arange(start,3650.,365.):
        hi=lo+365.;a=average(t,pn,lo,hi);c=average(t,pc,lo,hi);z=average(t,zc,lo,hi)
        led=(at(t,y[:,m.ntr:,0],hi)-at(t,y[:,m.ntr:,0],lo))/(hi-lo)
        annual.append(dict(year=int((hi-start)/365.),PN=float(a.sum()),PC=float(c.sum()),ZC=float(z.sum()),
            CN=ratio(c.sum(),a.sum()),geomean_N=spectrum(a,m.sizes,1,200)['geomean'],
            E_P_N=spectrum(a,m.sizes,1,200)['evenness'],E_P_C=spectrum(c,m.sizes,1,200)['evenness'],
            E_Z=spectrum(z,m.zsizes,2,1000)['evenness'] if m.nZ else None,
            ZN_PN=ratio(z.sum()/p.zoo_CN,a.sum()),ZC_PC=ratio(z.sum(),c.sum()),
            production_Z_C=float(led[Z_PROD]),NPP_C=float(led[GPP]-led[RESP_P]),
            TTE_flux=ratio(led[Z_PROD],led[GPP]-led[RESP_P]),
            P_N=spectrum(a,m.sizes,1,200),P_C=spectrum(c,m.sizes,1,200),
            Z=spectrum(z,m.zsizes,2,1000) if m.nZ else None,
            diet_CN=ratio(led[INGEST_C],led[INGEST_N]),yield_Z=ratio(led[Z_PROD],led[INGEST_C])))
    pnt=pn[sel].sum(axis=1);meanPN=pnmean.sum()
    peaks,_=find_peaks(pnt,prominence=.05*meanPN);periods=np.diff(xx[peaks])
    nlimited=np.array([np.sum(o['ingestion_C']*o['N_limited_Z']) for o in oo])
    flow=nav('flows');dinout=nav('dilution');inputN=nav('external_input')
    instantaneous_balance_error=0.
    for state,ti,o in zip(y[sel][::37],xx[::37],net[::37]):
        f=m.rhs(ti,state)
        derivative=np.array([f[0,0],(m.traits[:,4]*f[m.a,0]+f[m.nr,0]).sum(),
            f[m.z,0].sum()/p.zoo_CN,f[m.dn,0].sum(),f[3,0]])
        instantaneous_balance_error=max(instantaneous_balance_error,float(abs(o['derivative']-derivative).max()))
    # Same-state counterfactual isolates immediate retention, before feedbacks.
    carbon=np.array([o['carbon_potential_Z'] for o in oo])
    ingN=np.array([o['ingestion_N'] for o in oo])
    direct_high=average(xx,np.minimum(carbon,p.zoo_CN*.67*ingN).sum(axis=1),start)
    direct_low=average(xx,np.minimum(carbon,p.zoo_CN*.335*ingN).sum(axis=1),start)
    modes_series=np.array([spectrum(a,m.sizes,1,200)['modes'] for a in pn[sel]])
    bin_series=np.array([rebin(a,m.sizes,1,200)[1] for a in pn[sel]])
    bin_fraction=bin_series/np.maximum(pnt[:,None],1e-300)
    zbasal=p.zoo_respiration*p.Q10**((p.temperature-18.)/10.)*zmean.sum()
    row=dict(id=case['id'],family=case['family'],structure=case['structure'],PAR=case['config']['surface_PAR'],
        inflow_DIN=case['config']['inflow_DIN'],zoo_mortality=p.zoo_mortality,assimilation_N=p.assimilation_N,
        N_recycling_rate=p.det_lysis_N,step=case['max_step'],nP=m.nP,nZ=m.nZ,
        PN=float(pnmean.sum()),PC=float(pcmean.sum()),ZC=float(zmean.sum()),CN=ratio(pcmean.sum(),pnmean.sum()),
        ZN_PN=ratio(zmean.sum()/p.zoo_CN,pnmean.sum()),ZC_PC=ratio(zmean.sum(),pcmean.sum()),
        production_Z_C=float(ledger[Z_PROD]),NPP_C=float(ledger[GPP]-ledger[RESP_P]),
        TTE_flux=ratio(ledger[Z_PROD],ledger[GPP]-ledger[RESP_P]),
        net_secondary_production_C=float(ledger[Z_PROD]-zbasal),
        TTE_net_secondary=ratio(ledger[Z_PROD]-zbasal,ledger[GPP]-ledger[RESP_P]),
        ingestion_C=float(ledger[INGEST_C]),ingestion_N=float(ledger[INGEST_N]),
        diet_CN=ratio(ledger[INGEST_C],ledger[INGEST_N]),yield_Z=ratio(ledger[Z_PROD],ledger[INGEST_C]),
        N_limited_ingested_fraction=ratio(average(xx,nlimited,start),float(av('ingestion_C').sum())),
        DIN=average(t,y[:,0,0],start),DON=average(t,y[:,3,0],start),detritus_N=float(average(t,y[:,m.dn,0].sum(axis=1),start)),
        carbon_reserve_fraction=ratio(average(t,y[:,m.cr,0].sum(axis=1),start),pcmean.sum()),
        nitrogen_reserve_fraction=ratio(average(t,y[:,m.nr,0].sum(axis=1),start),pnmean.sum()),
        P_N=spectrum(pnmean,m.sizes,1,200),P_C=spectrum(pcmean,m.sizes,1,200),
        Z=spectrum(zmean,m.zsizes,2,1000) if m.nZ else None,
        modes_sensitivity={str(nb):spectrum(pnmean,m.sizes,1,200,nb)['modes'] for nb in [7,14,28]},
        instantaneous_multimodal_time_fraction=average(xx,(modes_series>=2).astype(float),start),
        physical_bin_occupancy_above_1percent=average(xx,(bin_fraction>=.01).astype(float),start).tolist(),
        spectra=dict(sizes=m.sizes.tolist(),PN=pnmean.tolist(),PC=pcmean.tolist(),A=average(t,y[:,m.a,0],start).tolist(),
            CR=average(t,y[:,m.cr,0],start).tolist(),NR=average(t,y[:,m.nr,0],start).tolist(),
            division_specific=av('division_specific').tolist(),grazing_specific=av('grazing_specific').tolist(),
            uptake_N=av('uptake_N').tolist(),growth_net_N=av('growth_net_N').tolist(),
            prey_flux_N=av('prey_flux_N').tolist(),starvation_specific=av('starvation_specific').tolist(),
            zsizes=m.zsizes.tolist(),ZC=zmean.tolist()),
        PN_CV=float(np.std(pnt)/np.mean(pnt)),PN_range_05_95=np.quantile(pnt,[.05,.95]).tolist(),
        PN_peaks_5percent_prominence=len(peaks),median_interpeak_days=float(np.median(periods)) if len(periods) else None,
        half_window_PN_drift=ratio(average(t,pn.sum(axis=1),2737.5)-average(t,pn.sum(axis=1),1825.,2737.5),
                                 average(t,pn.sum(axis=1),1825.,2737.5)),
        annual_blocks=annual,N_network=dict(pools=N_POOLS,edges=[e[0] for e in N_EDGES],flows=flow.tolist(),
            inventory=nav('inventory').tolist(),dilution=dinout.tolist(),external_input=inputN,
            mean_derivative=nav('derivative').tolist(),inventory_change_per_day=((net[-1]['inventory']-net[0]['inventory'])/(xx[-1]-xx[0])).tolist(),
            instantaneous_RHS_absolute_error=instantaneous_balance_error),
        same_state_etaN=dict(production_at_067=direct_high,production_at_0335=direct_low,
            relative_change=ratio(direct_low-direct_high,direct_high)),
        run_signature=meta['run_signature'],archive_sha256=meta['archive_sha256'],
        minimum_live_state=meta['minimum_live_state'],minimum_DIC=meta['minimum_DIC'],
        residual_CN=meta['maximum_relative_residual'],runtime_seconds=meta['runtime_seconds'])
    return row


def growth(high,low):return ratio(high-low,low) if high is not None and low is not None else None


def correlations(rows):
    result={}
    for currency in ['N','C']:
        ykey='ZN_PN' if currency=='N' else 'ZC_PC'
        valid=[r for r in rows if r[ykey] and r['Z'] is not None and r['P_'+currency]['evenness'] is not None]
        if len(valid)<5:result[currency]=dict(n=len(valid),pass_signature=False);continue
        x=np.array([r['P_'+currency]['evenness'] for r in valid]);z=np.array([r['Z']['evenness'] for r in valid])
        yy=np.log10([r[ykey] for r in valid]);condition=np.linalg.cond(np.c_[np.ones(len(x)),x-x.mean(),z-z.mean(),(x-x.mean())*(z-z.mean())])
        def corr(a,b,kind):
            if np.std(a)<1e-12 or np.std(b)<1e-12:return None
            return float((spearmanr if kind=='rank' else pearsonr)(a,b).statistic)
        rp=corr(x,yy,'rank');rz=corr(z,yy,'rank')
        fit=np.linalg.lstsq(np.c_[np.ones(len(x)),x-x.mean(),z-z.mean(),(x-x.mean())*(z-z.mean())],yy,rcond=None)[0]
        result[currency]=dict(n=len(valid),rho_P=rp,rho_Z=rz,pearson_P=corr(x,yy,'linear'),pearson_Z=corr(z,yy,'linear'),
            pass_signature=bool(rp is not None and rz is not None and rp<-.2 and rz>.2),
            interaction_coefficients=fit.tolist(),interaction_design_condition_number=float(condition),
            interpretation='Deterministic cohort association; no sampling p-values or causal diversity effect.')
    return result


def targets(rows):
    by={r['id']:r for r in rows};decisions={}
    for structure in ['full','own','flat','uniform','prior','noZ']:
        contrasts=[];quality=[]
        for light in [30.,300.]:
            lo=by[identifier(structure,light,.5)];hi=by[identifier(structure,light,8.)]
            c=dict(PAR=light,PN_change=growth(hi['PN'],lo['PN']),geomean_change=growth(hi['P_N']['geomean'],lo['P_N']['geomean']),
                ESD95_change=growth(hi['P_N']['ESD95'],lo['P_N']['ESD95']),small_fraction_high=hi['P_N']['small_fraction'])
            c['pass']=bool(c['PN_change']>.05 and c['geomean_change']>.05 and c['ESD95_change']>.2 and c['small_fraction_high']>=.01)
            contrasts.append(c)
        cohort=[by[identifier(structure,l,n)] for l in [3.,30.,300.] for n in [.5,2.,8.]]
        modal=[r for r in cohort if r['PAR'] in [30.,300.] and r['P_N']['modes']>=2]
        for n in [.5,2.,8.]:
            for l0,l1 in [(3.,30.),(30.,300.)]:
                lo=by[identifier(structure,l0,n)];hi=by[identifier(structure,l1,n)]
                q=dict(N=n,light=[l0,l1],PC_change=growth(hi['PC'],lo['PC']),CN_change=growth(hi['CN'],lo['CN']),
                    diet_CN_change=growth(hi['diet_CN'],lo['diet_CN']),production_change=growth(hi['production_Z_C'],lo['production_Z_C']),
                    TTE_flux_change=growth(hi['TTE_flux'],lo['TTE_flux']))
                q['pass']=bool(all(q[k] is not None for k in ['PC_change','CN_change','diet_CN_change','production_change']) and
                    q['PC_change']>.05 and q['CN_change']>.05 and q['diet_CN_change']>.05 and q['production_change']<q['PC_change'])
                quality.append(q)
        corr=correlations(cohort)
        decisions[structure]=dict(role='mechanism removal control' if structure=='noZ' else 'candidate structure',
            P1=dict(pass_pattern=all(c['pass'] for c in contrasts),contrasts=contrasts),
            P2=dict(pass_pattern=len(modal)>=2,cases_with_multiple_modes=[r['id'] for r in modal],
                modes_all_6={r['id']:r['P_N']['modes'] for r in cohort if r['PAR'] in [30.,300.]}),
            P3=dict(pass_pattern=len({q['N'] for q in quality if q['pass']})>=2,contrasts=quality),
            P4=dict(pass_pattern=all(c.get('pass_signature',False) for c in corr.values()),common_9=corr))
        if structure in ['full','own']:
            extended=cohort+[r for r in rows if r['structure']==structure and r['family']=='closure_gradient']
            decisions[structure]['P4']['extended_15']=correlations(extended)
        decisions[structure]['all_four_patterns']=all(decisions[structure][k]['pass_pattern'] for k in ['P1','P2','P3','P4'])
    return decisions


def refinements(rows):
    by={r['id']:r for r in rows};result=[]
    for r in rows:
        if r['family'] not in ['numerical_control','closure_density_control','fine_plane','fine_corner']:continue
        base=by[identifier('own',r['PAR'],r['inflow_DIN'])] if r['family']=='closure_density_control' else by[r['id'].rsplit('_',1)[0]]
        compare={}
        for key in ['PN','PC','CN','production_Z_C','diet_CN','TTE_flux','ZC']:
            compare[key]=growth(r[key],base[key])
        compare['geomean_N']=growth(r['P_N']['geomean'],base['P_N']['geomean'])
        compare['evenness_N']=growth(r['P_N']['evenness'],base['P_N']['evenness'])
        result.append(dict(id=r['id'],base=base['id'],family=r['family'],relative_changes=compare,
            modes=[base['P_N']['modes'],r['P_N']['modes']],annual_ranges={
                k:[min(b[k] for b in base['annual_blocks']),max(b[k] for b in base['annual_blocks'])]
                for k in ['PN','CN','production_Z_C','geomean_N']},
            warning_10percent=[k for k,v in compare.items() if v is not None and abs(v)>.1]))
    return result


def filtered_observations(rows):
    result={}
    for structure in ['full','own','flat','uniform','prior']:
        selected=[];coverage=[]
        for row in rows:
            if row['family']!='holdout_plane' or row['structure']!=structure:continue
            sp=row['spectra'];x=np.asarray(sp['sizes']);zx=np.asarray(sp['zsizes'])
            pi=x>=2.;zi=zx>=200.;pc=np.asarray(sp['PC'])[pi];pn=np.asarray(sp['PN'])[pi];zc=np.asarray(sp['ZC'])[zi]
            coverage.append(dict(id=row['id'],observed_PC=float(pc.sum()),observed_ZC=float(zc.sum()),
                observed_ZC_fraction=ratio(zc.sum(),row['ZC']),included=bool(pc.sum()>1e-8 and zc.sum()>1e-8)))
            if pc.sum()<=1e-8 or zc.sum()<=1e-8:continue
            rr=dict(row);rr.update(P_C=spectrum(pc,x[pi],2,200),P_N=spectrum(pn,x[pi],2,200),
                Z=spectrum(zc,zx[zi],200,1000),ZC_PC=ratio(zc.sum(),pc.sum()),ZN_PN=ratio(zc.sum()/5.5,pn.sum()))
            selected.append(rr)
        result[structure]=dict(filter='P class centres 2-200um; Z centres 200-1000um; both C stocks >1e-8',
            status='post-primary observational diagnostic; does not replace P4',coverage=coverage,associations=correlations(selected))
    return result


def recycling(rows):
    by={r['id']:r for r in rows};result=[]
    for n in [.5,8.]:
        for rate in [.05,.15,.45]:
            def row(eta):
                id=identifier('full',30,n)
                return by[id] if eta==.67 and rate==.15 else by[id+f'_a{eta:g}_r{rate:g}']
            lo,hi=row(.335),row(.67)
            result.append(dict(N=n,N_recycling_rate=rate,
                direct_same_state_change=hi['same_state_etaN']['relative_change'],
                total_trajectory_changes={k:growth(lo[k],hi[k]) for k in ['PN','PC','CN','production_Z_C','diet_CN','ZC','DIN']},
                size_change=growth(lo['P_N']['geomean'],hi['P_N']['geomean']),
                evenness_change=growth(lo['P_N']['evenness'],hi['P_N']['evenness']),
                baseline_N_flows=hi['N_network'],lower_eta_N_flows=lo['N_network']))
    return result


def memory_analysis(cases):
    loaded={c['id']:(c,*load(c)) for c in cases};groups={};raw=[]
    for id,(c,meta,d) in loaded.items():
        m=build(c);t=d['time'];y=d['state'];p=m.physiology
        pc=m.traits[:,3][None,:]*y[:,m.a,0]+y[:,m.cr,0]+p.reserve_N_CN*y[:,m.nr,0]
        pn=m.traits[:,4][None,:]*y[:,m.a,0]+y[:,m.nr,0]
        v={}
        for day in [6.,9.,16.,36.]:
            yy=at(t,y[:,:,0],day);obs=observe(m,yy[:,None],day)
            v[f'PN_d{day:g}']=float(at(t,pn.sum(axis=1),day))
            v[f'CN_d{day:g}']=ratio(at(t,pc.sum(axis=1),day),at(t,pn.sum(axis=1),day))
            v[f'NR_fraction_d{day:g}']=ratio(yy[m.nr].sum(),at(t,pn.sum(axis=1),day))
            v[f'Zproduction_d{day:g}']=float(obs['production_Z_C'].sum())
        for lo,hi in [(0.,6.),(6.,36.)]:
            v[f'integrated_Zproduction_{lo:g}_{hi:g}']=float(at(t,y[:,m.ntr+Z_PROD,0],hi)-at(t,y[:,m.ntr+Z_PROD,0],lo,side='left'))
        key=(c['config']['inflow_DIN'],c['phase'],c['max_step'])
        groups.setdefault(key,{})[c['arm']]=dict(values=v,initial_sha=meta['common_initial_state_sha256'],seed=c['seed'])
        raw.append(dict(id=id,values=v,residual_CN=meta['maximum_relative_residual'],archive_sha256=meta['archive_sha256']))
    comparisons=[]
    for (n,phase,step),arms in groups.items():
        shas={r['initial_sha'] for r in arms.values()}
        if len(shas)!=1:raise ValueError('Memory cases did not share identical starting arrays')
        metrics={}
        for k in arms['control']['values']:
            y=lambda arm:arms[arm]['values'][k]
            nl=y('NeLl')-y('Ne')-y('Ll')+y('control')
            ln=y('LeNl')-y('Nl')-y('Le')+y('control')
            metrics[k]=dict(control=y('control'),raw_sequence_difference=y('NeLl')-y('LeNl'),
                interaction_NL=nl,interaction_LN=ln,order_interaction=nl-ln,
                normalized_order_interaction=ratio(nl-ln,abs(y('control'))))
        comparisons.append(dict(N=n,phase=phase,max_step=step,seed=arms['control']['seed'],
            exact_common_initial_sha256=next(iter(shas)),metrics=metrics))
    return dict(contrasts=comparisons,cases=raw)


def garden_analysis(cases):
    result=[]
    for c in cases:
        meta,d=load(c);m=build(c);t=d['time'];y=d['state'];obs0=observe(m,y[0]);obs1=observe(m,y[-1])
        result.append(dict(id=c['id'],N=c['config']['inflow_DIN'],phase=c['phase'],max_step=c['max_step'],arm=c['arm'],
            initial_live_CN=meta['initial_total_CN'],transplanted_CN=meta['transplanted_CN'],fallback_receiver_N_fraction=meta['fallback_receiver_N_fraction'],
            initial_division=float(obs0['division'].sum()),initial_uptake=float(obs0['uptake_N'].sum()),
            initial_Zproduction=float(obs0['production_Z_C'].sum()),
            integrated_Zproduction=float(y[-1,m.ntr+Z_PROD,0]-y[0,m.ntr+Z_PROD,0]),
            integrated_uptake=float(y[-1,m.ntr+2,0]-y[0,m.ntr+2,0]),
            final_PN=float(obs1['phyto_N'].sum()),residual_CN=meta['maximum_relative_residual'],archive_sha256=meta['archive_sha256']))
    return result


def reintroduction_analysis(cases):
    result=[]
    for c in cases:
        meta,d=load(c);m=build(c);t=d['time'];y=d['state']
        pc=m.traits[:,3][None,:]*y[:,m.a,0]+y[:,m.cr,0]+4*y[:,m.nr,0]
        pn=m.traits[:,4][None,:]*y[:,m.a,0]+y[:,m.nr,0]
        cm=average(t,pc,365.);nm=average(t,pn,365.)
        parent=c['seed']['parent_id']
        with np.load(OUT/(parent+'.npz')) as z:original=z['state'][-1]
        parentPN=m.biomass(original)[1][:,0];rare=parentPN<1e-10*parentPN.sum()
        percap=np.full(m.nP,np.nan);good=(y[0,m.a,0]>0)&(y[-1,m.a,0]>0)
        percap[good]=(np.log(y[-1,m.a,0][good])-np.log(y[0,m.a,0][good]))/730.
        result.append(dict(id=c['id'],parent=parent,redistribute_fraction=c['redistribute_PN_fraction'],
            max_step=c['max_step'],duration=730.,initial_live_CN=meta['initial_total_CN'],
            DIC_redistribution=meta['DIC_redistribution'],N_redistributed=meta['N_redistributed'],
            PN=float(nm.sum()),PC=float(cm.sum()),CN=ratio(cm.sum(),nm.sum()),P_N=spectrum(nm,m.sizes,1,200),
            final_CN=ratio(pc[-1].sum(),pn[-1].sum()),final_PN_in_classes_rare_in_parent=float(pn[-1,rare].sum()),
            spectra=dict(sizes=m.sizes.tolist(),PN_mean=nm.tolist(),PN_initial=pn[0].tolist(),PN_final=pn[-1].tolist(),
                rare_in_parent=rare.tolist(),realized_log_abundance_growth_per_day=[float(v) if np.isfinite(v) else None for v in percap]),
            residual_CN=meta['maximum_relative_residual'],archive_sha256=meta['archive_sha256']))
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--community-only',action='store_true');a=parser.parse_args()
    reg=json.loads((DOC/'community_registry.json').read_text());rows=[];allcases=list(reg['cases'])
    if (DOC/'followup_registry.json').exists():
        follow=json.loads((DOC/'followup_registry.json').read_text())
        if follow['followup_protocol_sha256']!=hashlib.sha256((DOC/'FOLLOWUP_PROTOCOL.md').read_bytes()).hexdigest():
            raise ValueError('Changed followup protocol')
        allcases+=follow['cases']
    rarecases=[]
    if (DOC/'refinement_registry.json').exists():
        refine=json.loads((DOC/'refinement_registry.json').read_text())
        if refine['refinement_protocol_sha256']!=hashlib.sha256((DOC/'REFINEMENT_PROTOCOL.md').read_bytes()).hexdigest():
            raise ValueError('Changed refinement protocol')
        allcases+=[c for c in refine['cases'] if c['family']!='rare_reintroduction']
        rarecases=[c for c in refine['cases'] if c['family']=='rare_reintroduction']
    if (DOC/'rare_pool_registry.json').exists():
        pool=json.loads((DOC/'rare_pool_registry.json').read_text())
        if pool['rare_pool_protocol_sha256']!=hashlib.sha256((DOC/'RARE_POOL_PROTOCOL.md').read_bytes()).hexdigest():
            raise ValueError('Changed rare-pool protocol')
        rarecases+=pool['cases']
    analysis_hash=hashlib.sha256(''.join(inspect.getsource(f) for f in [summarize,average,at,ratio]).encode()).hexdigest()
    cache=OUT/'analysis_cache';cache.mkdir(exist_ok=True)
    for i,c in enumerate(allcases):
        cp=cache/(c['id']+'.json');cached=json.loads(cp.read_text()) if cp.exists() else None
        if cached and cached['analysis_operator_sha256']==analysis_hash and cached['row']['run_signature']==signature(c):
            meta,_=load(c)
            if cached['row']['archive_sha256']!=meta['archive_sha256']:raise ValueError('Stale observation cache')
            rows.append(cached['row'])
        else:
            row=summarize(c);rows.append(row)
            atomic_json(cp,dict(analysis_operator_sha256=analysis_hash,row=row))
        if (i+1)%10==0:print(f'Analyzed {i+1}/{len(allcases)}',flush=True)
    result=dict(protocol_sha256=reg['protocol_sha256'],baseline_commit=reg['baseline_commit'],analysis_operator_sha256=analysis_hash,rows=rows,
        pattern_decisions=targets(rows),numerical_controls=refinements(rows),N_retention_recycling=recycling(rows),
        filtered_observation_diagnostic=filtered_observations(rows))
    result['annual_pattern_decisions']=[]
    for block in range(5):
        annual_rows=[]
        for row in rows:
            rr=dict(row);rr.update(row['annual_blocks'][block]);annual_rows.append(rr)
        result['annual_pattern_decisions'].append(dict(year=block+1,decisions=targets(annual_rows)))
    fine={identifier('full',r['PAR'],r['inflow_DIN']):r for r in rows if r['family']=='fine_plane'}
    if fine:
        if len(fine)!=9:raise ValueError('Incomplete whole-plane refinement')
        combined=[]
        for row in rows:
            if row['family']=='holdout_plane' and row['structure']=='full':
                rr=dict(fine[row['id']]);rr.update(id=row['id'],family='holdout_plane');combined.append(rr)
            elif row['family'] not in ['fine_plane','fine_corner']:combined.append(row)
        decision=targets(combined)['full'];decision['P4'].pop('extended_15',None)
        result['fine_full_pattern_decisions']=dict(step=.00625,source_case_ids={id:r['id'] for id,r in fine.items()},decision=decision,annual=[])
        for block in range(5):
            yy=[]
            for row in combined:
                rr=dict(row);rr.update(row['annual_blocks'][block]);yy.append(rr)
            dd=targets(yy)['full'];dd['P4'].pop('extended_15',None)
            result['fine_full_pattern_decisions']['annual'].append(dd)
    if not a.community_only:
        for stage,fn,key in [('memory',memory_analysis,'memory'),('garden',garden_analysis,'reserve_common_garden')]:
            cases=json.loads((DOC/(stage+'_registry.json')).read_text())['cases'];result[key]=fn(cases)
        if rarecases:result['rare_reintroduction']=reintroduction_analysis(rarecases)
    # Diagnostic associations do not replace the prespecified P4 correlations.
    result['secondary_diversity_diagnostics']={}
    for structure in ['full','own','flat','prior']:
        rr=[r for r in rows if r['family']=='holdout_plane' and r['structure']==structure]
        cs={}
        for k in ['production_Z_C','net_secondary_production_C','TTE_flux','TTE_net_secondary','PN','ZC']:
            valid=[r for r in rr if r[k] is not None and r['Z']['evenness'] is not None]
            cs[k]=dict(n=len(valid),rho=float(spearmanr([r['Z']['evenness'] for r in valid],[r[k] for r in valid]).statistic)) if len(valid)>4 else dict(n=len(valid),rho=None)
        result['secondary_diversity_diagnostics'][structure]=cs
    allmeta=[json.loads(p.read_text()) for p in sorted(OUT.glob('*.json')) if p.name.startswith(('V_','T_','G_','B_'))]
    result['verification']=dict(community_cases=len(rows),memory_cases=len(result.get('memory',{}).get('cases',[])),
        garden_cases=len(result.get('reserve_common_garden',[])),rare_reintroduction_cases=len(result.get('rare_reintroduction',[])),
        total_cases=len(allmeta),maximum_relative_residual_CN=np.max([r['maximum_relative_residual'] for r in allmeta],axis=0).tolist(),
        minimum_live_state=min(r['minimum_live_state'] for r in allmeta),minimum_DIC=min(r['minimum_DIC'] for r in allmeta),
        maximum_N_network_RHS_absolute_error=max(r['N_network']['instantaneous_RHS_absolute_error'] for r in rows),
        archive_check='SHA256 and complete ZIP CRC for every analyzed archive; signatures include frozen scientific source files.')
    atomic_json(DOC/'summary.json',result)
    atomic_json(DOC/'PATTERN_DECISIONS.json',dict(primary_registered_plane=result['pattern_decisions'],
        fine_reference_plane=result.get('fine_full_pattern_decisions'),general_ecological_validation=False,
        note='Primary results are retained; whole-plane refinement qualifies the ecological verdict.'))
    simple=[]
    for r in rows:
        d={k:v for k,v in r.items() if isinstance(v,(str,int,float)) or v is None}
        for group,prefix in [('P_N','P_N'),('P_C','P_C'),('Z','Z')]:
            if r[group]:
                for key in ['geomean','ESD95','evenness','modes','small_fraction']:d[prefix+'_'+key]=r[group][key]
        simple.append(d)
    fields=list(dict.fromkeys(k for r in simple for k in r))
    with (DOC/'summary.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');writer.writeheader();writer.writerows(simple)
    print(json.dumps(dict(verification=result['verification'],decisions={k:{p:v[p]['pass_pattern'] for p in ['P1','P2','P3','P4']} for k,v in result['pattern_decisions'].items()})),flush=True)

if __name__=='__main__':main()
