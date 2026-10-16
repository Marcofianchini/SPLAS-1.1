"""Audit and analyze every registered stoichiometry run, retaining counterexamples."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import csv,hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from splas2.model import GPP,RESP_P,Z_PROD,INGEST_C,INGEST_N
from splas2.stoichiometry import size_metrics,time_mean
from scripts.run_stoichiometry import registry,build,signature,OUT,ROOT
from scripts.run_stoichiometry_controls import registry as control_registry
from scripts.run_stoichiometry_robustness import registry as robustness_registry

DOC=ROOT/'research/stoichiometry'
PALETTE=['#167d9a','#df8e36','#8052a0','#4d985b','#bb5766']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
    'axes.prop_cycle':plt.cycler(color=PALETTE),'figure.dpi':150,'savefig.dpi':190,
    'axes.titleweight':'bold','font.family':'DejaVu Sans'})

def ratio(a,b,threshold=1e-14):
    return float(a/b) if b>threshold else None

def load(id):
    with np.load(OUT/(id+'.npz')) as z:return {k:z[k] for k in z.files}

def extract(case):
    id=case['id'];meta=json.loads((OUT/(id+'.json')).read_text());path=OUT/(id+'.npz')
    if signature(case)!=meta['run_signature']:raise ValueError('Stale signature '+id)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=meta['archive_sha256']:raise ValueError('Corrupt archive '+id)
    m=build(case);d=load(id);t=d['time'];start=t[-1]-min(365.,t[-1]/2) if not case['segments'] else 730.
    mean=lambda v:time_mean(t,v,start)
    pc=mean(d['phyto_C']);pn=mean(d['phyto_N']);zc=float(mean(d['state'][:,m.z,0].sum(axis=1)))
    first=np.flatnonzero(t>=start-1e-9)[0];ledger=(d['state'][-1,m.ntr:,0]-d['state'][first,m.ntr:,0])/(t[-1]-start)
    npp=float(ledger[GPP]-ledger[RESP_P]);ic=float(ledger[INGEST_C]);inn=float(ledger[INGEST_N]);prod=float(ledger[Z_PROD])
    avg=lambda k:float(mean(d[k].sum(axis=1)))
    nlimited=float(mean((d['ingestion_C']*d['N_limited_Z']).sum(axis=1)))
    nmean=pc.sum()/pn.sum() if pn.sum()>1e-14 else None
    pnseries=d['phyto_N'].sum(axis=1);zcseries=d['state'][:,m.z,0].sum(axis=1)
    mask=t>=start;seg=pnseries[mask]
    half=start+(t[-1]-start)/2
    firsthalf=(t>=start)&(t<half);secondhalf=(t>=half)
    # Means of the two physical half windows use interpolation at the split.
    def interval_mean(v,lo,hi):
        tt=np.r_[lo,t[(t>lo)&(t<hi)],hi]
        vv=np.interp(tt,t,v)
        return float(np.trapezoid(vv,tt)/(hi-lo))
    early=interval_mean(pnseries,start,half);late=interval_mean(pnseries,half,t[-1])
    long_window=None
    if case['days']>=3650:
        lo=float(t[np.flatnonzero(t>=t[-1]-1825)[0]])
        long_mean=lambda v:time_mean(t,v,lo)
        lp=long_mean(d['phyto_N']);lc=long_mean(d['phyto_C'])
        ix=np.flatnonzero(t>=lo);peaks,_=find_peaks(pnseries[ix],prominence=.05*lp.sum())
        spacing=np.diff(t[ix[peaks]])
        long_window=dict(days=[lo,float(t[-1])],phyto_N=float(lp.sum()),CN=float(lc.sum()/lp.sum()),
            size_N=size_metrics(lp,d['sizes']),late_peaks_5percent_prominence=int(len(peaks)),
            median_interpeak_days=float(np.median(spacing)) if len(spacing)>=2 else None,
            phyto_N_CV=float(np.std(pnseries[ix])/np.mean(pnseries[ix])))
    return dict(id=id,family=case['family'],PAR=case['config']['surface_PAR'],inflow_DIN=case['config']['inflow_DIN'],
        assimilation_N=case['physiology']['assimilation_N'],nP=m.nP,nZ=m.nZ,
        phyto_C=float(pc.sum()),phyto_N=float(pn.sum()),zoo_C=zc,CN=nmean,
        ZC_PC=ratio(zc,pc.sum()),ZN_PN=ratio(zc/m.physiology.zoo_CN,pn.sum()),
        NPP_C=npp,production_Z_C=prod,ingestion_C=ic,ingestion_N=inn,diet_CN=ratio(ic,inn),
        yield_Z_per_ingested_C=ratio(prod,ic),TTE_flux=ratio(prod,npp),
        N_limited_ingested_C_fraction=ratio(nlimited,avg('ingestion_C')),
        C_reserve_fraction=ratio(avg('carbon_reserve'),pc.sum()),N_reserve_fraction=ratio(avg('nitrogen_reserve'),pn.sum()),
        DIN=float(mean(d['state'][:,0,0])),minimum_DIC=meta['minimum_DIC'],
        phyto_N_CV=float(np.std(seg)/np.mean(seg)) if np.mean(seg)>1e-14 else None,
        phyto_N_halfwindow_drift=ratio(late-early,early),
        minimum_live_state=meta['minimum_live_state'],residual_CN=meta['maximum_relative_residual'],
        size_N=size_metrics(pn,d['sizes']),size_C=size_metrics(pc,d['sizes']),
        spectra=dict(sizes=d['sizes'].tolist(),PN=pn.tolist(),PC=pc.tolist()),long_window=long_window,
        final_window_days=[float(start),float(t[-1])],runtime_seconds=meta['runtime_seconds'],
        run_signature=meta['run_signature'],archive_sha256=meta['archive_sha256'])

def savefig(fig,name):
    fig.savefig(OUT/(name+'.png'),bbox_inches='tight',facecolor='white')
    fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight',facecolor='white');plt.close(fig)

def figures(rows):
    by={r['id']:r for r in rows};lights=[15,45,135,405,1215];ns=[.25,1,4,16]
    fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    for ax,(key,title) in zip(axes.flat,[('CN','C:N del fitoplancton'),('yield_Z_per_ingested_C','Resa: produzione Z / C ingerito'),('zoo_C','Biomassa Z [mmol C m$^{-3}$]'),('phyto_C','Biomassa P [mmol C m$^{-3}$]')]):
        a=np.array([[by[f'M_L{l}_N{n:g}_a0.67'][key] or 0 for n in ns] for l in lights])
        im=ax.imshow(a,origin='lower',aspect='auto',cmap='viridis')
        ax.set(xticks=range(4),xticklabels=ns,yticks=range(5),yticklabels=lights,xlabel='DIN entrante [mmol N m$^{-3}$]',ylabel='PAR [µmol fotoni m$^{-2}$ s$^{-1}$]',title=title)
        for i in range(5):
            for j in range(4):ax.text(j,i,f'{a[i,j]:.3g}',ha='center',va='center',color='white' if a[i,j]<.65*a.max() else '#122b35',fontsize=9)
        fig.colorbar(im,ax=ax,shrink=.8)
    savefig(fig,'01_minimal_resource_plane')
    fig,axes=plt.subplots(1,3,figsize=(11,3.5),layout='constrained')
    for eta in [.335,.67,1.]:
        rr=[by[f'M_L{l}_N1_a{eta:g}'] for l in lights]
        for ax,key in zip(axes,['CN','production_Z_C','yield_Z_per_ingested_C']):
            ax.plot(lights,[r[key] for r in rr],'o-',label=f'Assim. N={eta:g}')
    for ax,title in zip(axes,['C:N delle prede','Produzione Z [mmol C m$^{-3}$ d$^{-1}$]','Resa del C ingerito']):ax.set(xscale='log',xlabel='PAR',title=title)
    axes[-1].axhline(.33,color='grey',ls=':',lw=1);axes[0].legend(fontsize=8)
    savefig(fig,'02_quantity_quality_feedback')
    fig,axes=plt.subplots(2,3,figsize=(11,6.5),layout='constrained')
    keys=['CN','phyto_N','zoo_C','size_N.geomean','size_N.evenness','yield_Z_per_ingested_C']
    titles=['C:N medio di comunità','P [mmol N m$^{-3}$]','Z [mmol C m$^{-3}$]','Taglia geometrica in N [µm]','Evenness, 14 intervalli comuni','Resa: produzione Z / C ingerito']
    for light,color in zip([45,135,405],PALETTE):
        rr=[by[f'C_L{light}_N{n}_a0.67'] for n in [1,4,16]]
        low=[by[f'C_L{light}_N{n}_a0.335'] for n in [1,4,16]]
        for ax,key in zip(axes.flat,keys):
            get=lambda r:r[key] if '.' not in key else r[key.split('.')[0]][key.split('.')[1]]
            ax.plot([1,4,16],[get(r) for r in rr],'o-',color=color,label=f'PAR{light}')
            ax.plot([1,4,16],[get(r) for r in low],'x--',color=color,alpha=.6)
    for ax,title in zip(axes.flat,titles):ax.set(xscale='log',xticks=[1,4,16],xticklabels=[1,4,16],xlabel='DIN entrante',title=title)
    axes[0,0].legend(fontsize=8);fig.supxlabel('Linee continue: assimilazione N .67; tratteggiate: .335',fontsize=9)
    savefig(fig,'03_community_patterns')
    fig,axes=plt.subplots(2,3,figsize=(11,6.5),layout='constrained')
    for row,light in enumerate([45,405]):
        ids=[f'A_L{light}_Z{z}' for z in [0,1,2,4,8,16]]
        for k,id in enumerate(ids):
            r=by[id];pn=np.array(r['spectra']['PN']);pc=np.array(r['spectra']['PC']);x=r['spectra']['sizes']
            axes[row,0].plot(x,pn/pn.sum(),label=f'{r["nZ"]} Z',color=plt.cm.viridis(k/5))
            axes[row,1].plot(x,np.divide(pc,pn,out=np.full_like(pc,np.nan),where=pn>1e-8),color=plt.cm.viridis(k/5))
        axes[row,0].set(xscale='log',xlabel='ESD P [µm]',ylabel='Frazione di PN per classe',title=f'Assemblaggio, PAR{light}')
        axes[row,1].set(xscale='log',xlabel='ESD P [µm]',ylabel='C:N della classe',title='Solo classi PN > 10$^{-8}$')
        axes[row,2].plot([by[id]['nZ'] for id in ids],[by[id]['size_N']['geomean'] for id in ids],'o-',label='Taglia geo. in N')
        axes[row,2].set(xlabel='Numero predatori introdotti',ylabel='ESD [µm]',title='Rilascio dalle taglie piccole')
    axes[0,0].legend(fontsize=7,ncols=2)
    savefig(fig,'04_matched_assembly')
    fig,axes=plt.subplots(2,2,figsize=(10,6.2),layout='constrained')
    for ax,(light,nutrient) in zip(axes.flat,[(45,1),(405,1),(135,4),(405,16)]):
        for suffix,label in [('a0.67','Selettivo, aN=.67'),('a0.335','Selettivo, aN=.335'),('uniform','Preferenza uniforme'),('noZ','Senza Z')]:
            r=by[f'C_L{light}_N{nutrient}_{suffix}'];pn=np.array(r['spectra']['PN'])
            ax.plot(r['spectra']['sizes'],pn/pn.sum(),label=label)
        ax.set(xscale='log',xlabel='ESD P [µm]',ylabel='Frazione di PN per classe',title=f'PAR{light}, DIN entrante{nutrient}')
    axes[0,0].legend(fontsize=7)
    savefig(fig,'05_foodweb_controls')
    fig,axes=plt.subplots(2,2,figsize=(10,6.3),layout='constrained')
    for name,color in zip(['N_then_light','light_then_N','N_only','light_only','control'],PALETTE):
        d=load('T_'+name);sel=d['time']>=730;t=d['time'][sel]-730
        pc=d['phyto_C'][sel].sum(axis=1);pn=d['phyto_N'][sel].sum(axis=1)
        for ax,v in zip(axes.flat,[pc/pn,d['nitrogen_reserve'][sel].sum(axis=1)/pn,
            d['production_Z_C'][sel].sum(axis=1),pn]):ax.plot(t,v,label=name,color=color,lw=1.4)
    for ax,title in zip(axes.flat,['C:N delle prede','N di riserva / PN','Produzione Z [mmol C m$^{-3}$ d$^{-1}$]','PN [mmol N m$^{-3}$]']):
        ax.set(xlabel='Giorni dalla prima perturbazione',title=title);ax.axvline(6,color='grey',ls=':',lw=.8)
    axes[0,0].legend(fontsize=7)
    savefig(fig,'06_sequence_memory')
    fig,axes=plt.subplots(1,3,figsize=(11,4),layout='constrained')
    for ax,id in zip(axes,['C_L45_N1_a0.67','C_L405_N1_a0.67','C_L405_N16_a0.67']):
        d=load(id);pn=d['phyto_N'];fraction=pn/np.maximum(pn.sum(axis=1)[:,None],1e-30)
        im=ax.pcolormesh(d['time']/365,d['sizes'],np.log10(np.maximum(fraction.T,1e-6)),shading='auto',vmin=-4,vmax=0,cmap='magma')
        ax.set(yscale='log',xlabel='Tempo [anni]',ylabel='ESD P [µm]',title=id.replace('C_','').replace('_a0.67',''))
    fig.colorbar(im,ax=axes,label='log$_{10}$ frazione di PN')
    savefig(fig,'07_size_succession')
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for ax,id in zip(axes,['C_L405_N1_a0.67','C_L405_N16_a0.67']):
        d=load(id);t=d['time'];pc=time_mean(t,d['phyto_C'],3285.);pn=time_mean(t,d['phyto_N'],3285.)
        food=time_mean(t,d['prey_flux_C'],3285.)
        ax.plot(d['sizes'],pc/pc.sum(),'o-',label='C presente',markersize=3)
        ax.plot(d['sizes'],food/food.sum(),'s-',label='C ingerito',markersize=3)
        ax.set(xscale='log',xlabel='ESD P [µm]',ylabel='Frazione del C totale',title=id.replace('C_','').replace('_a0.67',''))
        q=ax.twinx();q.plot(d['sizes'],np.divide(pc,pn,out=np.full_like(pc,np.nan),where=pn>pn.sum()*.001),':',color='#8052a0')
        q.set_ylabel('C:N per classe',color='#8052a0')
    axes[0].legend(fontsize=8)
    savefig(fig,'08_diet_filtering')
    fig,axes=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
    for id,label in [('C_L405_N16_a0.67','h=.05'),('R_L405_N16_h025','h=.025'),('Q_L405_N16_h0125','h=.0125')]:
        d=load(id);mask=d['time']>=2920;t=d['time'][mask]-2920
        pn=d['phyto_N'][mask].sum(axis=1);pc=d['phyto_C'][mask].sum(axis=1)
        axes[0].plot(t,pn,label=label,lw=1)
        axes[1].plot(pn,pc/pn,label=label,lw=.7,alpha=.8)
    axes[0].set(xlabel='Giorni negli ultimi due anni',ylabel='PN [mmol N m$^{-3}$]',title='Variabilità con forzanti costanti')
    axes[1].set(xlabel='PN [mmol N m$^{-3}$]',ylabel='C:N di comunità',title='Traiettorie e fase dipendono dal passo')
    axes[0].legend(fontsize=8)
    savefig(fig,'09_autonomous_variability')

def contrasts(rows):
    by={r['id']:r for r in rows};summary={}
    summary['minimal_light_sweeps']=[]
    for n in [.25,1,4,16]:
        rr=[by[f'M_L{l}_N{n:g}_a0.67'] for l in [15,45,135,405,1215]]
        summary['minimal_light_sweeps'].append(dict(inflow_DIN=n,ids=[r['id'] for r in rr],
            CN=[r['CN'] for r in rr],PC=[r['phyto_C'] for r in rr],ZC=[r['zoo_C'] for r in rr],
            Zproduction=[r['production_Z_C'] for r in rr],yield_C=[r['yield_Z_per_ingested_C'] for r in rr]))
    summary['community_eta_contrasts']=[]
    for l in [45,135,405]:
        for n in [1,4,16]:
            a=by[f'C_L{l}_N{n}_a0.67'];b=by[f'C_L{l}_N{n}_a0.335']
            summary['community_eta_contrasts'].append(dict(PAR=l,inflow_DIN=n,reference=a['id'],lower_eta=b['id'],
                delta_PC_percent=100*(b['phyto_C']/a['phyto_C']-1),delta_PN_percent=100*(b['phyto_N']/a['phyto_N']-1),
                delta_Zprod_percent=100*(b['production_Z_C']/a['production_Z_C']-1) if a['production_Z_C']>1e-14 else None,
                delta_CN=b['CN']-a['CN'],delta_geomean_percent=100*(b['size_N']['geomean']/a['size_N']['geomean']-1),
                delta_evenness=b['size_N']['evenness']-a['size_N']['evenness']))
    summary['resolution_checks']=[]
    for l,n in [(45,1),(405,1),(405,16)]:
        a=by[f'C_L{l}_N{n}_a0.67']
        for suffix in ['h025','bins70','bins140']:
            b=by[f'R_L{l}_N{n}_{suffix}']
            summary['resolution_checks'].append(dict(reference=a['id'],control=b['id'],
                delta_CN_percent=100*(b['CN']/a['CN']-1),delta_Zprod_percent=100*(b['production_Z_C']/a['production_Z_C']-1),
                delta_geomean_percent=100*(b['size_N']['geomean']/a['size_N']['geomean']-1),
                delta_evenness=b['size_N']['evenness']-a['size_N']['evenness'],
                CN=b['CN'],geomean=b['size_N']['geomean'],evenness=b['size_N']['evenness']))
    # Exact symmetric N-weight decomposition of the ratio of mean C and N stocks.
    summary['community_light_decomposition']=[]
    for n in [1,4,16]:
        a=by[f'C_L45_N{n}_a0.67'];b=by[f'C_L405_N{n}_a0.67']
        na=np.array(a['spectra']['PN']);nb=np.array(b['spectra']['PN']);ca=np.array(a['spectra']['PC']);cb=np.array(b['spectra']['PC'])
        qa=np.divide(ca,na,out=np.zeros_like(ca),where=na>1e-300);qb=np.divide(cb,nb,out=np.zeros_like(cb),where=nb>1e-300)
        wa=na/na.sum();wb=nb/nb.sum();within=float(((wa+wb)/2*(qb-qa)).sum());composition=float(((qa+qb)/2*(wb-wa)).sum())
        summary['community_light_decomposition'].append(dict(inflow_DIN=n,delta_CN=b['CN']-a['CN'],within_class=within,composition=composition,
            identity_error=within+composition-(b['CN']-a['CN'])))
    da=load('T_N_then_light');db=load('T_light_then_N');summary['sequence_memory']=[]
    for t in [736,739,746,766]:
        ia=np.flatnonzero(abs(da['time']-t)<1e-7)[-1];ib=np.flatnonzero(abs(db['time']-t)<1e-7)[-1]
        pa=da['phyto_N'][ia];pb=db['phyto_N'][ib];ca=da['phyto_C'][ia];cb=db['phyto_C'][ib]
        za=float(da['production_Z_C'][ia].sum());zb=float(db['production_Z_C'][ib].sum())
        summary['sequence_memory'].append(dict(days_after_final_environment=t-736,
            CN_N_then_light=float(ca.sum()/pa.sum()),CN_light_then_N=float(cb.sum()/pb.sum()),
            PN_N_then_light=float(pa.sum()),PN_light_then_N=float(pb.sum()),
            Zprod_N_then_light=za,Zprod_light_then_N=zb,
            delta_Zprod_percent=100*(za/zb-1),
            composition_half_L1=float(.5*abs(pa/pa.sum()-pb/pb.sum()).sum())))
    summary['paired_eta_robustness']=[]
    for l,n in [(45,1),(405,1),(405,16)]:
        for suffix,reference in [('h0125',f'Q_L{l}_N{n}_h0125'),('bins70',f'R_L{l}_N{n}_bins70')]:
            a=by[reference];b=by[f'B_L{l}_N{n}_a0335_{suffix}']
            summary['paired_eta_robustness'].append(dict(reference=reference,lower_eta=b['id'],
                delta_Zprod_percent=100*(b['production_Z_C']/a['production_Z_C']-1),
                delta_geomean_percent=100*(b['size_N']['geomean']/a['size_N']['geomean']-1),
                delta_evenness=b['size_N']['evenness']-a['size_N']['evenness'],delta_CN=b['CN']-a['CN']))
    summary['fine_step_sequence_memory']=[]
    da=load('B_T_N_then_light');db=load('B_T_light_then_N')
    for t in [736,739,746,766]:
        ia=np.flatnonzero(abs(da['time']-t)<1e-7)[-1];ib=np.flatnonzero(abs(db['time']-t)<1e-7)[-1]
        a=da['production_Z_C'][ia].sum();b=db['production_Z_C'][ib].sum()
        summary['fine_step_sequence_memory'].append(dict(days_after_final_environment=t-736,delta_Zprod_percent=float(100*(a/b-1))))
    return summary

def main():
    rows=[extract(c) for c in registry()+control_registry()+robustness_registry()]
    data=dict(cases=rows,contrasts=contrasts(rows),
        audit=dict(n_cases=len(rows),max_relative_residual_C=max(r['residual_CN'][0] for r in rows),
            max_relative_residual_N=max(r['residual_CN'][1] for r in rows),minimum_live_state=min(r['minimum_live_state'] for r in rows),
            minimum_DIC=min(r['minimum_DIC'] for r in rows)),
        analysis_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (DOC/'summary.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    keys=[k for k in rows[0] if k not in ['spectra','size_C','size_N','residual_CN','final_window_days','long_window']]
    with (DOC/'summary.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=keys+['ESDgeo_N','ESD95_N','evenness_N']);writer.writeheader()
        for r in rows:writer.writerow({**{k:r[k] for k in keys},'ESDgeo_N':r['size_N']['geomean'],'ESD95_N':r['size_N']['ESD95'],'evenness_N':r['size_N']['evenness']})
    figures(rows)
    print(json.dumps(data['audit'],indent=2));print(json.dumps(data['contrasts'],indent=2))

if __name__=='__main__':main()
