"""Scientific figures for the prespecified POM challenge; no synthetic data."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scripts.run_pom_validation import DOC,OUT,identifier,build

COLORS=['#147d92','#d78330','#8058a6','#50834b','#ba5863']
LABELS={'full':'Riferimento','own':'Chiusura per classe','flat':'Riserva C piatta','uniform':'Incontro uniforme','prior':'Fisiologia a priori','noZ':'Senza consumatori'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,
    'axes.prop_cycle':plt.cycler(color=COLORS),'figure.dpi':130,'savefig.dpi':200,'axes.titleweight':'bold'})


def save(fig,name):
    fig.savefig(OUT/(name+'.png'),facecolor='white',bbox_inches='tight')
    fig.savefig(OUT/(name+'.svg'),facecolor='white',bbox_inches='tight');plt.close(fig)


def main():
    s=json.loads((DOC/'summary.json').read_text());rows=s['rows'];by={r['id']:r for r in rows}
    get=lambda structure,l,n:by[identifier(structure,l,n)]
    structures=list(LABELS)+(['full_fine'] if 'fine_full_pattern_decisions' in s else [])
    decision=lambda st:s['fine_full_pattern_decisions']['decision'] if st=='full_fine' else s['pattern_decisions'][st]
    fig,ax=plt.subplots(figsize=(9,5),layout='constrained')
    a=np.array([[decision(st)[k]['pass_pattern'] for k in ['P1','P2','P3','P4']] for st in structures],int)
    ax.imshow(a,aspect='auto',cmap=ListedColormap(['#f2ded3','#cfe3db']),vmin=0,vmax=1)
    ax.set(xticks=range(4),xticklabels=['P1\nN e taglia','P2\nGruppi di taglia','P3\nQuantità e qualità','P4\nDiversità e Z:P'],
           yticks=range(len(structures)),yticklabels=[LABELS.get(st,'Riferimento / passo fine') for st in structures])
    for i,st in enumerate(structures):
        for j,k in enumerate(['P1','P2','P3','P4']):
            annual=sum(b[k]['pass_pattern'] for b in s['fine_full_pattern_decisions']['annual']) if st=='full_fine' else sum(b['decisions'][st][k]['pass_pattern'] for b in s['annual_pattern_decisions'])
            text=('Sì' if a[i,j] else 'No')+f'\n{annual}/5 blocchi'
            if st=='noZ' and k in ['P3','P4']:text='Non applicabile'
            ax.text(j,i,text,ha='center',va='center',fontsize=10,color='#183f43')
    ax.tick_params(length=0);ax.set_title('Compatibilità congiunta: criteri fissati prima dei nuovi run',pad=16)
    save(fig,'01_POM_matrix')

    fig,axes=plt.subplots(2,3,figsize=(11,6.7),layout='constrained')
    for i,light in enumerate([30,300]):
        for j,st in enumerate(['full','own','flat']):
            ax=axes[i,j]
            for n in [.5,2,8]:
                r=get(st,light,n);a=r['P_N'];x=a['bin_centres'];q=np.array(a['bin_fraction'])
                ax.plot(x,q,'o-',ms=3,label=f'N entrante {n:g}')
            ax.axvline(5.85,color='#999999',lw=.8,ls=':')
            ax.set(xscale='log',xlabel='ESD P [µm]',ylabel='Frazione di PN / intervallo log',
                   title=f'{LABELS[st]}, PAR {light}')
            ax.set_xticks([1,5,20,100]);ax.set_xticklabels([1,5,20,100])
    axes[0,0].legend(fontsize=8);save(fig,'02_resource_spectra')

    fig,axes=plt.subplots(2,3,figsize=(11,6.6),layout='constrained')
    quantities=[('PC','Stock P [mmol C m$^{-3}$]'),('CN','C:N della comunità'),('diet_CN','C:N della dieta ingerita'),
                ('production_Z_C','Allocazione a Z [mmol C m$^{-3}$ d$^{-1}$]'),('TTE_net_secondary','Produzione secondaria netta / NPP'),('PN_CV','Variabilità PN: CV')]
    for ax,(key,title) in zip(axes.flat,quantities):
        for n,color in zip([.5,2,8],COLORS):
            for st,style in [('full','-'),('own','--')]:
                rr=[get(st,l,n) for l in [3,30,300]]
                ax.plot([3,30,300],[r[key] for r in rr],style,marker='o',color=color,label=f'N {n:g}'+(' / full' if st=='full' else ' / own'))
        ax.set(xscale='log',xticks=[3,30,300],xticklabels=[3,30,300],xlabel='PAR [µmol fotoni m$^{-2}$ s$^{-1}$]',title=title)
    axes[0,0].legend(fontsize=7,ncols=2);save(fig,'03_quantity_quality')

    fig,axes=plt.subplots(2,2,figsize=(10,6.8),layout='constrained')
    for i,currency in enumerate(['N','C']):
        ykey='ZN_PN' if currency=='N' else 'ZC_PC'
        for st,color in zip(['full','own'],COLORS):
            rr=[get(st,l,n) for l in [3,30,300] for n in [.5,2,8]]
            x=np.array([r['P_'+currency]['evenness'] for r in rr]);z=np.array([r['Z']['evenness'] for r in rr]);y=np.log10([r[ykey] for r in rr])
            for j,v in enumerate([x,z]):
                ax=axes[i,j];ax.scatter(v,y,s=36,color=color,label=LABELS[st])
                fit=np.polyfit(v,y,1);grid=np.linspace(v.min(),v.max(),30);ax.plot(grid,np.polyval(fit,grid),color=color,lw=1,alpha=.7)
                rho=s['pattern_decisions'][st]['P4']['common_9'][currency]['rho_P' if j==0 else 'rho_Z']
                ax.text(.03,.96-(0 if st=='full' else .09),f'{st}: ρ={rho:+.2f}',transform=ax.transAxes,va='top',color=color,fontsize=9)
        for j in [0,1]:axes[i,j].set(xlabel=('Equità delle prede' if j==0 else 'Equità dei predatori')+f' (valuta {currency})',
                                      ylabel=f'log$_{{10}}$ Z{currency}:P{currency}',title='Intero piano di 9 ambienti')
    axes[0,0].legend(fontsize=8,loc='lower left');save(fig,'04_diversity_transfer')

    fig,axes=plt.subplots(2,3,figsize=(11,6.5),layout='constrained')
    for i,n in enumerate([.5,8]):
        rr=sorted([r for r in s['N_retention_recycling'] if r['N']==n],key=lambda r:r['N_recycling_rate']);x=[r['N_recycling_rate'] for r in rr]
        axes[i,0].plot(x,[100*r['direct_same_state_change'] for r in rr],'o--',label='Effetto sullo stesso stato')
        axes[i,0].plot(x,[100*r['total_trajectory_changes']['production_Z_C'] for r in rr],'s-',label='Sistema riorganizzato')
        axes[i,1].plot(x,[100*r['size_change'] for r in rr],'o-',label='Taglia geometrica')
        axes[i,1].plot(x,[100*r['evenness_change'] for r in rr],'s-',label='Equità di PN')
        for eta,color,style in [('baseline_N_flows',COLORS[0],'-'),('lower_eta_N_flows',COLORS[1],'--')]:
            for edge,marker in [('Z_excretion','o'),('remineralization','s')]:
                v=[r[eta]['flows'][r[eta]['edges'].index(edge)] for r in rr]
                axes[i,2].plot(x,v,style,marker=marker,color=color,label=edge+(' / .67' if eta=='baseline_N_flows' else ' / .335'))
        for j in [0,1]:axes[i,j].axhline(0,color='grey',lw=.7)
        for ax in axes[i]:ax.set(xscale='log',xticks=x,xticklabels=x,xlabel='Lisi e remineralizzazione N [d$^{-1}$]',title=f'DIN entrante {n:g}')
    axes[0,0].set_ylabel('Variazione allocazione Z [%]');axes[1,0].set_ylabel('Variazione allocazione Z [%]')
    axes[0,1].set_ylabel('Variazione struttura [%]');axes[1,1].set_ylabel('Variazione struttura [%]')
    axes[0,2].set_ylabel('Ritorno a DIN [mmol N m$^{-3}$ d$^{-1}$]');axes[1,2].set_ylabel('Ritorno a DIN [mmol N m$^{-3}$ d$^{-1}$]')
    for ax in axes[0]:ax.legend(fontsize=7)
    save(fig,'05_retention_recycling')

    fig,axes=plt.subplots(2,2,figsize=(10,6.4),layout='constrained')
    for ax,(st,n) in zip(axes.flat,[('full',.5),('full',8),('flat',8),('own',8)]):
        r=get(st,300,n);sp=r['spectra'];pc=np.array(sp['PC']);pn=np.array(sp['PN'])
        x=np.array(sp['sizes']);prey=np.array(sp['prey_flux_N'])
        ax.plot(x,pc/pc.sum(),label='Quota nello stock C')
        ax.plot(x,pn/pn.sum(),label='Quota nello stock N')
        ax.plot(x,prey/prey.sum(),label='Quota nel flusso di N ingerito')
        ax.set(xscale='log',xlabel='ESD P [µm]',ylabel='Frazione per classe',title=f'{LABELS[st]}, N {n:g}')
        ax.text(.98,.95,f'C:N stock {r["CN"]:.1f}\nC:N dieta {r["diet_CN"]:.1f}',transform=ax.transAxes,ha='right',va='top',fontsize=9)
    axes[0,0].legend(fontsize=7);save(fig,'06_diet_filter')

    fig,axes=plt.subplots(2,3,figsize=(11,6.8),layout='constrained')
    for i,st in enumerate(['full','own']):
        ids=[identifier(st,300,8),identifier(st,300,8)+'_h0125',identifier(st,300,8)+'_bins70']
        labels=['35 classi, h .025','35 classi, h .0125','70 classi, h .025']
        for id,label,color in zip(ids,labels,COLORS):
            r=by[id]
            with np.load(OUT/(id+'.npz')) as z:
                case=json.loads((OUT/(id+'.json')).read_text())['case'];m=build(case);y=z['state'];t=z['time']
                pn=m.traits[:,4][None,:]*y[:,m.a,0]+y[:,m.nr,0]
                k=t>=3285;axes[i,0].plot(t[k]-3285,pn[k].sum(axis=1),color=color,lw=.9,label=label)
            q=r['P_N'];axes[i,1].plot(q['bin_centres'],q['bin_fraction'],color=color,label=label)
            a=r['annual_blocks'];j=ids.index(id)
            axes[i,2].plot([j]*5,[b['CN'] for b in a],'o',color=color,alpha=.55,ms=4)
            axes[i,2].plot(j,r['CN'],'D',color=color,ms=7)
        axes[i,0].set(xlabel='Giorni nell’ultimo anno',ylabel='PN [mmol N m$^{-3}$]',title=LABELS[st]+': dinamica')
        axes[i,1].set(xscale='log',xlabel='ESD P [µm]',ylabel='Frazione di PN / intervallo',title='Spettro medio su 5 anni')
        axes[i,2].set(xticks=range(3),xticklabels=['h .025','h .0125','70 classi'],ylabel='C:N della comunità',title='Blocchi annuali e media lunga')
    axes[0,0].legend(fontsize=7);save(fig,'07_regime_robustness')

    fig,axes=plt.subplots(2,2,figsize=(9.7,6.3),layout='constrained')
    for ax,(n,phase) in zip(axes.flat,[(.5,'min'),(.5,'max'),(2,'min'),(2,'max')]):
        rr=sorted([r for r in s['memory']['contrasts'] if r['N']==n and r['phase']==phase],key=lambda r:-r['max_step'])
        for key,label in [('integrated_Zproduction_0_6','Produzione integrata giorni 0-6'),('integrated_Zproduction_6_36','Produzione integrata giorni 6-36')]:
            ax.plot([r['max_step'] for r in rr],[100*r['metrics'][key]['normalized_order_interaction'] for r in rr],'o-',label=label)
        ax.axhline(0,color='grey',lw=.7);ax.set(xscale='log',xticks=[.025,.0125,.00625],xticklabels=['.025','.0125','.00625'],
            xlabel='Passo massimo [d]',ylabel='Interazione d’ordine / controllo [%]',title=f'N {n:g}, fase {phase}')
        ax.invert_xaxis()
    axes[0,0].legend(fontsize=8);save(fig,'08_sequence_memory')

    fig,axes=plt.subplots(2,2,figsize=(9.7,6.3),layout='constrained')
    gardens=s['reserve_common_garden']
    for ax,(n,phase) in zip(axes.flat,[(.5,'min'),(.5,'max'),(2,'min'),(2,'max')]):
        for key,label in [('initial_uptake','Uptake N iniziale'),('integrated_uptake','Uptake N integrato 1d'),('integrated_Zproduction','Allocazione Z integrata 1d')]:
            values=[]
            for h in [.025,.0125,.00625]:
                a=next(r for r in gardens if r['N']==n and r['phase']==phase and r['max_step']==h and r['arm']=='NeLl')
                b=next(r for r in gardens if r['N']==n and r['phase']==phase and r['max_step']==h and r['arm']=='LeNl')
                values.append(100*(a[key]/b[key]-1))
            ax.plot([.025,.0125,.00625],values,'o-',label=label)
        ax.axhline(0,color='grey',lw=.7);ax.set(xscale='log',xticks=[.025,.0125,.00625],xticklabels=['.025','.0125','.00625'],
            xlabel='Passo massimo [d]',ylabel='(NeLl / LeNl - 1) [%]',title=f'N {n:g}, fase {phase}')
        ax.invert_xaxis()
    axes[0,0].legend(fontsize=8);save(fig,'09_reserve_common_garden')
    if 'fine_full_pattern_decisions' in s:
        fig,axes=plt.subplots(2,3,figsize=(11,6.4),layout='constrained')
        for j,(light,n) in enumerate([(3,.5),(300,.5),(300,8)]):
            ids=[identifier('full',light,n)+suffix for suffix in ['', '_h0125','_h00625','_h003125']]
            rr=[by[id] for id in ids]
            for i,key in enumerate(['CN','geomean_N']):
                ax=axes[i,j];values=[r['CN'] if key=='CN' else r['P_N']['geomean'] for r in rr]
                ax.plot([r['step'] for r in rr],values,'o-',color=COLORS[j])
                for r in rr:
                    ax.scatter([r['step']]*5,[a[key] for a in r['annual_blocks']],color=COLORS[j],alpha=.35,s=17)
                ax.set(xscale='log',xticks=[.025,.0125,.00625,.003125],xticklabels=['.025','.0125','.00625','.003125'],
                    xlabel='Passo massimo [d]',ylabel='C:N' if key=='CN' else 'Taglia geometrica PN [µm]',title=f'PAR {light}, N {n:g}')
                ax.invert_xaxis()
        save(fig,'10_fine_whole_plane')

    if s.get('rare_reintroduction'):
        fig,axes=plt.subplots(2,2,figsize=(10.4,6.5),layout='constrained')
        for i,parent in enumerate(['V_full_L300_N0.5','V_full_L300_N0.5_h0125']):
            base='B_'+parent[2:]
            for suffix,label,color,style in [('_resident','Residente',COLORS[0],'-'),('_seeded','Solo P: 10$^{-5}$',COLORS[1],'-'),
                ('_Zseeded','Solo Z: 10$^{-5}$',COLORS[2],'--'),('_PZseeded','P e Z: 10$^{-5}$',COLORS[3],'-'),
                ('_PZseeded_small','P e Z: 10$^{-7}$',COLORS[4],':')]:
                id=base+suffix;meta=json.loads((OUT/(id+'.json')).read_text());m=build(meta['case'])
                with np.load(OUT/(id+'.npz')) as z:
                    t=z['time'];y=z['state'];pn=m.traits[:,4][None,:]*y[:,m.a,0]+y[:,m.nr,0]
                    pc=m.traits[:,3][None,:]*y[:,m.a,0]+y[:,m.cr,0]+4*y[:,m.nr,0]
                    axes[i,0].plot(t,pc.sum(axis=1)/pn.sum(axis=1),style,color=color,lw=1.2,label=label)
                row=next(r for r in s['rare_reintroduction'] if r['id']==id);a=row['P_N']
                axes[i,1].plot(a['bin_centres'],a['bin_fraction'],style,color=color,label=label)
            axes[i,0].set(xlabel='Giorni dalla redistribuzione',ylabel='C:N della comunità',title='Stato parentale h '+('.025' if i==0 else '.0125'))
            axes[i,1].set(xscale='log',xlabel='ESD P [µm]',ylabel='Frazione PN / intervallo',title='Spettro medio giorni 365-730')
        axes[0,0].legend(fontsize=7,ncols=2);save(fig,'11_rare_trait_pool')

    if (DOC/'invasion_diagnostic.json').exists():
        inv=json.loads((DOC/'invasion_diagnostic.json').read_text())
        fig,axes=plt.subplots(3,2,figsize=(10.2,8.2),layout='constrained')
        for i,row in enumerate(inv['rows']):
            for j,key,limit in [(0,'phyto',10.),(1,'zoo',50.)]:
                d=row[key];x=np.array(d['sizes']);v=np.array(d['invasion_rate_per_day']);rare=np.array(d['rare_in_parent']);sel=x<=limit
                ax=axes[i,j]
                for israre,color,label in [(False,COLORS[0],'Residente'),(True,COLORS[1],'Rara nel residente')]:
                    kk=sel&(rare==israre)
                    ax.plot(x[kk],v[kk],'o',color=color,label=label,ms=4)
                ax.axhline(0,color='grey',lw=.8);ax.set(xscale='log',xlabel='ESD [µm]',ylabel='Crescita specifica [d$^{-1}$]',
                    title=('Produttori' if j==0 else 'Consumatori')+' / '+row['id'].replace('V_full_L300_N0.5','h .025').replace('h .025_h0125','h .0125').replace('h .025_h00625','h .00625'))
        axes[0,0].legend(fontsize=8);save(fig,'12_local_invasion')
    print('Created scientific figures (PNG + editable SVG).')

if __name__=='__main__':main()
