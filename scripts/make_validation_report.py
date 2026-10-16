"""Render a compact Italian scientific development report from actual results."""
from pathlib import Path
import json,csv,subprocess,html,sys,re
import numpy as np
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,Image,KeepTogether
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
import matplotlib
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent
FIG=ROOT/'results/figures';OUTPUT=ROOT/'output/pdf/SPLAS_Sviluppo_Verifica_2026-10-16.pdf'
S=json.load(open(ROOT/'research/validation_summary.json'));CULT=json.load(open(ROOT/'results/cultures/validation.json'))
ASSESS=json.load(open(ROOT/'research/acceptance_assessment.json')) if (ROOT/'research/acceptance_assessment.json').exists() else {}
DATA=json.load(open(ROOT/'data/processed/cultures.json'));MECH=json.load(open(ROOT/'results/mechanisms/mechanisms.json'))
REG=json.load(open(ROOT/'research/experiment_registry.json'))
for extra in ['supplement_registry.json','carbon_boundary_registry.json','recommended_registry.json','recommended_refinement_registry.json','export_intervention_registry.json']:
 if (ROOT/'research'/extra).exists():REG+=json.load(open(ROOT/'research'/extra))
CASES={}
for c in REG:
 p=ROOT/'results/thesis_replay'/(c['id']+'.json')
 if p.exists():
  r=json.load(open(p))
  if r['id'] not in S['missing_cases']:CASES[r['id']]=r
if S['missing_cases'] and '--allow-partial' not in sys.argv:raise RuntimeError('Do not present an incomplete run as final')
fontdir=Path(matplotlib.get_data_path())/'fonts/ttf'
for name,file in [('DV','DejaVuSans.ttf'),('DVB','DejaVuSans-Bold.ttf'),('DVI','DejaVuSans-Oblique.ttf')]:pdfmetrics.registerFont(TTFont(name,str(fontdir/file)))
pdfmetrics.registerFontFamily('DV',normal='DV',bold='DVB',italic='DVI',boldItalic='DVB')
styles=getSampleStyleSheet();styles.add(ParagraphStyle(name='BodyDV',fontName='DV',fontSize=9.2,leading=13.2,spaceAfter=7,textColor=colors.HexColor('#243540')))
styles.add(ParagraphStyle(name='TitleDV',fontName='DVB',fontSize=31,leading=35,spaceAfter=16,textColor=colors.HexColor('#14526a')))
styles.add(ParagraphStyle(name='HDV',fontName='DVB',fontSize=16,leading=21,spaceAfter=11,textColor=colors.HexColor('#14526a')))
styles.add(ParagraphStyle(name='SmallDV',fontName='DV',fontSize=7.7,leading=10.5,spaceAfter=6,textColor=colors.HexColor('#485b66')))
styles.add(ParagraphStyle(name='TableDV',fontName='DV',fontSize=7.6,leading=10))
styles.add(ParagraphStyle(name='FormulaDV',fontName='DV',fontSize=10.5,leading=16,spaceBefore=5,spaceAfter=9,leftIndent=10,textColor=colors.HexColor('#14526a')))
story=[];W=174*2.83465

def tidy(s):
 s=str(s)
 if any(x in s for x in ['http','python','.py','SHA256']):return s
 s=re.sub(r'(?<=[a-zàèéìòù])(?=\d)',' ',s)
 s=re.sub(r'(?<=\d)(?=(?:mmol|mol|mg|pg|µmol|µm|m²|m³|m⁻|°C|giorni|classi|celle|test|repliche|m\b|d\b))',' ',s)
 s=re.sub(r'(?<=[,;])(?=\S)',' ',s)
 return s.replace('m⁻²s⁻¹','m⁻² s⁻¹')
def linked(s):
 def wrap(match):
  token=match.group(0);value=token.rstrip(').,;]');tail=token[len(value):]
  target=value if value.startswith('http') else 'https://doi.org/'+value
  return '<link href="'+html.escape(target,quote=True)+'" color="#166b8c">'+html.escape(value)+'</link>'+tail
 return re.sub(r'https?://[^\s<>]+|(?<![\w/])10\.\d{4,9}/[A-Za-z0-9().;:/_-]+',wrap,tidy(s))
def p(s,small=False):story.append(Paragraph(linked(s),styles['SmallDV' if small else 'BodyDV']))
def h(s):story.append(Paragraph(s,styles['HDV']))
def page():story.append(PageBreak())
def formula(s):story.append(Paragraph(s,styles['FormulaDV']))
def tbl(rows,widths=None):
 data=[[Paragraph(linked(v),styles['TableDV']) for v in row] for row in rows]
 t=Table(data,colWidths=widths or [W/len(rows[0])]*len(rows[0]),repeatRows=1,hAlign='LEFT')
 t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#dfeef2')),('FONTNAME',(0,0),(-1,0),'DVB'),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),('LINEBELOW',(0,0),(-1,0),.7,colors.HexColor('#9abac7')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f4f7f8')])]))
 story.append(t);story.append(Spacer(1,9))
def fig(name,height=None):
 path=FIG/(name+'.png')
 if not path.exists():return
 from PIL import Image as PILImage
 wi,he=PILImage.open(path).size;hh=W*he/wi
 if height and hh>height:story.append(Image(str(path),width=height*wi/he,height=height))
 else:story.append(Image(str(path),width=W,height=hh))
 story.append(Spacer(1,7))
def n(x,d=3):return 'n.d.' if x is None else f'{x:.{d}g}'
def pct(x):return 'n.d.' if x is None else f'{100*x:+.2f}%'
def summ(key):return CASES[key]['summary']
def balance(key):return next((r for r in S.get('last_window_element_balance',[]) if r['id']==key),{})
def footer(canvas,doc):
 canvas.setStrokeColor(colors.HexColor('#b7cbd3'));canvas.line(18*2.83465,17*2.83465,192*2.83465,17*2.83465)
 canvas.setFont('DV',7);canvas.setFillColor(colors.HexColor('#59717c'));canvas.drawString(18*2.83465,12*2.83465,'SPLAS-CN 0.2 | verifica del 16 ottobre 2026');canvas.drawRightString(192*2.83465,12*2.83465,str(doc.page))

#1
story.append(Spacer(1,18));story.append(Paragraph('SPLAS-CN<br/>Sviluppo e verifica',styles['TitleDV']))
p('Dalle tesi De Pasquale → Fianchini a una formulazione conservativa con stechiometria cellulare dinamica.')
p('<b>Esito:</b> implementazione funzionante e verificata nei bilanci C/N. La valutazione fisiologica è parziale; la similarità dei pattern non equivale a una validazione biogeochimica di campo. Fallimenti e sensibilità numeriche sono parte del risultato.')
rows=[['Evidenza completata','Risultato e portata'],['Verifica del codice','20 test: bilanci indipendenti, positività, trasporto analitico, R*, DOP853, confronto esplicito/implicito.'],['Serie delle tesi',f"{S.get('decadal_case_count',S['current_cases'])} serie decennali e {S.get('annual_continuation_count',0)} confronti annuali; configurazioni e impronte del codice registrate."],['Fisiologia','4 colture; parametri stimati su Ot/Tw, Tp/Micromonas esclusi dall’ottimizzatore. Valutazione esplorativa condizionata ai dati iniziali.'],['Controlli di meccanismo','8 assemblaggi Monod,25 stati luce×DIN,16 fasi di coltura,18 esperimenti di mixing; recupero luce/N e chemostato.'],['Chiusura dei bilanci',f"Errore relativo massimo: C {S['mass_C_max_relative']:.2e}; N {S['mass_N_max_relative']:.2e}. Minimo stato: {S['minimum_state']:.2g}."],['Dominio non validato','P,Fe,Si,O₂, microbiologia esplicita, termodinamica completa, predazione e flussi oceanici osservati.']]
tbl(rows,[W*.32,W*.68])
p('Il codice storico rimane disponibile accanto al nuovo nucleo. I confronti sono reinterpretazioni degli esperimenti con flussi corretti, non riproduzioni del vecchio eseguibile.',small=True)
if S['missing_cases']:p('<b>BOZZA INCOMPLETA:</b> mancano '+str(len(S['missing_cases']))+' casi.',small=True)
#2
page();h('1. Che cosa è stato sviluppato')
p('Il modello distingue abbondanza cellulare, struttura, riserva di solo carbonio e pool organico di azoto che contiene anche carbonio. La stechiometria emerge dall’evoluzione di questi inventari; non viene ricostruita a posteriori da una conversione fissa.')
formula('B<sub>C,i</sub> = c<sub>0,i</sub>a<sub>i</sub> + R<sub>C,i</sub> + ρ<sub>R</sub>R<sub>N,i</sub><br/>B<sub>N,i</sub> = n<sub>0,i</sub>a<sub>i</sub> + R<sub>N,i</sub>')
p('ρ<sub>R</sub>=4 molC/molN è un’ipotesi ridotta per il pool organico labile. Con c₀/n₀≥4 e stati non negativi, il C:N cellulare non può scendere sotto4. Il pool non rappresenta esplicitamente proteine,RNA e nitrato vacuolare separati.')
tbl([['Processo','Nuova rappresentazione'],['Fissazione e divisione','Fotosintesi dipendente da PAR,DIC,N e capacità di accumulo; divisione co-limitata dalle riserve.'],['Assimilazione N','Trasferimento del carbonio associato all’azoto organico e costo respiratorio esplicito; uptake nullo con riserva C nulla.'],['Mantenimento e carestia','Costo basale; deficit pagato mediante attrizione cellulare con C respirato e N restituito al DIN.'],['Grazing','Disponibilità φᵀB_C; massimo d’ingestione rispettato; consumo con C:N effettivo della preda.'],['Consumatori','C:N omeostatico5.5; crescita vincolata separatamente da C e N; eccessi respirati/escreti.'],['Detrito e DOM','C e N separati, classi di particelle, lisi e remineralizzazione distinte.'],['Colonna1D','Volumi finiti, flussi uguali e opposti alle facce, sinking ed export espliciti.']],[W*.28,W*.72])
p('<b>Limite energetico:</b> il costo di assimilazione N è un bilancio di carbonio, non una chiusura di ATP,elettroni e redox. Il DIN unifica nitrato e ammonio; il costo1.5 molC/molN resta un prior. La correttezza atomica non dimostra da sola una cinetica biologica corretta.',small=True)
#3
page();h('2. Genealogia, unità e correzioni')
p('Moloney (1988) assume rapporti C:N di massa fissi per gruppi. L’allegato MoloneyField0D v0.1 è un benchmark di sviluppo conservativo, non un rilascio originale verificato del modello1991. SPLAS-CN è un’ulteriore ipotesi fisiologica, con provenienza dichiarata.')
tbl([['Punto','Conseguenza verificabile'],['Carbonio da volume','Montagnes etal1994:0.109 V^0.991 pgC/cellula. Nel codice storico compare0.0109: differenza×10.'],['Conversioni elementari','1 mmolC=12mgC;1 mmolN=14mgN. C:N di massa ×14/12 → molare. Quota per10⁸ cellule: pgC/120,pgN/140.'],['Bilancio storico del N','La crescita di P_N e l’uptake espresso tramite P/Q non chiudono automaticamente lo stesso inventario. Il nuovo nucleo usa un unico flusso per ogni trasferimento.'],['Contrazione del kernel','Matrice φ(preda,predatore): somma sull’asse delle prede. Un test confronta l’ingestione con un oracolo matriciale indipendente.'],['Diffusione','10⁻⁴m²/s=8.64m²/d. Per D variabile si usa la divergenza del flusso alle facce.'],['Bordi e sorgente N','40 celle reali da5m per200m. Sorgente uniforme0.1S, distinta dal flusso al fondo; non proporzionale all’indice dello strato.'],['Periodicità','Luce storica sin²(2πt/365):182.5d. Nutriente sin²(πt/365):365d. Si distinguono formule e interpretazioni verbali.']],[W*.29,W*.71])
p('Il coefficiente di Montagnes è una relazione empirica, non una legge universale. Per due specie con volume medio1.8µm³ il riferimento misurato C differisce di oltre un fattore2. Nelle colture uso le quote replete osservate; quindi la prova non è una previsione basata sulla sola taglia.',small=True)
#4
page();h('3. Verifica numerica e integratori')
p('Il riferimento esplicito usa SSPRK3 con stima locale embedded RK2, rifiuto dei passi e nessun clipping. Il calcolo lungo usa SSPRK3 positivo per i processi non stiff e un solve implicito dell’uptake N. La combinazione IMEX è di <b>primo ordine</b>; il controllo locale del riferimento esplicito non certifica l’errore di splitting.')
p('Il solve implicito conserva i medesimi flussi C e N e risolve un’equazione monotona per il DIN comune, includendo feedback della quota N e disponibilità energetica. Il limite CFL del trasporto e i controlli degli stadi proteggono la positività. La riduzione del passo resta necessaria.')
rr=MECH['fixed_CN_assembly'][0]
tbl([['Oracolo','Esito'],['R* Monod, nessun predatore',f"DIN numerico {rr['final_DIN']:.12g}; soluzione algebrica {rr['analytic_Rstar']:.12g} mmolN/m³; errore {rr['Rstar_absolute_error']:.2e}. Vince1µm."],['Trasporto passivo','Modo coseno con decadimento esponenziale discreto, pareti no-flux; controllo anche con D variabile.'],['Riferimento indipendente','Confronto con SciPy DOP853 a tolleranza stretta; convergenza esplicita e IMEX verso il riferimento.'],['Conservazione locale','200 stati casuali con bilanci ricostruiti da struttura/riserve/consumatori/detrito e ingressi/uscite.'],['Sfide fisiologiche','DIN0;DIC0;buio e riserva C0; riserve presenti; ingestione a biomassa elevata; tracce prossime all’estinzione.']],[W*.32,W*.68])
p('Durante lo sviluppo un integratore positivo ma senza controllo adeguato della stiffness produceva un R* oscillante spurio. Il bilancio di massa era comunque chiuso. Questo fallimento è conservato nella storia di sviluppo e motiva gli oracoli indipendenti.',small=True)
#5
page();h('4. Dati, stima e valutazione fisiologica')
p('Dati Liefer:81 repliche originali;6 escluse per contaminazione batterica nota. I16 valori medi C,N repleti/starvati concordano con la tabella pubblicata entro0.4%. L’archivio PANGAEA pubblicato nel2026 segnala l’assenza della consueta approvazione finale degli autori: il controllo della tabella non elimina questa riserva.')
p('Colture a85µmol fotoni m⁻²s⁻¹, ciclo12:12. Inizializzazione dal primo campione successivo alla diluizione con DIN misurato, senza inventare il fattore di diluizione. Ot/Tw stimano6 parametri comuni; Tp/Ms non entrano nell’ottimizzazione. Quote replete e stato iniziale sono osservati anche per le specie escluse.')
p('<b>Trasparenza:</b> i punteggi Tp/Ms erano già stati esaminati in un primo ciclo di sviluppo. La valutazione finale è esplorativa, non un holdout prospettico intatto. Omta2024 usa le stesse colture: non è un ulteriore dataset indipendente.')
rows=[['Specie','Uso','Errore log C:N<br/>dinamico','Redfield<br/>fisso','C:N iniziale<br/>costante']]
for m in CULT['metrics']:rows.append([m['species'],'stima' if m['split']=='calibration' else 'valutazione',n(m['CN_RMSLE']),n(m['CN_Redfield_RMSLE']),n(m['CN_initial_ratio_RMSLE'])])
tbl(rows,[W*.12,W*.20,W*.23,W*.20,W*.25])
formula('Errore log = √media[ln(predetto/osservato)²]; minore è meglio.')
p('Si escludono i punti iniziali già usati per inizializzare. Il C:N osservato è la media dei rapporti tra repliche. I controlli a rapporto fisso sono nulli descrittivi, non modelli dinamici con cinetiche perfettamente abbinate.')
p(f"Condizionamento locale del Jacobiano della stima: {CULT['jacobian_condition']:.0f}. Due parametri raggiungono un limite di ricerca. Non vengono attribuiti intervalli di confidenza affidabili o identificabilità globale ai valori stimati.",small=True)
#6
page();h('5. Dove la fisiologia funziona e dove fallisce')
fig('01_cultures')
p('Punti: media±SD tra3 colture; linee: previsioni ai tempi di campionamento. Gli scarti possono superare la variabilità tra repliche anche nelle specie usate per stimare i parametri. La risposta generale alla carenza N non equivale a una descrizione quantitativa soddisfacente di ogni traiettoria.',small=True)
rows=[['Specie esclusa','C:N finale predetto','Osservato','Interpretazione']]
for key in ['Tp','Ms']:
 m=next(x for x in CULT['metrics'] if x['species']==key);rows.append([key,n(m['final_CN_predicted']),n(m['final_CN_observed']),'Sovrastima; evitare trasferibilità universale.'])
tbl(rows,[W*.19,W*.22,W*.19,W*.40])
p('Su Micromonas la persistenza dei valori iniziali è più accurata per C,N e densità cellulare. La temperatura6°C è vicina all’optimum del ceppo: un Q10 universale non rappresenta l’adattamento termico. Differenze di accumulo e regolazione tra taxa restano una priorità scientifica, non un dettaglio da assorbire in una curva di taglia.')
if S.get('culture_phase_uncertainty'):
 q=S['culture_phase_uncertainty'];p('Spostando il ciclo diurno di0/6/12/18ore rispetto allo stato iniziale, il C:N finale di Tp varia tra '+n(q['Tp']['min_final_CN'])+' e '+n(q['Tp']['max_final_CN'])+'; quello di Micromonas tra '+n(q['Ms']['min_final_CN'])+' e '+n(q['Ms']['max_final_CN'])+'. Le quattro fasi non risolvono gli scarti rispetto ai dati.',small=True)
#7
page();h('6. De Pasquale: dal meccanismo alla comunità')
fig('09_predator_assembly',height=210)
p('Il controllo Monod ripercorre l’aggiunta progressiva dei predatori a partire da2µm, poi2.19µm e classi successive. La controparte a riserve usa la stessa famiglia di kernel ma differisce per fisiologia e riciclo. Gli spettri sono densità per intervallo logaritmico; non si confrontano valori per bin non normalizzati.')
tbl([['Famiglia','Confronto eseguito'],['2.6','N-P e assemblaggi progressivi; R* e taglia vincente; controllo Monod e quote.'],['3.3-3.4','Chiusure di mortalità Z; sorgenti N da0.001 a10; taglie, biomassa, evenness e flussi.'],['3.5.1-4','Selettività, mortalità/sinking P,35/70/140 classi, resa C e saturazione del grazing.'],['3.5.5','Sorgente sin² a fasi diverse; controllo costante con uguale apporto medio.'],['4','Colonna, luce/ombreggiamento, D costante/variabile, particelle ed export.']],[W*.26,W*.74])
p('Conservare le equazioni storiche per il controllo didattico non implica mantenerne gli errori di unità o interpretazione. I casi quota non sono una sostituzione uno-a-uno delle curve Monod.',small=True)
#8
page();h('7. Nutrienti, taglie e trasferimento trofico')
fig('03_bottom_up')
p('La relazione tra apporto N e taglie è non monotona. Anche nel controllo con DIC mantenuto elevato, taglia media e percentile95 diminuiscono nella parte più arricchita della serie. Il pattern storico di aumento delle taglie non è quindi recuperato in modo robusto. Selettività, qualità stechiometrica delle prede e mortalità sono ipotesi da discriminare; l’evenness per taglia non è biodiversità tassonomica.')
p('Le taglie sintetiche sono media geometrica e percentile95 pesati sullo stock N. Il percentile è un osservabile diverso dalla massima classe con biomassa positiva; evita di interpretare tracce numeriche come presenza ecologica.',small=True)
p('TTE di flusso=produzione C dello zooplancton/NPP_C nello stesso intervallo; Z:P è un rapporto tra stock. La loro relazione non viene imposta. Un’associazione tra diversità e Z:P non valida automaticamente un meccanismo di efficienza trofica.')
p(f"{len(S['drifting_cases'])} casi mostrano una differenza oltre10% tra la biomassa N media dell’ultimo intervallo annuale e quella del precedente. Dieci anni di integrazione non certificano l’equilibrio. Gli intervalli effettivi sono registrati nei JSON dei risultati.",small=True)
# Carbon boundary and a derived trophic hypothesis
page();h('8. Il confine del carbonio cambia l’esperimento')
fig('10_carbon_boundary',height=190)
p('L’azoto aggiunto non porta automaticamente un rifornimento di carbonio. Nei box con rilassamento DIC0.01/d il carbonio può diventare limitante. La serie aggiuntiva con reservoir1/d mantiene il DIC vicino al valore imposto2200; è un controllo ideale, non una parametrizzazione di scambio aria-mare.')
p('Anche la colonna aperta di riferimento riceve N ma non reintegra il C esportato. Il DIC scende localmente fino a circa1µmol/L: le differenze non possono essere attribuite alla sola limitazione N o alla luce. La configurazione riveduta aggiunge un confine C dichiarato e più classi di detrito.')
formula('Resa C del grazer ≤ min[0.33, 5.5 × 0.67 / (C:N della preda)]')
fig('11_food_quality_hypothesis',height=175)
p('La soglia C:N≈11.2 separa i due vincoli nel modello. È una previsione condizionata alle efficienze e all’omeostasi assunte; non un dato sperimentale. La taglia resta un tratto strutturale: volume e preferenza di preda non si aggiornano con l’accumulo di riserve.',small=True)
#9
page();h('9. Fianchini: aperto, chiuso e luce')
fig('04_vertical_profiles')
p('La controparte di riferimento usa200m,40 celle da5m e70+70 classi. La sorgente volumetrica è0.4mmolN/m³/d e quella al fondo20mmolN/m²/d, per un apporto totale100mmolN/m²/d. Sono forzanti idealizzate: nessun sito reale è associato alla colonna.')
p('L’irradianza959.83W/m² viene convertita assumendo che sia PAR con4.57µmol/J. La conversione è un’ipotesi spettrale; una variante attribuisce il43% a PAR se il valore indicava radiazione totale. La fisiologia non include fotoinibizione: la risposta ad alte luci non è validata sperimentalmente qui.')
p('I profili C:N sono rapporti degli inventari medi, non media non pesata dei rapporti. Un massimo profondo di biomassa, da solo, non identifica il processo responsabile e non equivale necessariamente a un massimo di clorofilla.')
rows=[['Serie,7 livelli luce','ρ C:N','ρ evenness','ρ e-ratio','R² export~luce']]
labels={'Fianchini_constant_open':'Costante aperto','Fianchini_constant_closed':'Costante chiuso','Fianchini_semiannual_200m':'Stagionale200m','Fianchini_semiannual_100m':'Stagionale100m'}
for a in S.get('pattern_assessment',[]):
 if a['family'] not in labels:continue
 m=a['metrics'];rows.append([labels[a['family']],*[n(m.get(k,{}).get('spearman_rho')) for k in ['phyto_CN_inventory_ratio','phyto_evenness_N','e_ratio_at_bottom']],n(m.get('bottom_export_C_mmol_m2_d',{}).get('linear_R2'))])
if len(rows)>1:tbl(rows,[W*.32,W*.16,W*.16,W*.16,W*.20])
p('Controllo dei pattern della tesi: C:N crescente, evenness N decrescente, export circa lineare ed e-ratio decrescente. ρ è la correlazione di rango tra sette esperimenti deterministici, non un test statistico su osservazioni. R² descrive il fit lineare interno; nel sistema chiuso export ed e-ratio sono nulli. Il CSV riporta evenness con pesi N e C.',small=True)
for a in S.get('pattern_assessment',[]):
 if a['family']=='Fianchini_constant_open':
  q=a['metrics']['phyto_CN_inventory_ratio'];e=a['metrics']['phyto_evenness_N']
  if q['spearman_rho']<0 and e['spearman_rho']>0:p('Nel sistema aperto i segni di C:N ed evenness sono opposti al pattern storico: tra0.7 e1.3× luce, C:N cambia del '+n(100*(q['values'][-1]/q['values'][0]-1))+'%. Questo mancato recupero del pattern resta un risultato; la decomposizione successiva ne descrive le componenti senza identificare da sola la causa.',small=True)
#10
page();h('10. Stagionalità e attribuzione del C:N')
fig('05_equal_light_dose')
p('Forzante semestrale storica, annuale e costante hanno la stessa dose annua di PAR. La presenza di due massimi luminosi nella formula storica impedisce di usare il semplice conteggio dei bloom come prova indipendente di una dinamica stagionale emergente.')
if S['stoichiometric_decompositions']:
 rows=[['Luce0.7→1.3','Δ C:N','Entro classe','Composizione']]
 for d in S['stoichiometric_decompositions']:rows.append(['Aperto' if 'open' in d['reference'] else 'Chiuso',n(d['total_CN_change']),n(d['within_class_term']),n(d['composition_term'])])
 tbl(rows,[W*.28,W*.24,W*.24,W*.24])
p('La scomposizione simmetrica usa pesi di biomassa N e separa variazione dei rapporti entro le classi da redistribuzione dei pesi. È un’identità contabile; non una prova di causalità unica. Per una classe assente si usa la quota dell’altro caso come convenzione dichiarata.',small=True)
#11
page();h('11. Profondità critica: un test più circoscritto')
fig('08_mixing_mechanism')
p('Diciotto esperimenti annuali incrociano mixing prescritto25/50/100m, PAR40/120/360 e presenza/assenza di grazers. D=10⁻²m²/s nello strato mescolante e10⁻⁶ sotto; luce e mixing variano indipendentemente.')
p('Il confronto storico200→100m cambia anche l’apporto areale N:0.4×H+20 passa da100 a60mmolN/m²/d, oltre a troncare l’inventario iniziale. Questo confronto di domini non isola causalmente la profondità di mixing.')
p('Il diagnostico Dcr cerca il primo attraversamento positivo→negativo dell’integrale verticale del tasso biologico di crescita netta C, a stato istantaneo congelato e senza interpretare il trasporto come crescita. Un valore assente è registrato come tale. Una soglia locale della crescita non viene chiamata profondità critica integrata.')
p('<b>Portata:</b> il diagnostico non soddisfa automaticamente tutte le ipotesi di Sverdrup: uniformità della popolazione nello strato, tempi di rimescolamento, acclimatazione e perdita. Le colture non validano l’innesco di bloom oceanici. Zeu è ottica; se raggiunge il fondo del dominio è censurata, e non viene identificata con la profondità del mixing.')
p('I lavori osservativi recenti su luce e stratificazione sostengono la necessità di inferire il mixing indipendentemente. Qui si verifica la distinzione dei meccanismi; non si riproduce una serie osservata di un bacino.',small=True)
#12
page();h('12. Robustezza numerica: che cosa regge')
fig('07_refinement',height=230)
rows=[['Confronto','Δ P_N','Δ C:N','Δ NPP','Δ export C']]
for c in S['numerical_comparisons']:
 rows.append([c['alternative'].replace('D07_size_bins_','Taglie').replace('F07_layers_','Celle').replace('F08_detritus_bins_','Detrito').replace('F09_dt_','dt1D=').replace('D12_dt_','dt0D='),*[('n.d.' if c.get(k+'_relative_change') is None else f"{100*c[k+'_relative_change']:+.1f}%") for k in ['phyto_N_mmol_m2','phyto_CN_inventory_ratio','NPP_C_mmol_m2_d','bottom_export_C_mmol_m2_d']]])
tbl(rows,[W*.30,W*.175,W*.175,W*.175,W*.175])
p('Differenze nell’ultimo intervallo annuale possono combinare errore del metodo e cambiamento di fase di dinamiche persistenti. Non sono intervalli di confidenza. Se non diminuiscono con il raffinamento, non è dimostrata la convergenza della specifica metrica; una curva simile non basta.',small=True)
#13
page();h('13. Export: numerica e ipotesi biologiche')
fig('06_structural_sensitivity')
p('Export gravitazionale al fondo e perdita verso livelli trofici superiori sono contatori distinti. e-ratio=flusso C al fondo/NPP_C integrata sulla stessa colonna e sullo stesso intervallo; non si trasferisce automaticamente il valore a100m se il fondo è200m. Export da100/200m non misura sequestro a lungo termine.')
p('Due classi detritiche chiudono la massa ma interpolano in modo grossolano le velocità: la convergenza in numero di particelle va verificata. Il numero di classi non valida la relazione velocità-taglia, la mappatura delle egesta o la remineralizzazione. Mancano aggregazione,ballast,O₂ e batteri espliciti.')
p('La maggiore complessità della fisiologia può aumentare NPP senza aumentare export. La letteratura2025-2026 distingue inoltre export superficiale e trasferimento in profondità. Le forzanti ideali non autorizzano un confronto quantitativo diretto con medie globali.',small=True)
# Revised reference and annual interventions
page();h('14. La configurazione riveduta alla prova')
p('Il riferimento R00 integra dieci anni con DIC rilassato in ogni cella a2200mmolC/m³, tempo100giorni,16 classi detritiche e passo massimo0.0125d. Usa70+70 classi di plankton,40 celle e la luce costante storica. È un esperimento ideale con ingressi/uscite C contabilizzati, non una ricostruzione di scambio aria-mare.')
if 'R00_buffered_reference' in CASES:
 r=summ('R00_buffered_reference')
 tbl([['Indicatore, ultimo intervallo annuale','Valore'],['P_N',n(r['phyto_N_mmol_m2'])+' mmolN/m²'],['C:N fitoplancton',n(r['phyto_CN_inventory_ratio'])],['NPP C',n(r['NPP_C_mmol_m2_d'])+' mmolC/m²/d'],['Export C a200m',n(r['bottom_export_C_mg_m2_d'])+' mgC/m²/d'],['DIC minimo campionato',n(r['last_year_min_DIC'])+' mmolC/m³'],['DIN massimo campionato',n(r['last_year_max_DIN'])+' mmolN/m³'],['Apporto C netto dal reservoir',n(balance('R00_buffered_reference').get('C_external_net_mmol_m2_d'))+' mmolC/m²/d'],['Accumulo N nel sistema',n(balance('R00_buffered_reference').get('N_inventory_accumulation_mmol_m2_d'))+' mmolN/m²/d'],['Deriva P_N tra ultimi intervalli annuali',n(100*r['last_year_N_drift_fraction'])+'%']],[W*.66,W*.34])
p('L’apporto N supera le uscite e il DIN continua ad accumularsi. La quasi-stazionarietà della biomassa non è equilibrio dell’intero sistema. Il reservoir C è un controllo volumetrico ideale; non rappresenta uno scambio gassoso alla superficie.',small=True)
p('I confronti successivi partono tutti dallo stesso stato fisico finale R00 e durano365giorni. I contatori vengono azzerati. Il raddoppio della griglia verticale prolunga conservativamente le concentrazioni;16→32 classi di particelle conserva C e N mediante mappatura logaritmica.')
rows=[['Prova annuale','Δ P_N','Δ C:N','Δ NPP','Δ export','Δ ESD','Esito']]
for g in ASSESS.get('recommended_numerical_gates',[]):
 changes=g['relative_changes'];label={'time_step':'Passo temporale','vertical_grid':'40 vs80 celle','vertical_grid80_to160':'80 vs160 celle','particle_grid':'Classi detritiche'}[g['test']]
 rows.append([label,*[pct(changes[k]) for k in ['phyto_N_mmol_m2','phyto_CN_inventory_ratio','NPP_C_mmol_m2_d','bottom_export_C_mmol_m2_d','N_weighted_geomean_ESD_um']], 'passa' if g['pass_5pct_bulk_10pct_size'] else 'non passa'])
tbl(rows,[W*.22,*[W*.125]*5,W*.155])
p('Scarti rispetto al riferimento più fine; soglie euristiche dichiarate:5% per stock N,C:N,NPP ed export,10% per taglia media geometrica. Non sono errori osservativi. Un solo confronto annuale non dimostra convergenza asintotica, né trasferisce automaticamente l’esito alla traiettoria decennale.',small=True)
if ASSESS:p('<b>Giudizio numerico:</b> '+('tutte le soglie annuali sono rispettate.' if ASSESS.get('recommended_bulk_convergence_demonstrated') else 'la convergenza complessiva della configurazione riveduta non è dimostrata. I confronti che falliscono sono mantenuti e limitano le inferenze quantitative.'))
page();h('15. Export e luce: interventi controllati')
fig('12_annual_interventions',height=260)
p('Sette interventi annuali riprendono le domande di DePasquale4.5 e5: apporto N×0.25/2, lisi e remineralizzazione×2, sinking del fitoplancton×2, minore resa C dei consumatori e attenuazione dell’acqua0.015/0.06m⁻¹. Il controllo è R01 con identico stato iniziale, passo e finestra annuale.')
p('I cambiamenti includono risposte transitorie della comunità, non nuovi equilibri decennali. NPP ed export possono avere segno e ampiezza differenti: la quantità di carbonio fissato non determina da sola quanto attraversa il fondo. Nessun intervento convalida la legge di sinking o la remineralizzazione contro osservazioni oceaniche.')
p('Gli interventi sull’apporto N partono da un grande inventario DIN ereditato dal riferimento. Una risposta biologica debole nel corso di un anno riguarda quel regime già ricco di N; le serie inizializzate separatamente e le colture in carenza esplorano condizioni diverse.',small=True)
if S.get('annual_intervention_comparisons'):
 rows=[['Intervento','Δ NPP','Δ export C']]
 for a in S['annual_intervention_comparisons']:
  rows.append([a['alternative'].replace('R04_','').replace('R05_',''),pct(a['relative_changes']['NPP_C_mmol_m2_d']),pct(a['relative_changes']['bottom_export_C_mmol_m2_d'])])
 tbl(rows,[W*.56,W*.22,W*.22])
#14
page();h('16. Confronto critico con la letteratura recente')
tbl([['Studio primario','Implicazione per il modello'],['Omta2024 -10.1016/j.jtbi.2024.111883','Accumulo e allocazione dipendono anche dal gruppo; i parametri di storage sono incerti. Le stesse colture non forniscono validazione indipendente.'],['Harcourt2024 -10.1371/journal.pone.0292337','Differenze tra ceppi e adattamento termico limitano una fisiologia universale basata sulla sola taglia.'],['Mollica2025 -10.1002/lno.70058','Le piccole frazioni possono dominare la stechiometria della comunità: occorrono pesi e pool confrontabili.'],['Karakuş2025 -10.1029/2024GB008403','Flusso di export e trasferimento verso il mare profondo rispondono a processi diversi.'],['Dogliotti2026 -10.5194/bg-23-1859-2026','Luce e mixing devono essere misurati o prescritti indipendentemente per interpretare la fenologia.'],['Rogerson2026 -10.5194/bg-23-3179-2026','La complessità delle quote modifica NPP ed export in modo non uniforme: servono confronti controllati.'],['Papapostolou2026 -10.5194/gmd-19-6627-2026','Taglia e costi metabolici sono utili, ma mixotrofia,ontogenesi e altri percorsi possono cambiare le reti trofiche.'],['Sharoni2026 -10.1038/s41558-026-02598-w','Allocazione intracellulare e redistribuzione verticale possono produrre risposte diverse del C:N integrato.'],['Bernish2026 -10.3389/feart.2026.1839921','L’interazione Fe-luce-temperatura non è sostituibile con il solo parametro di saturazione DIN.']],[W*.45,W*.55])
p('Questa è una comparazione mirata con studi primari, verificata al16/10/2026; non una revisione sistematica. Il supporto teorico a un meccanismo non equivale alla sua validazione nel nuovo SPLAS.',small=True)
#15
page();h('17. Configurazione operativa e decisioni')
p('<b>Utilizzabile adesso:</b> il codice per esperimenti di meccanismo, controlli di bilancio, confronto di ipotesi e valutazione fisiologica circoscritta. Ogni esecuzione registra metodo,passo,parametri,forzanti,finestre e hash del nucleo.')
p('<b>Da non considerare dimostrato:</b> validità universale delle allometrie e del Q10; fit entro l’incertezza delle repliche; convergenza di tutte le metriche; risposta realistica dell’export o previsione di bloom per un sito.')
tbl([['Criterio','Regola operativa'],['C:N e unità','Usare inventari molari espliciti; conservare carbonio anche nel pool organico di N.'],['Integratore','Riferimento esplicito per casi piccoli; IMEX solo con raffinamento del passo e riscontri per la metrica di interesse.'],['Griglia verticale','La configurazione40celle fallisce il confronto con80; consultare anche il controllo80→160 prima di usare export quantitativo.'],['Risoluzione','Conservare inventari iniziali e profilo fisico DIN sotto raffinamento; normalizzare spettri in Δlog ESD.'],['Particelle','Scegliere classi dopo il confronto2/4/8/16; la chiusura della massa non certifica export convergente.'],['Diagnostici','TTE di flusso e Z:P separati; e-ratio alla profondità dichiarata; soglie critiche integrate, Zeu e mixing distinti.'],['Uso empirico','Stimare su dati dedicati e verificare su un nuovo dataset non usato nello sviluppo, con forzanti e pool confrontabili.']],[W*.29,W*.71])
p('Priorità scientifiche: un controllo a rapporto fisso con cinetiche abbinate; dati indipendenti su quota e crescita; diagnostica di N-speciation e dei costi energetici; successivamente P,Fe,Si e processi delle particelle se necessari alla domanda. Una componente aggiuntiva va giustificata dalla capacità predittiva, non dalla sola disponibilità di nuovi stati.')
#16
page();h('18. Riproducibilità e materiali')
p('Repository: https://github.com/Marcofianchini/SPLAS-1.1<br/>Ramo: codex/dynamic-cn-development<br/>Proposta: https://github.com/Marcofianchini/SPLAS-1.1/pull/1<br/>Base storica:76f992c9924420c715a660d017387882efe992dd')
p('Impronta SHA256 del nucleo verificato:<br/>'+S['model_sha256'],small=True)
p('Installazione e prove:<br/><font face="DV">python -m pip install -c research/requirements-reference.txt -e &quot;.[fast,research]&quot;<br/>python -m unittest discover -s tests -v<br/>python -m splas2 --days 30 --output results/example.npz</font>',small=True)
p('Sequenza riproducibile: prepare_cultures.py → calibrate_cultures.py → run_experiments.py → run_mechanisms.py → run_physiological_protocols.py → run_supplement.py → run_carbon_boundary.py → run_recommended.py → run_recommended_refinement.py → run_export_interventions.py → analyze_results.py → assess_validation.py → make_validation_report.py → package_report_data.py. Il PDF è generato dai JSON delle esecuzioni.')
p('Il PDF distribuito contiene un allegato ZIP con i risultati numerici e un manifesto SHA256. Per estrarlo: python scripts/extract_report_results.py rapporto.pdf --output results. I sorgenti e le configurazioni versionate restano nel repository.',small=True)
p('Materiali principali: splas2/model.py;test_conservation.py;research/MODEL.md;LITERATURE.md;experiment_registry.json;supplement_registry.json;experiment_summary.csv;validation_summary.json;reference_environment.json. I dati grezzi CC-BY4.0 e la loro estrazione sono identificati da DOI e hash.')
refs=[('Moloney1988','https://hdl.handle.net/11427/22132'),('DePasquale2018','https://arts.units.it/bitstream/11368/2929831/2/M_DePasquale_TESI_FINAL.pdf'),('Fianchini2021','PDF corretto fornito dall’utente'),('Montagnes1994','https://doi.org/10.4319/lo.1994.39.5.1044'),('Liefer2019','https://doi.org/10.3389/fmicb.2019.00763'),('Liefer etal2026 dataset','https://doi.org/10.1594/PANGAEA.900248'),('Hagstrom2024','https://doi.org/10.1029/2023GB007991')]
for name,url in refs:p('<b>'+name+'</b> - '+url,small=True)
p('Il file LITERATURE.md contiene i riferimenti completi dei restanti studi e una matrice evidenza→decisione→limite. I PDF delle tesi e gli articoli non sono redistribuiti come parte del codice.')
if 'R02_grid160' in CASES:
 page();h('Appendice. Il raffinamento verticale adattivo')
 fig('13_nested_vertical_refinement',height=310)
 p('Il mancato superamento del criterio40→80 celle ha motivato un ulteriore confronto a160 celle. Le tre simulazioni annuali partono dallo stesso stato fisico R00, prolungato conservativamente. Ogni confronto mantiene parametri, forzanti, intervallo e passo massimo0.0125d.')
 p('Il confronto aggiuntivo è adattivo, dichiarato dopo il fallimento40→80; le soglie5%/10% sono rimaste invariate. Tre punti di una dinamica ecologica non lineare non dimostrano da soli un ordine asintotico di convergenza.')
 for g in ASSESS.get('recommended_numerical_gates',[]):
  if g['test']=='vertical_grid80_to160':
   d=g['relative_changes'];p('Rispetto alle160 celle, le80 celle danno scarti di '+pct(d['phyto_N_mmol_m2'])+' per P_N, '+pct(d['NPP_C_mmol_m2_d'])+' per NPP e '+pct(d['bottom_export_C_mmol_m2_d'])+' per export C. Il miglior accordo dell’export non basta a superare il criterio complessivo sugli stock e sui flussi.')
 p('La prima formulazione usa advezione upwind del primo ordine e luce valutata al centro delle celle con clorofilla diagnostica. Diffusione numerica, discretizzazione della luce e riassetto biologico sono contributi da discriminare. Il bilancio atomico può restare chiuso durante una risposta sensibile alla griglia.')
 p('L’export quantitativo richiede una verifica spaziale dedicata. Un futuro schema di trasporto di ordine superiore va sottoposto agli stessi controlli indipendenti di positività, conservazione e confronto con soluzioni analitiche; cambiare schema o aumentare la risoluzione non convalida le cinetiche delle particelle.',small=True)
OUTPUT.parent.mkdir(parents=True,exist_ok=True)
doc=SimpleDocTemplate(str(OUTPUT),pagesize=A4,rightMargin=18*2.83465,leftMargin=18*2.83465,topMargin=18*2.83465,bottomMargin=23*2.83465,title='SPLAS-CN: sviluppo e verifica',author='Codex - analisi per il progetto SPLAS')
doc.build(story,onFirstPage=footer,onLaterPages=footer)
print(OUTPUT)
