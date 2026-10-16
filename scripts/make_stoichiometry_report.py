"""Italian research notebook: results, counterexamples, process plots and limits."""
from pathlib import Path
import hashlib,html,json,re
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,Image
from pypdf import PdfReader,PdfWriter

ROOT=Path(__file__).resolve().parents[1];DOC=ROOT/'research/stoichiometry';FIG=ROOT/'results/stoichiometry'
OUT=ROOT/'output/pdf/SPLAS_Patterns_Stechiometria_2026-10-16.pdf';W=A4[0]-92
NAVY=colors.HexColor('#163746');TEAL=colors.HexColor('#167d9a');MUTED=colors.HexColor('#526771')
pdfmetrics.registerFont(TTFont('DV','/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
pdfmetrics.registerFont(TTFont('DVB','/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'))
pdfmetrics.registerFontFamily('DV',normal='DV',bold='DVB',italic='DV',boldItalic='DVB')
ST={
 'body':ParagraphStyle('body',fontName='DV',fontSize=10.3,leading=14.5,spaceAfter=8,textColor=NAVY),
 'small':ParagraphStyle('small',fontName='DV',fontSize=8.6,leading=11.7,spaceAfter=7,textColor=MUTED),
 'cell':ParagraphStyle('cell',fontName='DV',fontSize=8.6,leading=11.7,textColor=NAVY),
 'head':ParagraphStyle('head',fontName='DVB',fontSize=8.6,leading=11.7,textColor=colors.white),
 'h1':ParagraphStyle('h1',fontName='DVB',fontSize=22,leading=27,spaceAfter=13,textColor=NAVY),
 'h2':ParagraphStyle('h2',fontName='DVB',fontSize=12,leading=16,spaceBefore=8,spaceAfter=8,textColor=TEAL),
 'cover':ParagraphStyle('cover',fontName='DVB',fontSize=32,leading=38,spaceAfter=18,textColor=NAVY),
 'kicker':ParagraphStyle('kicker',fontName='DVB',fontSize=8.3,leading=12,spaceAfter=9,textColor=TEAL)}
story=[]

def clean(s):
    # Preserve links, identifiers and fingerprints while spacing scientific prose.
    chunks=re.split(r'(<[^>]+>|&(?:[A-Za-z]+|#x[0-9A-Fa-f]+|#\d+);|\b[0-9a-f]{40,64}\b|\S*_\S*|\b[01]D\b)',str(s))
    for i in range(0,len(chunks),2):
        chunks[i]=re.sub(r'([A-Za-zÀ-ÿ])(\d)',r'\1 \2',chunks[i])
        chunks[i]=re.sub(r'(\d)([A-Za-zÀ-ÿµ])',r'\1 \2',chunks[i])
        chunks[i]=re.sub(r';(?=\S)','; ',chunks[i])
    return ''.join(chunks).replace('assimilazioneN','assimilazione N').replace('produzioneZ','produzione Z').replace('ingestioneC','ingestione C').replace('vedere<link','vedere <link')

def inline(s):
    s=html.escape(s)
    s=re.sub(r'\[([^\]]+)\]\(([^)]+)\)',r'<link href="\2" color="#167d9a">\1</link>',s)
    s=re.sub(r'\*\*(.+?)\*\*',r'<b>\1</b>',s)
    s=re.sub(r'`([^`]+)`',r'<font name="Courier" size="8.7">\1</font>',s)
    return clean(s)

def p(s,kind='body',raw=False):story.append(Paragraph(clean(s) if raw else inline(s),ST[kind]))
def page(title,kicker):
    if story:story.append(PageBreak())
    p(kicker.upper(),'kicker');p(title,'h1')
def table(rows,widths=None):
    widths=widths or [W/len(rows[0])]*len(rows[0])
    data=[[Paragraph(inline(str(v)),ST['head' if i==0 else 'cell']) for v in row] for i,row in enumerate(rows)]
    t=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT')
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),NAVY),('VALIGN',(0,0),(-1,-1),'TOP'),
        ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#edf4f5')]),
        ('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7),
        ('TOPPADDING',(0,0),(-1,-1),7),('BOTTOMPADDING',(0,0),(-1,-1),7),
        ('LINEBELOW',(0,0),(-1,0),.6,NAVY),('LINEBELOW',(0,1),(-1,-1),.3,colors.HexColor('#d8e4e7'))]))
    story.extend([t,Spacer(1,9)])
def figure(name,caption):
    img=Image(str(FIG/(name+'.png')));img.drawHeight=img.imageHeight*W/img.imageWidth;img.drawWidth=W
    story.extend([img,Spacer(1,6)]);p(caption,'small')

text=(DOC/'RESULTS.md').read_text()
sections={}
for chunk in re.split(r'^##\s*',text,flags=re.M)[1:]:
    title,body=chunk.split('\n',1);sections[title.strip()]=body.strip()

def section(prefix):
    body=next(v for k,v in sections.items() if k.startswith(prefix))
    blocks=re.split(r'\n\s*\n',body)
    for block in blocks:
        if block.startswith('|'):
            rows=[[v.strip() for v in line.strip().strip('|').split('|')] for line in block.splitlines() if not re.match(r'^\|[-:| ]+\|$',line)]
            table(rows)
        elif block.startswith('```'):
            for line in block.splitlines()[1:-1]:p(line,'small')
        elif block.startswith('`produzione_Z_C'):
            p('g<sub>Z,C</sub> = min(0.33 I<sub>C</sub>, 5.5 a<sub>N</sub> I<sub>N</sub>)','h2',raw=True)
            p('I flussi I indicano il materiale ingerito; aN è l’efficienza di assimilazione dell’azoto.','small')
        else:p(' '.join(block.splitlines()))

data=json.loads((DOC/'summary.json').read_text());by={r['id']:r for r in data['cases']};audit=data['audit']
archive=json.loads((DOC/'data_manifest.json').read_text())
p('Quaderno di ricerca /16 ottobre2026','kicker')
p('SPLAS<br/>stechiometria<br/>e pattern emergenti','cover',raw=True)
p('Dalle riserve cellulari alla dieta e all’organizzazione della comunità','h2')
p('166 configurazioni ideali0D, senza nuova calibrazione.142 casi pianificati prima delle integrazioni e24 controlli successivi, registrati separatamente.')
section('Risultato centrale')
table([['Esito','Interpretazione POM'],['Più C non implica più produzione Z','Quantità, qualità e riciclo possono compensarsi.'],
 ['La dieta differisce dalla biomassa media','Il kernel dimensionale filtra la composizione elementare disponibile.'],
 ['Utilizzo N dei consumatori → struttura delle taglie','Contrazione di taglia ed evenness conservata nei confronti appaiati.'],
 ['Memoria temporale non dimostrata','Segni e ampiezze dipendono da passo, fase iniziale ed età dell’impulso.']],[W*.36,W*.64])
p('Studio orientato ai pattern, secondo la progressione De Pasquale → Fianchini. I valori numerici definiscono i regimi del toy world; non sono una previsione di un ecosistema localizzato.','small')

page('Dal protocollo agli esperimenti','01 / Disegno e confrontabilità')
section('Disegno e valuta')
table([['Blocco','Casi'],['Piano minimo1P+1Z e senza Z',80],['Comunità: N-assimilazione, phi uniforme e senza Z',36],
 ['Assemblaggio chiuso con inventari uguali',12],['Risoluzione temporale e dimensionale pianificata',9],['Sequenze di impulsi pianificate',5],
 ['Controlli fisiologia, accumulo e passo successivi',16],['Controlli appaiati e sequenze a passo fine',8]],[W*.83,W*.17])
p('I24 controlli successivi sono esplorativi. Il disegno iniziale è registrato prima delle nuove integrazioni, senza una preregistrazione pubblica esterna. Tutti gli esiti rimangono nel dataset.','small')

page('Quantità e qualità del cibo','02 / Caso minimo leggibile')
figure('02_quantity_quality_feedback','DIN entrante1. La resa diminuisce con il C:N, mentre la produzione Z dipende anche da disponibilità, ingestione e riciclo. La linea orizzontale indica la resa C-limite prescritta0.33.')
section('1. Più C')

page('Il gradiente luce × azoto','03 / Il regime cambia la risposta')
figure('01_minimal_resource_plane','Riferimento1P+1Z, assimilazioneN0.67. Medie degli ultimi365 giorni. La luce aumenta il C:N in ogni colonna, ma la biomassa Z non risponde nello stesso modo a tutti i livelli di N.')
p('A DIN entrante0.25, più luce aumenta PC di circa19% e diminuisce lievemente la produzione Z. A DIN1, PC aumenta circa38% e produzione Z rimane quasi invariata. A DIN4, produzione Z cresce: il vincolo stechiometrico non basta a imporre un paradosso energetico generalizzato.')
p('La plausibilità del disaccoppiamento quantità/qualità è confrontata con Goldenberg et al.(2024); non si calibra la soglia di resa sui loro mesocosmi. La diminuzione di resa prevista dalla formula è un vincolo, mentre il regime raggiunto dalla rete è il risultato da spiegare.','small')

page('Biomassa presente e dieta ingerita','04 / Filtraggio trofico emergente')
figure('08_diet_filtering','Frazioni di C presente e C ingerito per classe; linea viola: C:N delle classi rilevanti. Il caso povero di N separa fortemente il serbatoio di C dalla frazione più consumata.')
section('2. La comunità')

page('L’uso dell’N riorganizza la comunità','05 / Funzione e struttura insieme')
figure('03_community_patterns','Linee continue: assimilazioneN0.67; tratteggiate:0.335. Le classi restano tratti dimensionali. Evenness: distribuzione della biomassa su14 intervalli logaritmici fisici comuni.')
section('3. L\'uso')

page('Che cosa sostiene la coesistenza?','06 / Controlli sulla rete')
figure('05_foodweb_controls','Spettri in N normalizzati. Selettività, rendimento del consumatore e assenza di predazione producono organizzazioni diverse dalle stesse quantità iniziali totali di C e N.')
p('Senza predatori vince la classe1µm in tutti i nove ambienti. La stechiometria flessibile da sola non sostiene una comunità ricca di taglie. Il grazing selettivo modifica il vantaggio dei competitori più piccoli e apre spazio alle classi maggiori.')
p('Con phi=1 scompaiono le preferenze dimensionali: cambia anche la quantità di preda incontrata, mentre la massima ingestione per consumatore rimane uguale. Questo è un confronto strutturale dell’incontro, non un’ablazione isolata della qualità alimentare.')
p('L’evenness sulla griglia comune misura la distribuzione per dimensione, non la diversità tassonomica. La ricostruzione conservativa di una classe larga può distribuirla su due intervalli comuni: un valore positivo non prova la presenza di due specie.','small')

page('Assemblaggio a inventari corretti','07 / Riprendere De Pasquale')
figure('04_matched_assembly','Sistemi chiusi: C e N totali uguali fra casi. Le curve riportano classi con biomassa rilevante; l’aumento del C:N aggregato combina riserve e selezione della taglia vincente.')
section('4. Assemblaggio')

page('Plasticità e composizione si oppongono','08 / Riprendere la domanda di Fianchini')
section('5. Plasticità')
figure('07_size_succession','Storia dimensionale completa: frazione di PN in log10, a luce e ingresso N costanti. La successione del transiente e la banda dominante tardiva sono distinte dalla sola media finale.')

page('Oscillazioni senza stagionalità imposta','09 / Un pattern, non una previsione di fase')
figure('09_autonomous_variability','Regime ricco di N, stesse condizioni esterne. L’esistenza della variabilità persiste al raffinamento; traiettorie e composizione dettagliata differiscono.')
section('6. Variabilità')

page('Sequenze e memoria: esito non risolto','10 / Conservare il controesempio')
figure('06_sequence_memory','Impulso N e tre giorni di luce alta in ordine invertito, più controlli. Linea verticale: ritorno allo stesso ambiente esterno. Le differenze fra traiettorie non identificano automaticamente una memoria delle riserve.')
section('7. La prova')

page('Che cosa resiste ai controlli?','11 / Robustezza delle inferenze')
table([['Contrasto appaiato','Δ taglia N','Δ evenness','Δ produzione Z'],
 *[[x['lower_eta'].replace('B_',''),f"{x['delta_geomean_percent']:.1f}%",f"{x['delta_evenness']:.3f}",f"{x['delta_Zprod_percent']:.1f}%"] for x in data['contrasts']['paired_eta_robustness']]],[W*.43,W*.19,W*.19,W*.19])
p('Confronto assimilazioneN0.335 contro0.67 a identico passo/griglia. Il segno negativo congiunto resta riconoscibile nei sei controlli. La convergenza quantitativa dell’intero modello non è dimostrata da questi confronti.','small')
section('Verifica, portata')

page('Valutazione POM dei pattern','12 / Regole imposte e organizzazione emergente')
patterns=json.loads((DOC/'PATTERNS.json').read_text())['patterns']
italian=[('Quantità e qualità del cibo','Riconoscibile, dipende dal regime'),
    ('Filtraggio trofico della stechiometria','Riconoscibile con tratti e incontri dichiarati'),
    ('Uso dell’N e struttura della comunità','Firma qualitativa appaiata meglio sostenuta'),
    ('Assemblaggio e rilascio competitivo','Riconoscibile a inventari iniziali uguali'),
    ('Plasticità e sostituzione fra classi','Esempio informativo; segno aggregato non risolto'),
    ('Variabilità interna persistente','Riconoscibile; dinamica dettagliata non convergente'),
    ('Memoria stechiometrica di sequenza','Non dimostrata dal disegno attuale')]
table([['Pattern','Stato della spiegazione'],*[[x['id']+' — '+label,status] for x,(label,status) in zip(patterns,italian)]],[W*.46,W*.54])
p('Non si sommano sette conferme indipendenti. Resa e C:N della dieta sono legati algebricamente; spettri e loro medie non costituiscono osservazioni indipendenti; più simulazioni deterministiche non sono repliche empiriche.')
p('La forza esplicativa viene dalla congiunzione: una sola ipotesi di utilizzo N influenza produzione dei consumatori, rilascio competitivo e struttura delle taglie. Le alternative fisiologiche mostrano che la firma C:N media è meno generale di quella dimensionale.')
p('La parsimonia qui consiste nel rendere discriminabili le retroazioni esistenti. Prima di aggiungere preferenza comportamentale per qualità o nuovi elementi, occorre stabilire quali pattern il solo filtro dimensionale riesce a spiegare e dove fallisce.')

page('Confronto con letteratura corrente','13 / Plausibilità, alternative e portata')
refs=[
 ('Liefer et al.,2019','10.3389/fmicb.2019.00763','Riserve e perdita di materiale ricco in N motivano C:N flessibile. Non identificano un’unica allometria. Il fit Ot/Tw usa parte degli stessi dati: non è una conferma indipendente.'),
 ('Goldenberg et al.,2024','10.1038/s41598-024-54345-w','Mesocosmi: più biomassa/produttività non garantisce un buon trasferimento ai grazers. Pertinente alla distinzione quantità/qualità; non rende le classi SPLAS equivalenti alle diatomee osservate.'),
 ('Neun et al.,2025','10.1111/1365-2745.70161','Luce, nutrienti e composizione interagiscono nella qualità delle risorse. Spettro luminoso, fosforo e acidi grassi sono spiegazioni reali aggiuntive, assenti in questa campagna C/N.'),
 ('Matsumoto et al.,2026','10.1016/j.ecolmodel.2025.111423','Modello0D/1D: scelta per qualità modifica comunità e flussi. In SPLAS la preferenza resta dimensionale: è una spiegazione alternativa da confrontare, non un meccanismo già implementato.'),
 ('Sharoni et al.,2026','10.1038/s41558-026-02598-w','Allocazione macromolecolare e distribuzione delle popolazioni cambiano con l’ambiente. Motiva la distinzione plasticità/composizione; il toy case non verifica proiezioni climatiche.'),
 ('Diehl,2007','10.1086/516655','La teoria ammette risposte diverse all’arricchimento luminoso e minerale. SPLAS non deve essere forzato a produrre collasso, isteresi o multistabilità in ogni regime.')]
for name,doi,comment in refs:
    p(f'<b>{name}</b> — <link href="https://doi.org/{doi}" color="#167d9a">{doi}</link>','h2',raw=True);p(comment)
p('Fonti primarie consultate al16 ottobre2026. Per Matsumoto e Sharoni: abstract e sezioni/testo indicizzato disponibili. Titoli completi e distinzione fra evidenza sperimentale e risultato di modello sono nel documento LITERATURE.md.','small')

page('Riprodurre e sviluppare lo studio','14 / Codice, stati e provenienza')
section('Riproducibilità')
p('Dataset: '+archive['file_name'],'h2')
p(f"{archive['cases']} casi; {archive['bytes']/1024**2:.1f}MiB. Tutti gli stati, le griglie e i contatori integrati dei flussi sono conservati esattamente; gli osservabili istantanei si rigenerano dai sorgenti corrispondenti.")
p('SHA256 dataset: '+archive['sha256'],'small')
p('Modello numerico invariato: '+hashlib.sha256((ROOT/'splas2/model.py').read_bytes()).hexdigest(),'small')
p('<link href="https://github.com/Marcofianchini/SPLAS-1.1/pull/1" color="#167d9a">Codice e revisione di sviluppo — GitHub PR1</link>','body',raw=True)
p('<link href="https://github.com/Marcofianchini/SPLAS-1.1/tree/codex/dynamic-cn-development/research/stoichiometry" color="#167d9a">Registri, sintesi numerica e valutazione dei pattern</link>','body',raw=True)
p('Allegati a questo PDF: tabella completa dei risultati, sintesi JSON, valutazione POM e impronta del dataset. Il dataset numerico è fornito separatamente con il rapporto.','small')

def footer(canvas,doc):
    canvas.setFillColor(TEAL);canvas.rect(46,A4[1]-31,30,3,fill=1,stroke=0)
    canvas.setFont('DV',8);canvas.setFillColor(MUTED);canvas.drawString(83,A4[1]-31,'SPLAS | Stechiometria e pattern emergenti')
    canvas.setStrokeColor(colors.HexColor('#d8e4e7'));canvas.line(46,38,A4[0]-46,38)
    canvas.setFont('DV',7.5);canvas.drawString(46,25,'16 ottobre 2026 | Esplorazione POM in un mondo ecologico ideale')
    canvas.drawRightString(A4[0]-46,25,str(doc.page))

OUT.parent.mkdir(parents=True,exist_ok=True);temporary=ROOT/'tmp/stoichiometry-report-base.pdf';temporary.parent.mkdir(exist_ok=True)
SimpleDocTemplate(str(temporary),pagesize=A4,leftMargin=46,rightMargin=46,topMargin=52,bottomMargin=51,
    title='SPLAS — Stechiometria e pattern emergenti',author='SPLAS research notebook').build(story,onFirstPage=footer,onLaterPages=footer)
reader=PdfReader(temporary);writer=PdfWriter();writer.clone_document_from_reader(reader)
for name in ['summary.csv','summary.json','PATTERNS.json','data_manifest.json']:writer.add_attachment(name,(DOC/name).read_bytes())
with OUT.open('wb') as f:writer.write(f)
manifest=dict(file_name=OUT.name,bytes=OUT.stat().st_size,pages=len(reader.pages),sha256=hashlib.sha256(OUT.read_bytes()).hexdigest(),
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),data_sha256=archive['sha256'])
(DOC/'report_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest,indent=2))
