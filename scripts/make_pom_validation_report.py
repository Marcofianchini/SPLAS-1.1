"""Critical Italian POM research report; exact machine-readable attachments."""
from pathlib import Path
import hashlib,html,json,re
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,Image
from pypdf import PdfReader,PdfWriter

ROOT=Path(__file__).resolve().parents[1];DOC=ROOT/'research/pom_validation_0d';FIG=ROOT/'results/pom_validation_0d'
OUT=ROOT/'output/pdf/SPLAS_Validazione_POM_0D_2026-10-16.pdf';W=A4[0]-88
NAVY=colors.HexColor('#173c47');TEAL=colors.HexColor('#177c8c');MUTED=colors.HexColor('#546a72')
pdfmetrics.registerFont(TTFont('PV','/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
pdfmetrics.registerFont(TTFont('PVB','/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'))
pdfmetrics.registerFontFamily('PV',normal='PV',bold='PVB',italic='PV',boldItalic='PVB')
ST={
 'body':ParagraphStyle('body',fontName='PV',fontSize=9.5,leading=13.3,spaceAfter=7.5,textColor=NAVY),
 'small':ParagraphStyle('small',fontName='PV',fontSize=8,leading=11,spaceAfter=6,textColor=MUTED),
 'cell':ParagraphStyle('cell',fontName='PV',fontSize=8.3,leading=11.5,textColor=NAVY),
 'head':ParagraphStyle('head',fontName='PVB',fontSize=8.3,leading=11.5,textColor=colors.white),
 'h1':ParagraphStyle('h1',fontName='PVB',fontSize=21,leading=26,spaceAfter=13,textColor=NAVY),
 'h2':ParagraphStyle('h2',fontName='PVB',fontSize=11.3,leading=15,spaceBefore=6,spaceAfter=7,textColor=TEAL),
 'cover':ParagraphStyle('cover',fontName='PVB',fontSize=31,leading=37,spaceAfter=18,textColor=NAVY),
 'kicker':ParagraphStyle('kicker',fontName='PVB',fontSize=8,leading=11,spaceAfter=9,textColor=TEAL)}
story=[]


def clean(s):
    chunks=re.split(r'(<[^>]+>|&(?:[A-Za-z]+|#x[0-9A-Fa-f]+|#\d+);|\b[0-9a-f]{40,64}\b|\S*_\S*|\b[01]D\b)',str(s))
    for i in range(0,len(chunks),2):
        chunks[i]=re.sub(r'([A-Za-zÀ-ÿ])(\d)',r'\1 \2',chunks[i])
        chunks[i]=re.sub(r'(\d)([A-Za-zÀ-ÿµ])',r'\1 \2',chunks[i])
    result=''.join(chunks).replace('−','-').replace('–','-').replace('—','-').replace('\u2011','-')
    result=re.sub(r'×\s*10-(\d+)',r'×10<super>-\1</super>',result)
    return re.sub(r'\b([md])-(\d+)\b',r'\1<super>-\2</super>',result)


def inline(s):
    s=html.escape(str(s));s=re.sub(r'\[([^\]]+)\]\(([^)]+)\)',r'<link href="\2" color="#177c8c">\1</link>',s)
    s=re.sub(r'\*\*(.+?)\*\*',r'<b>\1</b>',s);s=re.sub(r'`([^`]+)`',r'<font name="Courier" size="7">\1</font>',s)
    return clean(s)


def p(s,kind='body',raw=False):story.append(Paragraph(clean(s) if raw else inline(s),ST[kind]))
def page(title,kicker):
    if story:story.append(PageBreak())
    p(kicker.upper(),'kicker');p(title,'h1')
def table(rows,widths=None):
    data=[[Paragraph(inline(v),ST['head' if i==0 else 'cell']) for v in row] for i,row in enumerate(rows)]
    t=Table(data,colWidths=widths or [W/len(rows[0])]*len(rows[0]),repeatRows=1,hAlign='LEFT')
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),NAVY),('VALIGN',(0,0),(-1,-1),'TOP'),
        ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#eff5f4')]),
        ('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7),
        ('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),
        ('LINEBELOW',(0,1),(-1,-1),.3,colors.HexColor('#d3dfdf'))]))
    story.extend([t,Spacer(1,8)])
def figure(name,caption,maxheight=320):
    img=Image(str(FIG/(name+'.png')));scale=min(W/img.imageWidth,maxheight/img.imageHeight)
    img.drawWidth=img.imageWidth*scale;img.drawHeight=img.imageHeight*scale
    story.extend([img,Spacer(1,5)]);p(caption,'small')


text=(DOC/'RESULTS.md').read_text();sections={}
for chunk in re.split(r'^##\s*',text,flags=re.M)[1:]:
    title,body=chunk.split('\n',1);sections[title.strip()]=body.strip()
def blocks(prefix):
    body=next(v for k,v in sections.items() if k.startswith(prefix))
    return [' '.join(b.splitlines()) for b in re.split(r'\n\s*\n',body) if not b.startswith('|')]
def section(prefix,selection=None):
    bb=blocks(prefix)
    for b in bb if selection is None else bb[selection]:p(b)


s=json.loads((DOC/'summary.json').read_text());by={r['id']:r for r in s['rows']}
archive=json.loads((DOC/'data_manifest.json').read_text());verification=json.loads((DOC/'VERIFICATION.json').read_text())
p('Quaderno ecologico / campagna avviata il16 ottobre2026','kicker')
p('SPLAS<br/>Validazione POM<br/>dello0D','cover',raw=True)
p('Esito: supporto parziale, con discrepanze localizzate','h2')
p('221 nuovi run, criteri espliciti, confronti fra ipotesi e verifiche successive. Nessun nuovo fit sui pattern; equazioni e integratore congelati.')
table([['Bersaglio','Verdetto dopo il controllo dell’intero piano a passo fine'],
    ['Arricchimento e taglie','Compatibile: aumentano biomassa e taglia, con persistenza di taglie piccole.'],
    ['Gruppi dimensionali','Compatibile nella finestra studiata; non dimostra resistenza all’invasione.'],
    ['Quantità e qualità','Presente in un regime; non supera il criterio di generalità dichiarato.'],
    ['Diversità e Z:P','Firma congiunta non riprodotta. Il confronto diretto con i gruppi del campo rimane insufficiente.']],[W*.29,W*.71])
p('Il risultato decisivo non è il numero dei run: una comunità quasi stazionaria può essere riorganizzata dalla disponibilità di tratti rarissimi, senza aggiungere C o N. La memoria ha una componente nelle riserve, ma dipende dallo stato ecologico.')
p('Il modello rimane utile per domande esplicative. Questa campagna non autorizza a definirlo generalmente validato né a usare il passaggio allo spazio come soluzione implicita delle discrepanze.','small')

page('Una matrice congiunta, senza punteggio','01 / Decisione sui pattern')
figure('01_POM_matrix','Sì/No applica le soglie del protocollo. I numeri indicano quanti dei5 blocchi annuali soddisfano il criterio; non sono repliche indipendenti. L’ultima riga è il controllo dell’intero piano full a h0.00625.',maxheight=310)
section('Matrice congiunta')

page('Domande prima delle integrazioni','02 / Disegno e fonti')
section('Che cosa')
table([['Pattern','Evidenza esterna','Trasferimento ammesso al toy case'],
    ['P1','Matthiessen et al.2025','Sostituzione delle classi; non evoluzione o volume plastico.'],
    ['P2','Vascotto et al.2024','Gruppi dimensionali; ESD/N non equivale a MLD/conteggi.'],
    ['P3','Goldenberg et al.2024','Disaccoppiamento di rete; driver luce/N diverso da upwelling/Si:N.'],
    ['P4','García-Comas et al.2016','Sfida storica delle tesi; equità della biomassa diversa da entropia individuale.']],[W*.10,W*.29,W*.61])
p('Le5 alternative sono full, mortalità Z per classe, capacità C senza esponente, incontri uniformi e fisiologia a priori; il controllo senza Z verifica la competizione. Le soglie numeriche sono criteri operativi di riconoscibilità, non intervalli empirici stimati.','small')

page('L’arricchimento riorganizza le taglie','03 / Competizione e costruzione della rete')
figure('02_resource_spectra','Spettri medi nella campagna primaria: stessa griglia fisica di14 intervalli, confronti N0.5/2/8. La linea verticale distingue le taglie <5.85µm.',maxheight=290)
section('Arricchimento')

page('Il collo di bottiglia è condizionale','04 / Energia, alimento e produzione')
figure('03_quantity_quality','Linee continue: full; tratteggiate: chiusura per classe. Il rapporto in basso usa produzione secondaria netta, sottraendo la respirazione basale Z. Queste curve sono del piano primario, da confrontare con il controllo fine.',maxheight=305)
section('Quantità',slice(0,2))
p('Nel piano fine, N2 è un controesempio al criterio: l’allocazione Z aumenta appena più di PC. L’esistenza di un vincolo stechiometrico non garantisce che domini la risposta dell’intera rete.','small')

page('La ritenzione N agisce anche sul riciclo','05 / Una distinzione causale')
figure('05_retention_recycling','Riduzione etaN0.67→0.335. Le linee confrontano l’effetto sullo stesso stato con il sistema riorganizzato; lisi/remineralizzazione C restano fisse.',maxheight=305)
section('Assimilazione')

page('Diversità e trasferimento restano separati','06 / La discrepanza delle tesi')
figure('04_diversity_transfer','Tutti i9 ambienti del piano comune, senza scegliere sottocoorti. I coefficienti sono associazioni di un disegno deterministico, non prove di un effetto causale della diversità.',maxheight=285)
section('La firma',slice(0,2))
p('Il filtro più vicino al campionamento di campo lascia soltanto1/9 ambienti full e4/9 own con P≥2µm e Z≥200µm osservabili. La sfida storica fallisce, ma la relazione empirica non è stata confutata: i gruppi e le misure non sono equivalenti.','small')

page('Lo stock non rappresenta la dieta','07 / Riserve, selezione e allometria')
figure('06_diet_filter','Quote nello stock C, nello stock N e nel flusso di N ingerito. La preferenza è soltanto dimensionale: la qualità della dieta è filtrata dalla rete.',maxheight=310)
section('Quantità',slice(2,None))
p('Full, capacità C piatta e fisiologia a priori generano più pattern simili in condizioni differenti. L’insieme attuale non identifica una legge universale della capacità C. L’esponente trasferito dal fit tocca un limite di stima ed è poco identificato; la presenza dei pattern non lo valida. Il volume delle cellule rimane fisso anche quando cambiano le riserve.')

page('La chiusura deve avere un significato ecologico','08 / Dinamica e risoluzione delle taglie')
figure('07_regime_robustness','PAR300/N8: le fasi delle traiettorie possono differire; si esaminano spettri e blocchi annuali. Chiusura comunitaria e per classe non regolano nello stesso modo l’insieme dei consumatori.',maxheight=310)
section('Convergenza',slice(0,2))

page('Il controllo fine modifica la spiegazione','09 / Verifica dell’intero piano')
figure('10_fine_whole_plane','Medie lunghe e5 blocchi annuali ai quattro passi. Non si richiede concordanza punto per punto in un regime oscillante. A basso N e luce alta, invece, la differenza grossolana era fra comunità quasi stazionarie.',maxheight=300)
table([['Criterio full','Piano primario h0.025','Piano fine h0.00625'],
    ['P1 / arricchimento','Sì,5/5 blocchi','Sì,5/5 blocchi'],['P2 / gruppi','Sì,5/5 blocchi','Sì,5/5 blocchi'],
    ['P3 / quantità-qualità','Sì,4/5 blocchi','No,3/5 blocchi'],['P4 / diversità-Z:P','No,0/5 blocchi','No,0/5 blocchi']],[W*.38,W*.31,W*.31])
p('La concordanza delle medie nei due passi più fini a PAR300/N0.5 non rende automaticamente la comunità resistente alla colonizzazione. La bassa luce conserva una sensibilità non risolta della composizione e della produzione.','small')

page('Le classi rare cambiano la storia','10 / Controllo conservativo dell’assemblaggio')
figure('11_rare_trait_pool','La stessa quantità di PN/ZN è redistribuita fra tratti potenziali. Gli inventari C/N restano identici; nessun nuovo nutriente entra nel box. Il controllo con sola Z distingue la disponibilità dei tratti ai due livelli.',maxheight=300)
section('Convergenza',slice(2,None))

page('Dall’equilibrio apparente all’invasibilità','11 / Una diagnosi ecologica locale')
figure('12_local_invasion','Ambienti residenti congelati. P: quote per cellula rilassate per200d prima di stimare la crescita. Z: crescita specifica di una classe rara. L’assenza funzionale non equivale a un tasso di invasione negativo.',maxheight=365)
p('Le classi P rare con crescita positiva hanno tassi0.14-0.41d-1 e residui delle quote circa1×10-15d-1. Alcune classi grandi non sono pienamente rilassate in200d: non si usano per sostenere la conclusione. La diagnosi locale da sola non dimostra una stabilità globale, ma concorda con la riorganizzazione effettiva dopo reintroduzione.')
p('Una densità estremamente bassa può rendere lungo il recupero, anche con crescita potenziale positiva. Precisione finita, storia iniziale e assenza di una sorgente di organismi possono quindi influenzare lo spettro tardivo. È necessario distinguere successione, persistenza nella finestra e coesistenza resistente all’invasione.','small')

page('L’ordine degli eventi lascia memoria','12 / Confronti da stati identici')
figure('08_sequence_memory','Differenza fra interazioni d’ordine, dopo aver sottratto gli effetti additivi dei singoli impulsi precoci/tardivi. Segni ed effetti dipendono dalla fase iniziale. I passi cambiano soltanto la perturbazione, non il burn-in.',maxheight=315)
section('Memoria',slice(0,1))
p('I_NL=(Y_NeLl-Y_Ne)-(Y_Ll-Y0); I_LN=(Y_LeNl-Y_Nl)-(Y_Le-Y0). Si confronta I_NL-I_LN. Le due sequenze hanno la stessa dose N e luminosa e l’ambiente è comune dal giorno6. La sottrazione non rende indipendenti tutti i feedback: misura una interazione nella storia del sistema.','small')

page('Una componente della memoria è nelle riserve','13 / Trapianto in ambiente comune')
figure('09_reserve_common_garden','A parità di abbondanze e ambiente, si trasferiscono soltanto le quote CR/A e NR/A delle due storie. Le risposte indicano una componente fisiologica, con effetti diversi su P e Z.',maxheight=315)
section('Memoria',slice(1,None))
p('Una maggiore divisione iniziale non garantisce una maggiore produzione del consumatore. Lo scambio di riserve modifica l’inventario iniziale: ogni giardino è verificato rispetto al proprio bilancio, evitando un falso pareggiamento di C/N.','small')

page('Che cosa conservare e che cosa sviluppare','14 / Decisione ecologica')
section('Decisione')
p('Guislain et al.2025 osservano il predominio di specie grandi con riduzione della luce per cDOM, contro l’ipotesi iniziale. Questo controesempio impedisce di estendere una pendenza luce-taglia del toy case a una legge ecologica universale. Nutrienti aggiuntivi, lipidi, mixotrofia e comportamento alimentare restano alternative nel reale; non sono processi da aggiungere automaticamente.')
table([['Affermazione','Portata sostenibile'],['C e N sono contabilizzati','Prerequisito verificato.'],['Il modello produce organizzazione','Sì, entro ambienti e storie dichiarati.'],['La formulazione è identificata dai pattern','No: capacità C e fisiologia restano alternative.'],['Le comunità tardive sono protette da invasione','Non dimostrato; i controlli a basso N mostrano riorganizzazione.'],['La validazione POM è generale','No: il verdetto rimane parziale.']],[W*.44,W*.56])

page('Fonti, controlli e riproducibilità','15 / Evidenza consultabile')
section('Verifica')
for author,title,doi in [
    ('Grimm & Railsback2012','Pattern-oriented modelling','10.1098/rstb.2011.0180'),
    ('Matthiessen et al.2025','Taglia, biomassa e arricchimento','10.1111/oik.10910'),
    ('Vascotto et al.2024','Distribuzioni e tratti nel Golfo di Trieste','10.1007/s42974-024-00215-4'),
    ('Goldenberg et al.2024','Upwelling e qualità dell’alimento','10.1038/s41598-024-54345-w'),
    ('García-Comas et al.2016','Diversità di prede/predatori e trasferimento','10.1098/rspb.2015.2129'),
    ('Guislain et al.2025','cDOM, luce e cambiamento della comunità','10.1002/ecy.70043'),
    ('Liefer et al.2019','Riserve e starvation N; dati in parte già usati nel fit','10.3389/fmicb.2019.00763'),
    ('Neun et al.2025','Spettro luminoso, nutrienti e trasferimento','10.1111/1365-2745.70161'),
    ('Matsumoto et al.2026','Modello alternativo di selettività per qualità','10.1016/j.ecolmodel.2025.111423'),
    ('Sharoni et al.2026','Modello di allocazione biochimica e composizione','10.1038/s41558-026-02598-w')]:
    p(f'**{author}.** {title}. [DOI](https://doi.org/{doi})','small')
p('Tesi consultate: De Pasquale2018, pp.28-100 per lo0D; Fianchini2021, pp.15-39 per luce, stechiometria e controllo trofico. La ricostruzione operativa è nel protocollo POM del repository.','small')
p('[Codice, protocolli e risultati](https://github.com/Marcofianchini/SPLAS-1.1/pull/1). Il dataset contiene221 stati/serie float64 esatti, metadati, registri e fonti delle firme. Gli allegati di questo PDF riportano anche il piano fine e i controesempi.','small')
p(f'Archivio dati: {archive["cases"]} casi, {archive["bytes"]/1024**2:.1f}MiB; SHA256:','small')
p(archive['sha256'],'small')


def footer(canvas,doc):
    canvas.saveState();canvas.setStrokeColor(TEAL);canvas.setLineWidth(.6);canvas.line(44,A4[1]-32,A4[0]-44,A4[1]-32)
    canvas.setFont('PV',7);canvas.setFillColor(MUTED)
    canvas.drawString(44,25,'SPLAS / Validazione POM dello 0D / supporto parziale')
    canvas.drawRightString(A4[0]-44,25,str(doc.page));canvas.restoreState()


OUT.parent.mkdir(parents=True,exist_ok=True);temporary=OUT.with_name(OUT.stem+'_layout.pdf')
doc=SimpleDocTemplate(str(temporary),pagesize=A4,rightMargin=44,leftMargin=44,topMargin=48,bottomMargin=43,
    title='SPLAS: validazione POM dello 0D - supporto parziale',author='SPLAS / ricerca ecologico-modellistica')
doc.build(story,onFirstPage=footer,onLaterPages=footer)
reader=PdfReader(temporary);writer=PdfWriter();writer.clone_document_from_reader(reader)
attachments=['summary.csv','summary.json','PATTERN_DECISIONS.json','VERIFICATION.json','data_manifest.json','invasion_diagnostic.json']
for name in attachments:writer.add_attachment(name,(DOC/name).read_bytes())
with OUT.open('wb') as f:writer.write(f)
manifest=dict(file_name=OUT.name,bytes=OUT.stat().st_size,pages=len(reader.pages),sha256=hashlib.sha256(OUT.read_bytes()).hexdigest(),
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),data_sha256=archive['sha256'],
    attachments={name:hashlib.sha256((DOC/name).read_bytes()).hexdigest() for name in attachments},figures=12)
(DOC/'report_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest,indent=2))
