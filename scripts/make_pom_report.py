"""Build the Italian thesis-grounded POM study report from the consolidated evidence."""
from pathlib import Path
import hashlib
import json
import html
import re

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
                               PageBreak, Image, KeepTogether)

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "research/pom"
OUT = ROOT / "output/pdf/SPLAS_Protocollo_POM_2026-10-16.pdf"
W = A4[0] - 92
NAVY = colors.HexColor("#163746")
TEAL = colors.HexColor("#147d92")
MUTED = colors.HexColor("#526771")
LIGHT = colors.HexColor("#edf4f5")
BRANCH_URL = "https://github.com/Marcofianchini/SPLAS-1.1/tree/codex/dynamic-cn-development/research/pom"
pdfmetrics.registerFont(TTFont("DV", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"))
pdfmetrics.registerFont(TTFont("DVB", "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"))
pdfmetrics.registerFontFamily("DV", normal="DV", bold="DVB", italic="DV", boldItalic="DVB")
STYLE = {
    "body": ParagraphStyle("body", fontName="DV", fontSize=11, leading=15.7, spaceAfter=8, textColor=NAVY),
    "small": ParagraphStyle("small", fontName="DV", fontSize=9, leading=12.5, spaceAfter=6, textColor=MUTED),
    "cell": ParagraphStyle("cell", fontName="DV", fontSize=9, leading=12.6, textColor=NAVY),
    "headcell": ParagraphStyle("headcell", fontName="DVB", fontSize=9, leading=12.6, textColor=colors.white),
    "h1": ParagraphStyle("h1", fontName="DVB", fontSize=22, leading=27, spaceAfter=15, textColor=NAVY),
    "h2": ParagraphStyle("h2", fontName="DVB", fontSize=12.2, leading=16, spaceBefore=9, spaceAfter=6, textColor=TEAL),
    "kicker": ParagraphStyle("kicker", fontName="DVB", fontSize=8.3, leading=11.4, spaceAfter=10, textColor=TEAL),
    "cover": ParagraphStyle("cover", fontName="DVB", fontSize=32, leading=38, spaceAfter=18, textColor=NAVY),
}


def read(name):
    return json.loads((DOC / name).read_text())


def clean(text):
    text = str(text).replace("→", " -> ").replace("–", "-").replace("—", "-").replace("\u2011", "-")
    # Improve prose spacing while preserving hashes, model identifiers, variables,
    # units and link targets. Never normalize the underlying numerical records.
    protected = r"(\b[0-9a-f]{40,64}\b|\b(?:P\d{2}|E\d{2}|L\d|[RDF]\d{2}|ESD95|SHA256|m[23])\b|[01]D|\S*_\S*)"
    chunks = re.split(r"(<[^>]+>|&(?:[A-Za-z][A-Za-z0-9]*|#x[0-9A-Fa-f]+|#\d+);)", text)
    for i in range(0, len(chunks), 2):
        s = re.sub(r"([A-Za-zÀ-ÿ])([01]D)\b", r"\1 \2", chunks[i])
        bits = re.split(protected, s)
        for j in range(0, len(bits), 2):
            b = re.sub(r"([A-Za-zÀ-ÿ])(\d)", r"\1 \2", bits[j])
            b = re.sub(r"(\d)([A-Za-zÀ-ÿµ])", r"\1 \2", b)
            b = re.sub(r";(?=\S)", "; ", b)
            b = re.sub(r",(?=[A-Za-zÀ-ÿ])", ", ", b)
            bits[j] = b
        chunks[i] = "".join(bits)
    return "".join(chunks)


def p(text, kind="body"):
    story.append(Paragraph(clean(text), STYLE[kind]))


def title(text, kicker):
    if story:
        story.append(PageBreak())
    p(kicker.upper(), "kicker")
    p(text, "h1")


def table(rows, widths, font="cell"):
    data = [[Paragraph(html.escape(clean(x)), STYLE["headcell" if i == 0 else font])
             for x in row] for i, row in enumerate(rows)]
    t = Table(data, colWidths=widths, repeatRows=1, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), NAVY),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 7),
        ("RIGHTPADDING", (0, 0), (-1, -1), 7),
        ("TOPPADDING", (0, 0), (-1, -1), 7),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, LIGHT]),
        ("LINEBELOW", (0, 0), (-1, 0), .5, NAVY),
        ("LINEBELOW", (0, 1), (-1, -1), .3, colors.HexColor("#d8e4e7")),
    ]))
    story.extend([t, Spacer(1, 10)])


def figure(path):
    img = Image(str(path))
    img.drawHeight = img.imageHeight * W / img.imageWidth
    img.drawWidth = W
    story.extend([img, Spacer(1, 9)])


def footer(canvas, doc):
    canvas.setFillColor(TEAL)
    canvas.rect(46, A4[1]-31, 30, 3, fill=1, stroke=0)
    canvas.setFont("DV", 8)
    canvas.setFillColor(MUTED)
    canvas.drawString(83, A4[1]-31, "SPLAS | Procedura di sviluppo e analisi orientata ai pattern")
    canvas.setStrokeColor(colors.HexColor("#d8e4e7"))
    canvas.line(46, 38, A4[0]-46, 38)
    canvas.setFont("DV", 7.5)
    canvas.drawString(46, 25, "16 ottobre 2026 | De Pasquale -> Fianchini -> nuova formulazione C/N")
    canvas.drawRightString(A4[0]-46, 25, str(doc.page))


def card(c):
    p(c["question"])
    for label, key in [("Mondo minimo", "minimal_world"), ("Disegno", "design"),
                       ("Controllo degli inventari", "inventory_control"),
                       ("Lettura del meccanismo", "interpretive_test"),
                       ("Spiegazione alternativa", "alternative_explanation"),
                       ("Lavoro necessario", "implementation_needed"), ("Decisione", "decision")]:
        p("<b>" + label + ".</b> " + html.escape(c[key]))
    p("<b>Output:</b> " + "; ".join(c["outputs"]) + ".", "small")
    p("Pattern " + ", ".join(c["patterns"]) + " | " + c["thesis"], "small")


story = []
baseline = read("BASELINE.json")
assessment = read("current_assessment.json")
catalogue = read("pattern_catalogue.json")
cards = read("next_experiments.json")["cards"]

p("QUADERNO OPERATIVO / STATO CONSOLIDATO", "kicker")
p("SPLAS\n<br/>dai meccanismi<br/>ai pattern", "cover")
p("Una procedura di ragionamento e sviluppo ricostruita dalle tesi<br/>De Pasquale (2018) e Fianchini (2021)", "h2")
p("""<b>Finalità.</b> Costruire un mondo ecologico relativamente semplice nel quale
struttura della comunità, controllo trofico, organizzazione verticale e dinamiche
temporali emergano dall'interazione di poche regole esplicite.""")
p("""Il criterio principale diventa la spiegazione congiunta di pattern riconoscibili
e delle condizioni in cui cambiano. Un valore di esportazione oceanica o la
previsione di un sito non costituiscono il bersaglio principale dello studio.""")
table([["Consolidato", "Prodotto di questo aggiornamento"],
       ["Codice e risultati", "Checkpoint identificato da333 hash;117 scenari di comunità conservati. Nessuna nuova simulazione ecologica; equazioni invariate."],
       ["Lettura delle tesi", "Ricostruzione di domande, contrasti, inferenze e limiti con pagine verificabili; distinzione fra testo e proposta attuale."],
       ["Protocollo POM", "18 voci centrali di diversa natura, una verifica fisiologica ausiliaria e8 schede sperimentali ordinate."],
       ["Diagnosi attuale", "Nucleo competitivo riconoscibile; relazioni diversità-trasferimento e robustezza della diversità sono domande prioritarie."]], [W*.27, W*.73])
p("""Il termine POM è una ricostruzione metodologica attuale, coerente con la finalità
esplicitata dall'autore. Non è attribuito retrospettivamente alle tesi come
protocollo formalmente dichiarato.""", "small")

title("La procedura di De Pasquale", "01 / dal caso minimo alla funzione")
p("""La tesi tratta la costruzione del modello come parte del contributo scientifico.
Il passaggio al livello successivo nasce da una domanda che il livello precedente
non può affrontare: prima capire i vincitori, poi la comunità, infine lo spazio e i flussi.""")
table([["Passaggio", "Ragionamento documentato", "DP, pagine stampate"],
       ["Regole ed equilibrio", "Allometrie e preferenza fissate; quota e requisiti di risorsa rendono interpretabili crescita e perdite.", "14-33"],
       ["Riduzione", "NPZ chiuso e Monod per seguire la competizione; prima senza predatori, poi aggiungendoli in ordine.", "34-49"],
       ["Rete e chiusura", "Spettri, limitazioni e diversità insieme; due chiusure possono avere Z:P simile e comunità molto diverse.", "50-64"],
       ["Pattern congiunti", "Apporto N, mortalità Z, diversità-biomassa e legami diversità-trasferimento; spiegazione degli scostamenti.", "65-79"],
       ["Comportamento e robustezza", "Selettività, perdite, risoluzione, assimilazione e Kz; distinguere taglie esatte e risposta qualitativa.", "80-95"],
       ["Tempo e spazio", "Forzante periodica, poi colonna, detrito, luce, shading e mixing; osservare che cosa cambia nella spiegazione.", "96-140"],
       ["Funzione e sviluppo", "Export come prima esplorazione della funzione; una risposta troppo rigida ai parametri è anche un limite da indagare.", "141-153"]], [W*.23, W*.62, W*.15])
p("""<b>Trasferimento al nuovo studio.</b> Conservare l'ordine dei quesiti, i controlli
ridotti e la lettura attraverso spettri e flussi. Il modello a C:N dinamico non
deve conservare ogni taglia vincente o coefficiente numerico della tesi.""")
p("""Nella numerazione De Pasquale la pagina PDF è maggiore di uno rispetto alla
pagina stampata. La tesi distingue spesso analisi di comportamento, confronto
con pattern e prospettiva applicativa futura: queste funzioni rimangono distinte.""", "small")

title("La procedura di Fianchini", "02 / seguire le ramificazioni di una forzante")
p("""Fianchini parte da un sistema noto, costruisce due archetipi e introduce la
regolazione della luce. Una variazione semplice viene seguita nelle sue
conseguenze su nutrienti, stechiometria, taglie, grazing, bloom e trasferimento.""")
table([["Passaggio", "Uso nel ragionamento", "F, pagine"],
       ["Finalità e riferimento", "Esplorazione di sistemi generati da leggi semplici; archetipi aperto/chiuso come metro di confronto ideale.", "5-6,15-22"],
       ["Luce costante", "Variazioni0.7-1.3 volte il riferimento; da profili e indicatori a ipotesi di effetti non lineari.", "23-32"],
       ["Luce periodica", "Due massimi luminosi annuali dichiarati; successioni, fase e riorganizzazioni tardive.", "33-39"],
       ["Bilancio e profondità", "Relazioni fra luce, proxy ottica del mixing e crescita; confronto con dominio meno profondo.", "35-44"],
       ["Limite produttivo", "La taglia del detrito e il suo legame con la stechiometria motivano lo sviluppo successivo.", "45-47"]], [W*.24, W*.61, W*.15])
p("""<b>Tre correzioni al confronto precedente.</b> I due massimi annuali non sono un
errore accidentale: sono la scelta di p.33. L'efficienza di pp.39-40 normalizza
l'export pesato per biomassa rispetto all'irradianza e non coincide con export/NPP.
Il corpo della tesi descrive risposte miste della diversità nel sistema aperto
(p.32), più articolate della sintesi generale di p.46.""")
p("""La lettura fedele conserva anche queste tensioni fra sintesi e risultati locali.
Il quesito diventa individuare regimi e catene di effetti; una sola pendenza non
può rappresentare tutta la tesi.""")
p("""Le pagine stampate e PDF di Fianchini coincidono. Formule e figure decisive alle
pp.26,32,33,36,40,46 sono state controllate visivamente.""", "small")

title("Il ciclo operativo del nuovo studio", "03 / una domanda, un contrasto, una spiegazione")
for n, head, text in [
    (1, "Definire il pattern", "Identificare fonte, scala, misura e regime; distinguere osservazione, teoria, risultato di un modello e ipotesi di sviluppo."),
    (2, "Costruire il caso minimo", "Esplicitare allometrie, kernel, vincoli C/N e input; scrivere quali proprietà sono prescritte e quali devono emergere."),
    (3, "Leggere il riferimento", "Seguire spettri e serie, poi crescita, grazing, riserve, riciclo e uscite; nelle colonne includere profili e flussi alle facce."),
    (4, "Discriminare una spiegazione", "Modificare un processo o confrontare due strutture nominate; conservare inventari e forzanti quando richiesto dal quesito."),
    (5, "Valutare la firma congiunta", "Confrontare coesistenza, forma, direzione, posizione relativa e fase dei pattern; verificare i contrasti a risoluzione e soglie alternative."),
    (6, "Decidere il passo successivo", "Mantenere la spiegazione più semplice sostenuta dall'insieme dei pattern; conservare controesempi e strutture indistinguibili, aggiungendo solo complessità motivata.")]:
    p(str(n) + ". " + head, "h2")
    p(text)
p("""Il POM descritto da Grimm e Railsback 2012 usa più pattern per guidare struttura,
selezione dei sottomodelli e parametri. Qui la logica viene applicata a
compartimenti e classi dimensionali; la scelta di uno scopo esplicativo è propria
di questo progetto.""", "small")
p("""<b>Regola di decisione:</b> nessun punteggio totale automatico. Non compensare un
meccanismo decisivo mancante contando molte metriche correlate. Il paesaggio dei
regimi fa parte del risultato.""")

short_status = {
    "P01": "Oracolo ridotto e vincitore riconoscibile", "P02": "Sostegno qualitativo",
    "P03": "Sensibilità presente; ablazione esatta da fare", "P04": "Condizionale; controesempio in altro regime",
    "P05": "Sostegno nella coda dimensionale", "P06": "Massimo interno con rimbalzi",
    "P07": "Firma congiunta non riprodotta nella coorte", "P08": "Risposta presente, non monotona globalmente",
    "P09": "Serie totali disponibili; successioni da salvare", "P10": "DBM possibile; controllo uniforme mancante",
    "P11": "Diagnosi parziale", "P12": "Multi-massimo non dimostrato nei profili medi",
    "P13": "Meccanismo da discriminare", "P14": "Feedback esplicito; composizione rilevante",
    "P15": "Controlli disponibili; transizioni da caratterizzare", "P16": "Bilancio parzialmente operazionalizzato",
    "P17": "Disaccoppiamento presente; limiti numerici", "P18": "Sensibilità presente; ipotesi da selezionare",
    "P19": "Verifica ausiliaria parziale"}
for first, last, heading, kicker in [(0, 9, "Catalogo: competizione e comunità", "04 / voci P01-P09"),
                                    (9, 19, "Catalogo: spazio, tempo e funzione", "05 / voci P10-P19")]:
    title(heading, kicker)
    p("""Le voci comprendono pattern ecologici, oracoli e ipotesi: non sono conferme
osservative indipendenti. Lo stato è una diagnosi retrospettiva; i criteri dei
prossimi esperimenti sono registrati separatamente.""", "small")
    rows = [["ID", "Pattern / domanda", "Stato corrente"]]
    for item in catalogue["patterns"][first:last]:
        rows.append([item["id"], item["pattern"], short_status[item["id"]]])
    table(rows, [W*.08, W*.59, W*.33])
    p("""Il catalogo completo associa a ogni voce natura dell'evidenza, pagine,
osservabili, controfattuale, limite e criterio del nuovo studio. P19 è ausiliario:
le colture vincolano la plausibilità senza decidere da sole il successo ecosistemico.""", "small")

title("Che cosa resta riconoscibile nel 0D", "06 / riuso dei risultati, nessun nuovo run")
figure(ROOT / "results/pom/01_assembly_resource_resolution.png")
p("""<b>A.</b> Senza predatori domina la classe1µm; introducendo predatori si aprono
spazi a taglie meno consumate. Il confronto Monod/riserve cambia anche le cinetiche:
mostra continuità qualitativa, non l'effetto isolato della stechiometria.""")
p("""<b>B.</b> La coda ESD95 cambia in modo diverso sotto due condizioni di C.
Occorre identificare il regime di risorsa e il ruolo dei consumatori; non imporre
una relazione N-taglia monotona sull'intero paesaggio di condizioni.""")
p("""<b>C.</b> La dominante resta circa20µm fra35,70,140 classi, ma la distribuzione
su14intervalli fisici comuni cambia. Robustezza della dominante e robustezza della
diversità sono conclusioni diverse.""")

title("Diversità e trasferimento come domanda", "07 / la firma di pattern va spiegata insieme")
figure(ROOT / "results/pom/02_diversity_trophic_patterns.png")
p("""<b>A-B.</b> Nella scansione di N compaiono massimi interni di evenness, insieme
a rimbalzi: il pattern non richiede una parabola imposta. Le linee collegano scenari
ordinati per apporto e non rappresentano repliche statistiche.""")
p("""<b>C-D.</b> La coorte di24casi0D usa famiglie di perturbazione analoghe a DP
pp.75-79. Le associazioni cambiano con la famiglia: la firma congiunta negativa
per le prede e positiva per i predatori non è ancora riprodotta in questa coorte.""")

title("Una diagnosi più fedele delle metriche", "08 / valuta, coorte e finestra fanno parte della domanda")
p("""Il confronto storico usa biomassa in N. Con C:N dinamico occorre calcolare
ZN:PN dal carbonio dei consumatori e dalla loro stechiometria; ZC:PC e il rapporto
di produzione restano osservabili aggiuntivi. Si usano rapporti di stock medi
dell'ultima finestra, non la media dei rapporti istantanei del decennio.""")
table([["Coorte / misura", "Prede", "Predatori"],
       ["24 casi: Pearson, ZN:PN", "-0.229", "-0.138"],
       ["24 casi: Spearman, ZN:PN", "+0.155", "+0.382"],
       ["24 casi: Pearson, ZC:PC", "-0.314", "-0.232"],
       ["24 casi: Pearson, produzioneZ_C/NPP_C", "+0.019", "-0.035"],
       ["5 casi mortalitàZ: Pearson, ZN:PN", "+0.973", "+0.986"],
       ["9 casi apportoN: Pearson, ZN:PN", "-0.637", "-0.525"]], [W*.64, W*.18, W*.18])
p("""La coorte comprende D00,D01,D02,D04,D05,D06: forzanti costanti e chiusura
di comunità. Non coincide con l'identico insieme storico. Non vengono calcolati
p-value o inferenze sulla popolazione degli ecosistemi: sono associazioni fra
esperimenti deterministici.""")
p("""Il cambio di segno fra misure e famiglie segnala una domanda meccanicistica
prioritaria: quali combinazioni di selettività, chiusura, risorse e qualità delle
prede generano il legame diversità-funzione? Scegliere retrospettivamente soltanto
una coorte o trasformazione favorevole eliminerebbe proprio questa informazione.""")
p("""Le quattro colture restano una sfida fisiologica ausiliaria. Il risultato
peggiore della persistenza per Micromonas limita la trasferibilità dei parametri;
non impone automaticamente tassonomia o maggiore complessità al toy case.""")

title("La nuova stechiometria nel ragionamento", "09 / leggere quote, composizione e flussi")
p("""Nel sistema chiuso, aumentando la luce da0.7 a1.3 volte il riferimento, il C:N
aggregato aumenta di4.847. Il termine entro classe è-2.798, quello di composizione
+7.645: la riorganizzazione della comunità sovrasta la variazione delle quote.
La decomposizione è un'identità contabile e non una causa identificata.""")
p("""La domanda nuova è seguire luce/N, riserve, dieta ingerita, crescita del
consumatore e riciclo. Un C:N comunitario non misura automaticamente il C:N della
dieta e non descrive la risposta di ogni cellula. Il vincolo omeostatico dei
consumatori e la composizione delle riserve sono ipotesi dichiarate.""")
table([["Stati correnti", "Lettura POM"],
       ["F00: massima interna62.5m; F02:52.5m", "Esistenza di regimi con DBM; isolare il gradiente che li sostiene."],
       ["F05 senza shading:132.5m", "Confronto con F02 a forzante abbinata: la retroazione ottica riorganizza il profilo."],
       ["R00: massima nella prima cella2.5m", "Il campo chiamato DBM_depth indica una posizione del massimo; questo caso ha una massima superficiale."],
       ["Profili medi ispezionati: un solo picco interno o nessuno", "Massimi multipli non dimostrati con prominenza2,5,10%; richiedono anche output istantanei."],
       ["R00: N trattenuto83.51mmol/m2/d", "Regime di arricchimento ammesso; da distinguere dall'equilibrio e da una transizione puramente endogena."],
       ["Raddoppio lisi+remin: export-89.6%, NPP+0.073%", "Disaccoppiamento nel toy case; separare i due stadi e i tempi prima di attribuire un meccanismo unico."]], [W*.47, W*.53])
p("""Una parametrizzazione astratta di C può essere un controllo sperimentale utile.
Non serve trasformarla in un oceano realistico per questa domanda; servono
inventari, flussi e condizioni di risorsa espliciti.""")

title("E01: ricominciare dai vincitori", "10 / primo blocco del nuovo studio")
card(cards[0])
p("""Le simulazioni esistenti costituiscono il punto di confronto. Il nuovo
protocollo abbinato evita che l'aggiunta del predatore sia confusa con
l'aggiunta di nutrienti. La scelta riprende la logica di DP, rendendo esplicito un
controllo in più reso utile dalla nuova valuta C/N.""", "small")

title("E02: interrogare le regole della rete", "11 / primo blocco del nuovo studio")
card(cards[1])
p("""Le regole imposte sono gli oggetti dell'esperimento. Un aumento della
complessità non è il criterio di avanzamento: avanzare significa sapere quali
regole sono sufficienti a sostenere l'insieme dei fenomeni scelti.""", "small")

title("E03: il paesaggio dei pattern0D", "12 / primo blocco del nuovo studio")
card(cards[2])
p("""Il criterio discriminante è la firma congiunta: espansione dello spettro,
coesistenza, diversità e rapporto trofico nello stesso insieme dichiarato di
condizioni. Gli esiti contrari delimitano il regime e suggeriscono la prova
successiva, senza essere rimossi dalla diagnosi.""", "small")

title("Quando aggiungere tempo e spazio", "13 / sviluppo successivo ordinato")
table([["Scheda", "Domanda e contrasto", "Condizione per avanzare"],
       ["E04 / robustezza", "Ripetere i contrasti0D decisivi con35/70/140 classi e passo dimezzato; diversità su intervalli comuni.", "Identificare quali pattern sono stabili alla scala della domanda."],
       ["E05 / successioni", "N periodico e costante a media uguale; due fasi; registrare spettri e flussi nel tempo.", "Distinguere ciclo, transiente lungo e memoria."],
       ["E06 / colonna", "PAR uniforme, attenuazione, poi shading e mixing separati; stessi inventari e input areali.", "Spiegare comparsa/scomparsa di DBM e differenziazione delle taglie."],
       ["E07 / luce", "Intensità0.7-1.3; cicli annuali/semestrali a dose uguale; mixing distinto da Zeu.", "Seguire le ramificazioni da riserve e dieta fino a bloom e organizzazione."],
       ["E08 / detrito", "Pool unico/classi; mapping, lisi, remineralizzazione e sinking separati; flussi a piani fissi.", "Aggiungere capacità esplicativa della ripartizione riciclo-export."]], [W*.2, W*.5, W*.3])
p("""<b>Output da aggiungere.</b> Snapshot dimensionali e verticali durante
transiente e cicli tardivi, crescita e perdite per classe, ingestione C/N e dieta,
flussi di riciclo e particelle per origine, inventari e flussi a piani fissati.
Le sole serie totali e gli spettri medi finali non ricostruiscono la successione.""")
p("""<b>Confronto dei tempi di riciclo.</b> Un solo stadio a0.15/d ha un tempo medio
di6.67 giorni. Due stadi consecutivi a0.15/d ciascuno hanno media13.33 giorni.
Per confrontare la forma del ritardo a media abbinata si può usare0.30/d per
ciascuno dei due stadi, dichiarando che la distribuzione dei tempi resta diversa.""")
p("""Tutte le schede di questa sezione sono protocolli proposti. Le famiglie già
disponibili ne coprono parti, non costituiscono esecuzione automatica del nuovo studio.""", "small")

title("Controlli necessari alla spiegazione", "14 / precisione pertinente allo scopo")
table([["Classi P/Z", "Dominante P", "Evenness nativa", "Evenness su14bin comuni"],
       ["35", "19.31µm", "0.151", "0.311"], ["70", "19.98µm", "0.155", "0.323"],
       ["140", "20.31µm", "0.299", "0.388"]], [W*.2, W*.24, W*.27, W*.29])
p("""La dominante è stabile, la distribuzione e l'equiripartizione meno. Il
rebinning conserva la massa e usa gli stessi intervalli log-ESD, ma non elimina
tutta la differenza. Una sola percentuale di variazione della biomassa non
verifica questo pattern.""")
p("""Nella colonna i confronti40/80/160celle non dimostrano ancora la convergenza
complessiva. Per il POM occorre verificare anche presenza dei massimi, nicchie,
segni e fase dei contrasti. Le verifiche annuali da stato comune non dimostrano
che l'attrattore decennale sia identico a risoluzione diversa.""")
p("""<b>Coerenza.</b> Unità, positività, conservazione di C e N e trasferimenti
consistenti rimangono obbligatori. Una conversione errata non è una semplificazione
del sistema. Il checkpoint soddisfa i controlli di inventario registrati.""")
p("""<b>Parsimonia.</b> Contare stati, regole, parametri liberi e ipotesi aggiunte;
confrontarli con i pattern spiegati e con la leggibilità dei meccanismi. Una
griglia più fitta non è una nuova ipotesi biologica. Nuovi nutrienti o gruppi
entrano se diventano necessari a una domanda definita.""")
p("""<b>Esiti ammessi.</b> Sostegno descrittivo; meccanismo discriminato; compatibilità
condizionale; controesempio; confronto non omogeneo; non ancora valutato. La
revisione deve conservare questi esiti e non sostituirli con un unico voto.""")

title("Provenienza e materiali riutilizzabili", "15 / stato consolidato del progetto")
p("<b>Checkpoint:</b> " + baseline["baseline_commit"], "small")
p("<b>SHA256 del modello:</b> " + baseline["model_sha256"], "small")
p("""Manifest del checkpoint:333 file verificati;117 firme modello/caso e archivi di
comunità controllati. Nessuna nuova integrazione ecologica in questo aggiornamento.""")
p('''Documentazione completa ed eseguibile: <link href="'''+BRANCH_URL+'''" color="#147d92">cartella research/pom nel repository SPLAS</link>.
Contiene ricostruzione delle argomentazioni, protocollo, catalogo, schede,
manifest e diagnosi numerica riproducibile.''')
p('''<link href="https://github.com/Marcofianchini/SPLAS-1.1/pull/1" color="#147d92">Pull request di sviluppo SPLAS-CN</link>.
Il rapporto quantitativo precedente è conservato come documento del checkpoint;
questa analisi ne aggiorna la lettura rispetto alla finalità POM.''')
p("Ambiente della nuova analisi: " + "; ".join(k+" "+v for k,v in assessment["analysis_environment"].items()) + ". I run congelati conservano il proprio ambiente di riferimento.", "small")
p("Fonti principali", "h2")
sources = [
    ("De Pasquale2018, tesi di dottorato, copia165pagine. Pagina PDF=stampata+1.", "https://arts.units.it/bitstream/11368/2929831/2/M_DePasquale_TESI_FINAL.pdf"),
    ("Fianchini2021, tesi magistrale, copia corretta51pagine fornita. Hash e nome della copia in SOURCE_PROVENANCE.json.", None),
    ("Grimm e Railsback 2012. Pattern-oriented modelling: a multi-scope for predictive systems ecology. Abstract verificato; riferimento di metodo.", "https://doi.org/10.1098/rstb.2011.0180"),
    ("Banas 2011. Adding complex trophic interactions to a size-spectral plankton model. Abstract dell'editore verificato; pattern e transitori.", "https://doi.org/10.1016/j.ecolmodel.2011.05.018"),
    ("Ward etal2013. Top-down and bottom-up controls in a global size-structured plankton food-web model. Fonte del quesito storico nelle tesi; metadati MIT verificati.", "https://doi.org/10.1093/plankt/fbt097"),
    ("Irigoien etal2004. Global biodiversity patterns of marine phytoplankton and zooplankton. Bersaglio storico di DP; metadati verificati, dati non rianalizzati.", "https://doi.org/10.1038/nature02593"),
    ("García-Comas etal2016. Prey size diversity hinders biomass trophic transfer and predator size diversity promotes it. Abstract e figura1 verificati.", "https://doi.org/10.1098/rspb.2015.2129"),
    ("Liefer etal2019. The macromolecular basis of phytoplankton C:N:P under nitrogen starvation. Vincolo fisiologico ausiliario del checkpoint.", "https://doi.org/10.3389/fmicb.2019.00763")]
for text, url in sources:
    if url:
        p('<link href="'+url+'" color="#147d92">'+html.escape(text)+'</link>', "small")
    else:
        p(html.escape(text), "small")

OUT.parent.mkdir(parents=True, exist_ok=True)
doc = SimpleDocTemplate(str(OUT), pagesize=A4, rightMargin=46, leftMargin=46,
                        topMargin=54, bottomMargin=52,
                        title="SPLAS: protocollo POM e ricostruzione delle due tesi",
                        author="Codex - analisi del progetto SPLAS",
                        subject="Consolidamento del checkpoint e procedura di sviluppo orientata ai pattern")
doc.build(story, onFirstPage=footer, onLaterPages=footer)
manifest = {
    "date": "2026-10-16", "pdf": str(OUT.relative_to(ROOT)), "pages": doc.page,
    "bytes": OUT.stat().st_size, "sha256": hashlib.sha256(OUT.read_bytes()).hexdigest(),
    "baseline_commit": baseline["baseline_commit"], "model_sha256": baseline["model_sha256"],
    "new_ecological_simulations": 0,
    "inputs_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                     for path in sorted(DOC.glob("*"))
                     if path.is_file() and path.suffix in {".json", ".md"} and path.name != "report_manifest.json"},
    "authoring_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
}
(DOC / "report_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False)+"\n")
print(json.dumps({"pdf": str(OUT), "bytes": manifest["bytes"],
                  "pages": manifest["pages"], "sha256": manifest["sha256"]}))
