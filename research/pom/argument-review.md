# Dalle tesi alla procedura di ragionamento di SPLAS

## Materiali, localizzatori e metodo

Sono state esaminate le copie identificate in `research/SOURCE_PROVENANCE.json`:
De Pasquale, *The role of plankton size in community structure, biodiversity and
biogeochemical fluxes: a modeling approach* (2018), e Fianchini, *Il ruolo della
luce nella struttura e biodiversità di comunità planctoniche e sui flussi
bio-geo-chimici: un approccio modellistico* (2021, copia corretta fornita).

Qui **DP p.** indica la pagina stampata di De Pasquale: la pagina PDF è maggiore di
uno. **F p.** indica la pagina stampata e PDF di Fianchini, che coincidono. Sono
stati letti i passaggi metodologici e interpretativi nelle introduzioni, nella
costruzione del modello, negli esperimenti, nella discussione e nelle conclusioni.
Sono state controllate visivamente le pp.26,32,33,36,40,46 di Fianchini per formule,
assi e qualificazioni che l'estrazione testuale non conserva completamente.

Le tesi non dichiarano un protocollo formale denominato Pattern Oriented
Modelling. La compatibilità con il POM è una mia interpretazione, coerente con la
finalità esplicitata dall'autore il 16 ottobre 2026. La ricostruzione del testo e
il protocollo proposto rimangono distinguibili. Nessun risultato delle nuove
simulazioni viene usato per riscrivere retrospettivamente le intenzioni degli autori.

## Ricostruzione dall'alto: domanda, contributo e portata

**De Pasquale.** La domanda è come le dimensioni degli organismi, attraverso tratti
allometrici e selettività della predazione, organizzino la comunità e influenzino
la funzione biogeochimica. Il contributo dichiarato comprende anche la costruzione
progressiva del modello, non soltanto il suo risultato finale (DP pp.14-17).
Il percorso 0D permette di rendere leggibili competizione e controllo trofico;
il passaggio 1D verifica come luce e trasporto cambino quella spiegazione. La
possibilità di future applicazioni realistiche è presente, ma le stime di export
sono esplicitamente preliminari (DP pp.16,102,141-153).

**Fianchini.** La domanda è come la modulazione della luce riorganizzi un sistema
di partenza noto. La semplicità relativa consente di esplorare interazioni e
risultati inattesi (F pp.5-6). Il caso di riferimento aperto/chiuso costruisce un
metro di confronto deliberatamente ideale (F pp.15-17). La procedura aggiunge
prima una variazione assoluta della luce, poi una dinamica periodica, e segue
ramificazioni verso stechiometria, predazione, distribuzione delle taglie,
fioriture e trasferimento del carbonio (F pp.23-47).

La finalità attuale riprende questo nucleo esplicativo. Le ambizioni applicative
future presenti nei testi non diventano requisiti obbligatori per il nuovo toy case.

## Ricostruzione dal basso: evidenze e passaggi inferenziali

| Affermazione e posizione | Supporto esaminato | Passaggio inferenziale e assunzioni | Limiti o contraddizioni | Trasferimento proposto |
|---|---|---|---|---|
| Poche regole di taglia possono organizzare una comunità, DP pp.14-17 | Due assunti generali, allometrie, preferenza preda-predatore; confronto 0D e 1D | Le proprietà collettive non sono prescritte come spettri di biomassa | Le allometrie e il dominio di taglia sono imposti; non tutto ciò che appare è indipendente dagli assunti | Dichiarare regole imposte e pattern emergenti separatamente |
| La spiegazione richiede prima un sistema più semplice, DP pp.34-35 | NPZ chiuso, crescita Monod, riciclo immediato, nessuna forzante esterna | Rimuovere quota e ritardi permette di seguire i trasferimenti di N | Quel sistema non ha le stesse cinetiche del modello a quota | Usarlo come controllo didattico e matematico; progettare a parte il confronto stechiometrico controllato |
| Senza predatori vince la classe con il minimo requisito, DP pp.28-39 | Soluzioni di equilibrio, caso NP, introduzione del primo predatore | La riduzione della risorsa esclude competitori più esigenti; la predazione apre spazio a taglie meno predate | Il risultato R* scalare richiede il regime a singola risorsa e le sue ipotesi | Ripetere NP, NPZ con un predatore, poi assemblaggio; nel C/N dinamico diagnosticare anche riserve e DIC |
| Picchi e lacune derivano dalla rete di preferenze, DP pp.40-49 | Introduzione progressiva di predatori, spiegazione delle classi vincenti e dei vuoti | La pressione distribuita sulle prede crea coesistenza e competizione indiretta | Gli estremi del dominio troncano alcune preferenze, DP pp.45-46 | Registrare nicchie, flussi e sensibilità ai limiti di taglia; non richiedere identiche classi vincenti |
| La chiusura dello zooplancton è una scelta strutturale decisiva, DP pp.58-64 | Confronto Zj per Ztot contro Zj al quadrato; spettri e diversità molto diversi a Z:P simile | Più indicatori distinguono comunità che una sola metrica aggregata non distingue | Riprodurre soltanto Z:P non seleziona una struttura della rete | Confrontare chiusure sull'intera firma di pattern e sul bilancio dei flussi |
| Più apporto di N consente la coesistenza di taglie grandi e piccole, DP pp.65-68 | Spettri, frazioni pico/nano/micro, limitazioni e biomasse sotto diversi apporti | L'arricchimento riduce il controllo di risorsa e sostiene predatori che limitano le piccole cellule | Non è soltanto una relazione tra N e taglia media; il ruolo del regime è essenziale | Valutare coda dimensionale e persistenza delle piccole classi insieme alle medie |
| Diversità-biomassa con massimo intermedio, DP pp.69-70 | Scansione di N e insieme di perturbazioni ambientali | Diverse organizzazioni della comunità generano una forma collettiva riconoscibile | Equità fra taglie non equivale a diversità tassonomica; distribuzione delle perturbazioni condiziona il grafico | Richiedere una regione con massimo interno, specificare valuta, coorte e scala; evitare una parabola imposta |
| Diversità di prede e predatori hanno associazioni opposte con Z:P, DP pp.75-79 | Confronto con García-Comas2016 e analisi di casi dispersi attorno alla relazione | Il modello propone una possibile connessione tra struttura e trasferimento | L'autore precisa che pendenza e causalità non sono stabilite; esamina già anche i rapporti di produzione | Conservare ZN:PN per il confronto storico e aggiungere ZC:PC e produzioneZ/NPP come osservabili distinti |
| La robustezza importante riguarda anche il pattern, DP pp.80-95 | Selettività, perdite, assimilazione, Kz e risoluzione di taglia | Taglie vincenti e ampiezze possono cambiare mentre la risposta qualitativa persiste | Una sola configurazione risolta più finemente non dimostra robustezza di tutte le risposte | Raffinare i contrasti e ricostruire gli stessi pattern su intervalli fisici comuni |
| La forzante periodica modifica la successione, DP pp.96-100 | Apporto sinusoidale, ampiezza e frequenza variate, spettri nel tempo | La risposta dinamica aggiunge informazione rispetto alle medie | Il testo attribuisce a sin²(pi*t/365) anche un periodo biennale: il periodo effettivo è365giorni | Conservare l'esperimento e dichiarare il periodo algebrico; confrontare medie di input uguali |
| Luce e risorsa rompono la simmetria verticale, DP pp.104-119 | Prima colonna quasi omogenea nelle taglie; poi gradiente luminoso e DBM | La luce aggiunge una differenziazione dell'ambiente lungo z | Un DBM è un controllo di comportamento e una proprietà possibile, non un bersaglio metrico universale | Confrontare luce uniforme e attenuata mantenendo geometria e apporto areale |
| Ombreggiamento e mixing modificano i profili, DP pp.123-140 | Variazione di attenuazione, self-shading, profili di diffusività e loro gradienti | Feedback di biomassa e sovrapposizione delle regioni attive spiegano cambiamenti di struttura | La seconda massima è discussa come ipotesi; alcune prove variano due coefficienti insieme | Ablazioni separate, conteggio di massimi con soglie dichiarate, diffusione conservativa completa |
| L'export serve anche a interrogare i limiti del modello, DP pp.141-152 | Un pool OM, rapporto Redfield, sensibilità piccola a diversi parametri | Una risposta troppo rigida può segnalare che il collegamento tra rete e funzione è debole | La vicinanza a una media globale non dimostra correttezza meccanicistica | Studiare la risposta qualitativa a produzione, dimensioni e riciclo; trattare i valori assoluti come diagnostici |
| Il riferimento è un mondo controllato, F pp.15-22 | Due configurazioni, profili, spettri, diversità, C:N e serie temporali | Conoscere il riferimento precede l'introduzione della nuova forzante | Aperto e chiuso differiscono in più processi; sono regimi compositi | Preservare entrambi come archetipi e aggiungere contrasti a un solo processo |
| La luce rivela una rete di effetti, F pp.24-32 | Luce +/-10,20,30%; TTE proxy, taglie, mappe di limitazione, DBM, bloom | Una variazione semplice dell'input può innescare risposte non lineari | Il corpo descrive diversità mista e C:N non strettamente monotono; le conclusioni pp.45-46 sono più generali | Catalogare pattern condizionali e collegamenti fra metriche, senza trasformare la sintesi in legge universale |
| La stagionalità è un esperimento dichiarato con due massimi, F pp.33-39 | Formula16, figura4.13, analisi nel tempo, transizioni tardive a luce alta | Fase, successione e cambi di organizzazione possono emergere sotto una forzante nota | I due massimi di luce sono prescritti; una transizione tardiva non dimostra bistabilità | Separare risposta forzata, transiente lungo e cambiamento persistente; eventuale isteresi richiede nuove prove |
| La profondità critica è esplorata con una proxy ottica, F pp.36-44 | MLD*=Zeu, r da variazioni di stock, relazioni con fioritura, dominio100m | Il testo stesso indica MLD* come tracciante; si cercano relazioni fra produzione, luce e perdite | Profondità ottica e mixing sono costrutti diversi; disuguaglianze e notazioni non sono uniformi fra pp.36,40,43 | Mantenere il quesito; prescrivere mixing e calcolare separatamente Zeu e bilancio biologico integrato |
| L'efficienza di export diminuisce con la luce, F pp.39-40,43,46 | Export pesato per biomassa e poi normalizzato per irradianza, figura4.22 | La normalizzazione separa in parte crescita della biomassa e rendimento rispetto all'input | Non è l'e-ratio export/NPP; il dettaglio della pesatura richiede codice/archivi originali | Non dichiarare fallimento o conferma mediante il segno dell'e-ratio; mantenere due diagnostiche nominate distintamente |
| La taglia del detrito è un'estensione motivata, F p.46 | Limite della conversione fissa e legame proposto taglia-riciclo-export | Una lacuna specifica suggerisce il prossimo meccanismo da esplicitare | La relazione detrito-taglia è una nuova ipotesi, non una conseguenza già verificata | Confrontare pool unico e classi di particelle con trasferimenti elementari conservati |

## Discrepanze fra formulazioni globali e risultati locali

**Finalità.** In DP pp.14-17 convivono sviluppo teorico e prospettiva applicativa;
in pp.111-119 e150-152 è esplicito il valore di un comportamento ragionevole nel
mondo ideale. In F pp.5-6,17 e45-47 la priorità dell'esplorazione qualitativa è
ancora più netta. La richiesta attuale seleziona questa finalità come centrale,
senza negare che nei testi siano presenti prospettive future di previsione.

**Luce, taglia e diversità.** F p.27 descrive un maggiore intervallo di taglie e una
distribuzione più equa in alcune simulazioni; p.32 mostra un effetto misto nel
sistema aperto e differente nel chiuso. P.46 sintetizza invece una diminuzione di
diversità. La procedura trasferibile è interrogare il regime e la composizione,
non esigere una pendenza negativa in ogni esperimento. Lo stesso vale per il C:N:
la figura4.4 a p.26 non è una retta monotona, mentre la conclusione p.45 parla di
aumento lineare generale.

**Forzante ed emergenza.** I due picchi annuali di F p.33 sono deliberati. In DP
p.96 la forma sinusoidale al quadrato produce un periodo diverso da quello
attribuito nel testo alla sinusoide di base. Sono due casi distinti: si conserva
l'intento sperimentale, si dichiara la funzione effettiva e si valuta che cosa
aggiunge la dinamica della comunità rispetto alla forzante.

**Proxy e processi.** DP pp.78-79 già confronta rapporti di produzione con Z:P:
non sarebbe corretto attribuire al lavoro storico l'ignoranza di questa distinzione.
F p.36 riconosce il carattere di tracciante della proxy ottica del mixing. Il nuovo
studio può rendere i processi più espliciti senza presentare questa esplicitazione
come una confutazione automatica del valore esplorativo dei lavori precedenti.

## Revisioni del giudizio precedente

| Prima | Dopo | Motivo e conseguenza | Verifica |
|---|---|---|---|
| Il mancato confronto di campo è il limite principale del progetto | Il successo principale è la spiegazione congiunta di pattern attraverso poche regole | Finalità dell'autore e lettura DP/F; la previsione di sito rimane un possibile progetto diverso | Documentata nei testi e nella richiesta attuale |
| Un e-ratio crescente contraddice l'efficienza decrescente di Fianchini | Le due quantità non sono confrontabili direttamente | Denominatori differenti; ricostruire la misura storica prima del confronto | F pp.39-40, figura4.22 controllata visivamente |
| Diversità crescente con luce è opposta alla tesi nel suo complesso | Il testo documenta effetti misti e dipendenti dal regime | Riprendere corpo e figure, qualificando la generalizzazione delle conclusioni | F pp.27,32,46 |
| Due massimi annuali segnalano una stagionalità accidentalmente duplicata | Sono una scelta esplicita del toy case | Il pattern da spiegare è la risposta ecologica alla forzante | F p.33, formula16 e figura4.13 controllate |
| Stabilità numerica e corrispondenza del valore di export giudicano l'intero modello | Coerenza e robustezza dei pattern sono requisiti distinti; valori di export sono secondari | Anche POM richiede che selezione e coesistenza non siano create dalla griglia | DP pp.89-95 e nuova diagnosi del checkpoint |

## Questioni aperte utili al nuovo studio

Restano da distinguere sostegno descrittivo e spiegazione del meccanismo per:
relazioni diversità-trasferimento, qualità alimentare in C/N, molteplicità dei
massimi verticali, successioni dimensionali e transizioni tardive. Gli archivi
ridotti attuali non conservano l'intera storia taglia-profondità-tempo dei dieci
anni. Serviranno nuovi output per queste domande, non soltanto più simulazioni.
Il catalogo e le schede del nuovo studio specificano le osservazioni interne
necessarie. Eventuali run originali delle tesi permetterebbero confronti storici
più esatti, ma non sono un prerequisito per il programma POM qui definito.

