# Protocollo operativo per il nuovo studio SPLAS

## Domanda e criterio di successo

Quale insieme relativamente semplice di allometrie, competizione per risorse,
preferenze di predazione e riciclo è sufficiente a generare **congiuntamente**
struttura dimensionale, diversità, controllo trofico, organizzazione verticale e
dinamiche temporali riconoscibili negli ecosistemi planctonici? Come cambia questa
capacità quando il C:N diventa una proprietà dinamica che partecipa ai flussi?

L'obiettivo è esplicativo e orientato ai pattern. Il toy case non deve riprodurre un
sito, un calendario geografico o un valore assoluto di export. Numeri, unità e
bilanci servono a definire regimi, osservabili e relazioni; permettono di
riconoscere se due configurazioni stanno davvero mettendo alla prova la stessa
ipotesi. Una configurazione ideale può essere deliberatamente estrema e utile.

Il protocollo riprende la progressione documentata nelle tesi. La formulazione
POM è attuale e proposta: non si attribuisce agli autori un'adozione formale non
dichiarata. Grimm e Railsback 2012 descrivono il POM come progettazione, selezione e
calibrazione mediante più pattern a scale e livelli diversi. Qui trasferiamo
questa logica a un modello di compartimenti e classi dimensionali; non occorre
trasformare SPLAS in un modello a individui.

## Regole di lavoro

1. Formulare la domanda prima di cambiare un processo. Scrivere che cosa è imposto,
   che cosa dovrebbe emergere e quale esito renderebbe la spiegazione insufficiente.
2. Ripartire dal caso minimo che conserva la domanda. Luce uniforme non significa
   assenza di energia per la fotosintesi; assenza di un gradiente e assenza di un
   processo sono contrasti diversi.
3. Tenere un riferimento leggibile per ogni livello di complessità. Un unico
   riferimento valido per tutti i quesiti non è necessario.
4. Seguire prima la dinamica interna e i flussi; poi gli indicatori aggregati.
   Una curva simile con trasferimenti diversi costituisce una spiegazione alternativa.
5. Usare più pattern e più scale. Evitare di contare come conferme indipendenti
   quantità algebricamente legate o misurate sugli stessi scenari.
6. Confrontare ipotesi alternative. I confronti fra formulazioni devono dichiarare
   anche le cinetiche, i costi e i trasferimenti che cambiano insieme.
7. Conservare i controesempi. Un pattern presente solo in un regime può essere un
   buon risultato se il regime e il suo meccanismo sono espliciti.
8. Aggiungere complessità soltanto per risolvere una domanda o una discrepanza
   definita. La letteratura aggiorna ipotesi e vincoli, non impone l'inclusione di
   tutti i processi conosciuti.

## La scala di costruzione, tratta dalle tesi

| Livello | Domanda | Configurazione e contrasto | Risultato da rendere leggibile | Riferimento |
|---|---|---|---|---|
| L0 - regole e valuta | Quali tratti e flussi sono prescritti? | Allometrie, kernel, C e N, definizione di crescita e perdite | Unità compatibili, scambi elementari, assunti e costi liberi | DP pp.18-33; F pp.6-13 |
| L1 - competizione | Chi vince senza predatori? | NP chiuso; controllo analitico ridotto e modello dinamico | Requisito di risorsa, esclusione, riserve | DP pp.28-39 |
| L2 - costruzione della rete | Come compaiono nicchie, picchi e lacune? | Aggiunta ordinata di predatori; kernel selettivo e ablazione uniforme | Bilanci di crescita/perdita delle classi che emergono | DP pp.35-49 |
| L3 - paesaggio 0D | Le stesse regole generano più pattern ambientali? | Apporto N, mortalità Z, selettività, perdite P; chiusura alternativa | Coesistenza, espansione della taglia, diversità-biomassa, Z:P | DP pp.50-95 |
| L4 - tempo | Che cosa aggiunge la dinamica alla media? | N periodico, ampiezza e frequenza; riferimento a media uguale | Successioni, fase, ritardi, transiente e memoria | DP pp.96-100 |
| L5 - spazio | Quali gradienti organizzano la colonna? | 1D con trasporto e detrito; PAR uniforme/attenuata, poi shading e mixing | DBM possibili, taglie per profondità, nutrienti e flussi | DP pp.104-140 |
| L6 - luce come nuova domanda | Come si ramificano gli effetti di una forzante semplice? | Archetipi aperto/chiuso, luce0,7-1,3 volte il riferimento, poi periodica | Stechiometria, qualità delle prede, controllo trofico, bloom e transizioni | F pp.15-39 |
| L7 - bilancio integrato | Quale bilancio nello strato mescolato sostiene la fioritura? | Mixing prescritto nello stesso dominio; luce e grazers separati | Zeu, crescita biologica integrata, trasporto e stock distinti | F pp.35-44 |
| L8 - funzione | Quale riorganizzazione raggiunge il comparto detritico e l'export? | Pool unico/classi particellari; ipotesi di riciclo e sinking | Ripartizione dei flussi e loro dipendenza dalla comunità | DP pp.141-152; F pp.39-46 |

La sequenza è un ordine di ragionamento. La nuova implementazione contiene già
molti di questi processi, ma la spiegazione deve poterli ricostruire per livelli.
Non occorre riscrivere né riattivare tutto insieme. Quando disattivare un termine
romperebbe il bilancio, occorre definire un controllo conservativo esplicito.

## Un ciclo operativo per ogni esperimento

**Scheda iniziale.** Registrare domanda; pattern interessati; fonti e loro natura
(osservazione, teoria, risultato di un altro modello, ipotesi della tesi); regime;
formulazione; riferimento; modifica; inventari iniziali; forzanti; osservabili;
alternative meccanicistiche; criteri di esito; output necessari. Le schede in
`next_experiments.json` applicano questa struttura.

**Lettura del riferimento.** Esaminare serie di N e C, spettri P e Z, crescita e
perdite per taglia, ingestione, riserve, riciclo e uscite. Nelle colonne aggiungere
profili di luce, nutrienti e flussi alle facce. Usare mappe taglia-tempo e
taglia-profondità prima di comprimere il risultato in una media.

**Contrasto.** Modificare un processo alla volta quando si vuole identificare il suo
effetto. Aperto/chiuso rimane un confronto fra archetipi ecologici: separatamente
si costruiscono confronti con sola variazione di riciclo, input, trasporto o
esportazione. Quando più modifiche costituiscono proprio l'ipotesi da confrontare,
chiamare il confronto strutturale e riportarle tutte.

**Spiegazione.** Individuare quale classe guadagna una crescita netta positiva,
quale risorsa o perdita cambia e dove passa la materia. Una decomposizione o una
correlazione propone una lettura; per discriminarla si usa il contrasto successivo.

**Verifica del pattern.** Valutare direzione, forma, coesistenza, organizzazione
spaziale e temporale, dipendenza dal regime. Registrare esito e controesempi;
verificare il contrasto a passo o griglia alternativi quando pertinente. Il
criterio principale non è una percentuale di errore sul valore di export.

**Decisione.** Mantenere il modello più semplice che spiega l'insieme di pattern
richiesti nel regime studiato. Se due strutture sono indistinguibili, dichiararlo e
individuare un nuovo pattern discriminante. Se un meccanismo aggiuntivo modifica
soltanto un numero secondario, la sua inclusione resta facoltativa.

## Contratto delle misure

| Oggetto | Misura primaria | Distinzione da conservare |
|---|---|---|
| Struttura delle taglie | Biomassa per intervallo log-ESD; frazioni pico/nano/micro, modi, quantili e medie | Una media non identifica l'espansione dello spettro; ESD95 non è la massima classe con tracce numeriche |
| Diversità | Evenness in N per il confronto storico; in C e su intervalli fisici comuni come diagnostiche | Il numero di classi potenziali cambia il denominatore; equità delle taglie non è ricchezza tassonomica |
| TTE e accoppiamento | ZN:PN come proxy storico; ZC:PC e produzioneZ_C/NPP_C separati | C:N dinamico rende i due rapporti di stock differenti; già DP pp.78-79 esamina rapporti di produzione |
| Luce e C:N | C:N per classe, comunità e dieta ingerita; decomposizione entro classe/composizione | Una comunità più ricca di C non implica che tutte le sue cellule lo siano |
| Fioritura | Crescita biologica netta e accumulo dello stock, fase e ampiezza | Trasporto e accumulo di nutrienti possono contribuire allo stock; numero di massimi luminosi è prescritto |
| Massimi verticali | Profili N, C e clorofilla separati; massimi interni con prominenza | Il massimo superficiale non è un DBM; massimi medi e istantanei rispondono a domande diverse |
| Profondità | Dominio, mixing prescritto, Zeu e profondità di bilancio integrato | Non usare Zeu come mixing fisico; non definire una soglia dalla crescita e poi usarla come unica prova della crescita |
| Export | Flussi gravitazionali C e N a piani fissati, export/NPP e perdite trofiche superiori | F pp.39-40 normalizza export pesato per biomassa rispetto all'irradianza: è un'altra metrica, da ricostruire nel codice storico |
| Tempo | Stessa finestra nei confronti, più fasi del transiente e cicli tardivi | Le tesi spesso mediano tutti i10anni; gli spettri attuali mediano gli ultimi364giorni |

I rapporti attuali sono rapporti degli stock medi sulla finestra. Una media dei
rapporti istantanei può differire: anche quest'operazione va abbinata quando si
ricostruisce esattamente una figura storica. Le analisi correnti non pretendono
di ricostruire quelle curve punto per punto.

Le formule di limitazione delle tesi sono diagnostiche del loro modello a quota.
Nel nuovo modello l'uptake, la divisione e l'accumulo di C sono distinti: non si
sommano automaticamente percentuali di limitazione come se fossero frazioni
indipendenti. Vanno salvati tassi e flussi effettivi e, quando utile, risposte a
controfattuali istantanei dichiarati.

## Confronto di formulazioni e parsimonia

Il riferimento Monod a C:N fisso è utile per l'oracolo R* e per rendere leggibile
l'assemblaggio. Non costituisce un confronto che isola la sola stechiometria,
perché cambiano anche crescita e uptake. Un confronto più specifico richiede una
formulazione di controllo con gli stessi processi cinetici, costi e trasferimenti,
e un vincolo stechiometrico dichiarato. Se mantenere tutte queste proprietà non è
possibile, la differenza va attribuita alle due strutture complessive.

Si propone una scala di ipotesi: conversione C/N implicita storica; controllo
ridotto fisso; riserve C/N con consumatore omeostatico; eventuali varianti di
qualità alimentare o detrito. La conversione implicita rimane riferimento storico,
non una variante a bilanci equivalenti. Per ogni ipotesi riportare stati,
regole, parametri liberi, prior e ipotesi nuove. Un numero elevato di classi non è
di per sé una maggiore complessità biologica; più parametri calibrabili e più
retroazioni indipendenti lo sono.

La selezione confronta i compromessi: pattern spiegati e loro robustezza,
leggibilità dei meccanismi, costo di stati e parametri aggiuntivi. Non assegniamo
oggi pesi retrospettivi per far vincere la nuova formulazione. I dati di coltura
restringono la plausibilità fisiologica, ma un errore su una specie non obbliga a
creare un gruppo tassonomico quando il quesito resta di organizzazione ecosistemica.

## Che cosa resta necessario nel toy case

**Coerenza.** Unità, positività e bilanci elementari devono essere corretti. Una
conversione errata non è una semplificazione ecologica. Inventari e input devono
essere confrontabili quando il quesito lo richiede. La correzione del coefficiente
C-volume è distinta dall'introduzione di un nuovo meccanismo biologico.

**Robustezza dei pattern.** Controllare sopravvivenza, classi dominanti, bande,
massimi, segni e fase rispetto a passo, risoluzione dimensionale, griglia verticale
e soglie di rilevanza. La precisione del flusso assoluto può essere secondaria;
una griglia che cambia il numero dei massimi o il regime della comunità coinvolge
direttamente l'inferenza POM. Gli attuali limiti del40/80/160 non invalidano per
estensione gli esperimenti0D.

**Regime.** Un accumulo persistente di DIN è ammesso in un toy case se fa parte
della domanda. Va distinto da uno stato stazionario o da un regime periodico
chiuso sul ciclo. Per attribuire una transizione alla rete, confrontare anche la
spiegazione alternativa dell'accumulo secolare di risorsa. Il serbatoio di C può
essere un controllo ideale dichiarato; non occorre trasformarlo in chimica dei
carbonati per svolgere questa funzione sperimentale.

**Vincolo empirico pertinente.** I pattern devono avere un riferimento ecologico,
specificando scala, misura e regime. Non è richiesto un hindcast di sito. Si
mantengono separati oracoli matematici, risultati di modelli, osservazioni e ipotesi
di sviluppo: il catalogo non è un elenco di conferme indipendenti.

## Output necessari nei prossimi esperimenti

Il checkpoint attuale conserva serie totali, spettri e profili medi dell'ultima
finestra, e stati iniziale/finale. Non permette di ricostruire retroattivamente
l'intera storia di successione dimensionale e tutti i massimi istantanei.

Per i nuovi run salvare: configurazione e hash; stato iniziale; inventari e
forzanti; snapshot dimensionali e verticali durante il transiente e i cicli
tardivi; flussi C/N per processo; termini di crescita e perdita per classe;
ingestione e dieta ponderata; particelle e flussi a profondità fissate. Produrre
medie sia dell'intero intervallo sia di finestre dichiarate, senza presentarle come
equivalenti. La cadenza va verificata rispetto ai processi da riconoscere, non
imposta uguale per ogni protocollo.

## Ordine del prossimo studio

Il primo blocco operativo è **E01-E03**, descritto nelle schede: competizione e
assemblaggio, chiusura/selettività, quindi paesaggio0D con i pattern congiunti.
Questo blocco decide quali regole e quali regimi costituiscano il riferimento
spiegabile della nuova formulazione. E04 verifica i suoi pattern a risoluzione
alternativa. Solo dopo si ricostruiscono tempo, colonna, luce/mixing e detrito
(E05-E08). La corrispondenza a un valore assoluto di export non decide l'ordine.

Gli esiti vengono classificati come sostegno descrittivo, meccanismo discriminato,
compatibilità condizionale, controesempio, confronto non omogeneo o non ancora
valutato. Nessun conteggio riassuntivo sostituisce la spiegazione della firma
congiunta. Il modello viene consolidato a ogni livello e ampliato quando il
livello successivo pone una domanda che quello corrente non può affrontare.

## Riferimenti di metodo

- Grimm e Railsback 2012, DOI https://doi.org/10.1098/rstb.2011.0180: definizione e uso multi-criterio dei pattern; abstract verificato in PubMed/PMC e indice istituzionale UFZ.
- Banas 2011, DOI https://doi.org/10.1016/j.ecolmodel.2011.05.018: modello spettrale semplice, pattern di diversità e transitori lunghi; abstract dell'editore verificato.
- De Pasquale2018: fonte principale della progressione; copia e hash in SOURCE_PROVENANCE.json.
- Fianchini2021: fonte principale dell'analisi ramificata della luce; copia corretta e hash in SOURCE_PROVENANCE.json.
- García-Comas etal2016, DOI https://doi.org/10.1098/rspb.2015.2129: associazioni di diversità e trasferimento con rapporto di biomassa come proxy; abstract e figura1 verificati in PubMed.

