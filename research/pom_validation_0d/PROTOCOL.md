# Validazione POM dello 0D — protocollo prima dei nuovi run

Data: 16 ottobre 2026. Riferimento congelato: commit
`8b14dda650cde76a9559e0d8ae73ae7bc3beb0a1`. I 117 scenari del primo audit e i
166 esperimenti stechiometrici sono materiale esplorativo già esaminato, non dati
di conferma. Nessun parametro sarà stimato sui nuovi confronti. La registrazione
è interna e verificabile nei file; non è una preregistrazione esterna.

## Domanda ecologica

Le stesse regole semplici di acquisizione, riserve, predazione dimensionale e
riciclo riproducono insieme più proprietà del plancton, in ambienti ideali non
ancora esplorati? Quali proprietà richiedono realmente ciascuna ipotesi e quali
rimangono indeterminate? Una validazione POM può respingere una spiegazione:
non è il nome da assegnare automaticamente a una campagna conclusa.

Si riprende la progressione De Pasquale: competizione → assemblaggio della rete
→ gradienti e chiusura trofica (pp.28–95), quindi dinamiche (pp.96–100).
Da Fianchini si riprende la lettura delle retroazioni luce–risorse–qualità–controllo
trofico (pp.15–39). Il nuovo sistema risolve C e N conservativamente e non
identifica un rapporto fra stock con un rapporto di produzione.

## Bersagli e loro indipendenza

Le soglie seguenti sono criteri operativi di riconoscibilità nel toy case, non
stime ricavate dai dati empirici. L'esito va riportato anche cambiando finestra,
valuta e risoluzione; non si sommano punti in un punteggio unico.

| ID | Evidenza esterna e bersaglio | Regola operativa prima dei run | Limite di trasferimento |
|---|---|---|---|
| P1 | Arricchimento minerale: biomassa e taglia di comunità possono aumentare insieme (Matthiessen et al.2025; DP pp.50–67). | N entrante 0.5→8: PN e taglia geometrica in N +5%, ESD95 +20%, sia a PAR30 sia a PAR300; conservare almeno1% di PN nelle taglie <5.85µm al livello alto. | Testa sostituzione delle classi, non evoluzione né variazioni del volume entro classe. P e Si non sono risolti. PAR3 è un limite energetico da riportare separatamente. |
| P2 | Spettri non monotoni e compresenza di gruppi dimensionali (Vascotto et al.2024; DP pp.35–51). | Almeno2 delle6 celle PAR30/300 con due massimi della densità PN su14 intervalli logaritmici1–200µm; separazione ≥2 intervalli, prominenza ≥10% del massimo, valle ≤50% del picco minore. Estremi ammessi; controllare7/28 intervalli. | Analogo strutturale, non riproduzione della distribuzione empirica di MLD o di abbondanze individuali. La scomparsa dei modi con phi=1 discrimina il meccanismo, non è un'ulteriore osservazione indipendente. |
| P3 | Quantità di produttori e alimento utile possono disaccoppiarsi (Goldenberg et al.2024). | Nei6 contrasti luminosi adiacenti3→30 e30→300, riconoscere almeno2 contesti con N diverso: PC +5%, C:N della comunità e della dieta +5%, incremento relativo della produzione Z inferiore a quello di PC. Riportare anche produzione Z/NPP e frazione dell'ingestione limitata da N. | Capacità condizionale: l'esperimento empirico manipola upwelling/Si:N, qui luce/N. Si verifica un pattern di sistema, non il valore o il driver preciso del mesocosmo. La legge min(C,N) del consumatore è prescritta, non è una conferma empirica. |
| P4 | Sfida storica: associazioni opposte fra diversità delle prede/predatori e Z:P (García-Comas et al.2016; DP pp.75–84; F pp.25–26). | Nell'intero piano comune di9 ambienti, Spearman(E_P,logZ:P)<−0.2 e Spearman(E_Z,logZ:P)>0.2, sia in N sia in C; riportare Pearson, interazione e il piano esteso di15 ambienti per le due chiusure. Nessuna selezione a posteriori di sottocoorti. | Gli autori misurano entropia continua dei biovolumi individuali, qui equità della biomassa su intervalli fisici comuni. È una sfida alla firma delle tesi; non si dichiara una replica metrica dei dati. Correlazione non prova un effetto causale della diversità. |

P1 e P2 condividono gli spettri; P3 e P4 condividono parte degli stock e dei
flussi. Non sono quattro conferme statisticamente indipendenti. Si richiede
compatibilità congiunta, non quattro p-value. Il fallimento di P4 resta nel
verdetto anche se P1–P3 riescono. La diversità tassonomica non viene dedotta
dall'equità delle taglie. Nessuna specifica ampiezza o periodicità delle
oscillazioni è un bersaglio empirico in questa campagna.

## Alternative mantenute semplici

Tutte le formulazioni conservano C e N e condividono costi, ingressi e processi.
Il file `splas2/model.py` rimane immutato.

- `full`: fisiologia del riferimento trasferita senza nuovo fit, predazione
  dimensionale, mortalità Z proporzionale allo stock totale ZN.
- `own`: stessa fisiologia; mortalità quadratica dipendente dalla sola classe Z.
  Stesso coefficiente, senza compensare l'aumento del numero di nicchie. Cambia
  quindi anche la regolazione totale: non si pretende di isolare la sola forma.
- `flat`: capacità massima di riserva C priva della dipendenza dimensionale,
  ancorata allo stesso valore a5µm. La quota minima N conserva la sua allometria.
- `uniform`: phi=1 per ogni incontro. Cambiano anche disponibilità effettiva e
  saturazione; è una diversa rete di incontro, non una sola rimozione di scelta.
- `prior`: fisiologia a priori non calibrata, con gli stessi tassi ecologici di
  mortalità/riciclo/assimilazione del riferimento. Testa trasferibilità della
  spiegazione, non costituisce una nuova stima.
- `noZ`: controllo della competizione, inventari iniziali C e N pareggiati
  trasferendo lo stock iniziale del consumatore a DIC e DIN.

## Disegno spazialmente omogeneo

Piano principale: PAR3,30,300µmol fotoni m−2s−1 × DIN entrante0.5,2,8mmol N m−3;
35 classi P e35 Z, dominio1–200 e2–1000µm, profondità1m; diluizione0.02d−1,
luce costante, nessun export gravitazionale. Durata3650d, passo massimo0.025d,
salvataggio ogni2d con estremi esatti delle finestre. Inventario iniziale comune
per tutti: N=1.77 e C=2205.0225mmol m−3. La scelta dei livelli è nuova rispetto
alla campagna esplorativa; non è un campione statistico dell'oceano.

Le5 strutture più il controllo noZ danno54 run. Per full/own si aggiungono
6 ambienti ciascuno a PAR30 con mortalità Z0.05/0.4 e i3 livelli di N:12 run.
Contrasto assimilazione/riciclo a PAR30, N0.5/8: etaN0.335/0.67 × lisi e
remineralizzazione N0.05/0.15/0.45d−1, cinetiche C fisse:10 run nuovi,2
riferimenti già presenti. Il confronto separa ritenzione nel consumatore e
tempo di ritorno di N; non rende indipendenti tutti i feedback del sistema.

Controlli numerici scelti prima degli esiti: full e own nei3 angoli (PAR3,N0.5),
(PAR300,N0.5),(PAR300,N8), con passo0.0125 oppure70 classi P/Z a passo0.025:
12 run. Totale iniziale88. Si possono aggiungere soltanto controlli dichiarati
come successivi, senza cancellare quelli sfavorevoli.

Finestre principali: giorni1825–3650, con due metà e5 blocchi annuali; anche
l'ultimo anno. Analizzare distribuzioni, medie e occupancy delle taglie, non
chiedere sovrapposizione di fase a traiettorie oscillanti. Tolleranza operativa
di robustezza: stessa classificazione ecologica; variazioni >10% degli indicatori
aggregati o cambi di segno richiedono qualificare il risultato. Non si interpreta
un ordine di convergenza da tre soluzioni in regimi sensibili alla fase.

## Anatomia dei meccanismi

Ricostruire il circuito N in5 comparti: DIN, P, Z, detrito, DON. Distinguere
ingestione, produzione, escrezione, mortalità, egestione, lisi, remineralizzazione
e diluizione. Verificare il diagramma di flussi contro il RHS indipendente.
Valutare la riduzione di etaN anche sullo stesso stato istantaneo: l'effetto
diretto sulla crescita non va confuso con la successiva riorganizzazione del
riciclo. Leggere crescita/perdite per taglia, riserve e dieta prima delle medie.

## Memoria: nuovo esperimento correttivo

Usare stati esatti del riferimento full a PAR3, N0.5 e2; per ciascun ambiente
due fasi prescelte con una regola indipendente dall'esito della perturbazione:
minimo e massimo PN nell'ultimo anno. Tutte le storie e tutti i passi partono
dallo stesso array, senza rifare il burn-in. PAR alto30 per3d, impulso N0.5.
Le7 storie sono: controllo, N precoce, N tardivo, luce precoce, luce tardiva,
N precoce+luce tardiva, luce precoce+N tardivo. Impulsi ai giorni0 o3; dose
luminosa e dose N identiche per le due sequenze, ambiente comune dal giorno6.
Simulare fino al giorno36 ai passi0.025/0.0125/0.00625. Registrare gli impulsi
nel libro degli ingressi, conservando stati prima/dopo evento.

Per ogni indicatore Y, oltre al contrasto grezzo di sequenza si calcolano:
I_NL=(Y_NeLl−Y_Ne)−(Y_Ll−Y0), I_LN=(Y_LeNl−Y_Nl)−(Y_Le−Y0),
Delta=I_NL−I_LN. Usare stock al giorno6/9/16/36 e produzione integrata
nei giorni0–6 e6–36. Il contrasto rimuove gli effetti additivi dell'età del
singolo impulso e della posizione del singolo evento luminoso. Non prova, da
solo, che la memoria sia nelle riserve intracellulari.

Controllo common garden: al giorno6 trasferire le sole riserve per cellula
delle due sequenze sulla stessa abbondanza P, stock Z e comparti abiotici del
controllo. Classi donatrici con PN <10−10 del totale usano le riserve del
controllo; riportare il peso ricevente interessato. I giardini condividono
l'ambiente e le abbondanze, ma NON gli inventari iniziali C/N: lo scambio di
riserve è un intervento biologico esplicito, quantificato. Simulare1d e verificare
ogni bilancio rispetto al proprio inventario. Una differenza attribuibile alle
riserve è evidenza meccanicistica interna, non un pattern esterno indipendente.

## Decisione

Conservare una matrice struttura×pattern con riuscite, fallimenti e
indeterminazioni. Nessuna struttura viene promossa perché riproduce un solo
grafico. Se nessuna soddisfa i bersagli congiunti, dichiarare validazione
parziale e localizzare la spiegazione mancante. La conservazione degli elementi
è un prerequisito, non una validazione ecologica. Il passaggio a1D richiede un
quesito verticale aggiuntivo e un riferimento0D con limiti espliciti; non serve
risolvere ogni ambiguità biologica, ma non si userà lo spazio per occultarla.
