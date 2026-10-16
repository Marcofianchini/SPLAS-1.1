# Esperimenti sulla stechiometria emergente — 16 ottobre 2026

Scheda scritta prima delle nuove integrazioni. È un'esplorazione POM, senza
calibrazione sui risultati né pretesa di riprodurre un sito. Il checkpoint POM e
il solutore conservativo restano immutati. Le nuove osservazioni registrano
l'intera storia delle classi e i flussi, mancanti nei vecchi output ridotti.

## Domande e tentativi di confutazione

1. **Qualità contro quantità.** A parità di apporto N, luce crescente può aumentare
   C:N e ridurre il rendimento dei consumatori. Cercare anche i controesempi:
   produzione e biomassa Z possono aumentare pur con rendimento inferiore.
   Un rendimento ridotto non dimostra da solo il collasso dei grazers.
2. **Retroazione alimentare.** Variare esclusivamente l'efficienza di assimilazione
   N del consumatore (.335, .67, 1), mantenendo quella C a .396 e C:N Z=5.5.
   L'N non assimilato va al detrito come già nel modello. Il controllo modifica
   anche il riciclo: non è un'ablazione della sola legge del minimo.
3. **Organizzazione della comunità.** Verificare se il gradiente luce × N e il
   feedback trofico riorganizzano taglie, biomassa ed evenness insieme. Contrasti
   senza Z e con preferenza esattamente uniforme phi=1; quest'ultima modifica
   anche la disponibilità incontrata, preservando la massima ingestione per Z.
4. **Assemblaggio dalle tesi.** Aggiunta dei primi 0,1,2,4,8,16 predatori della
   griglia storica, con N totale=8 e C totale=2246.375 mmol/m³ identici. La
   dinamica C:N deve essere seguita per classe, non dedotta dal solo C:N totale.
5. **Memoria di sequenza.** Impulso N prima/dopo 3 giorni di luce intensa, dose
   luminosa e N uguali, stesso ambiente finale; riferimenti N solo, luce sola e
   nessuna perturbazione. Differenze al ritorno non equivalgono a multistabilità.

## Regimi, misure e controlli

Scatole 0D illuminate uniformemente, temperatura fissa, senza sinking. Piano
minimo 1P+1Z: ESD 5 e42 µm, diluizione .05/d, PAR15–1215, DIN entrante .25–16.
Comunità: 35P+35Z, diluizione .02/d, PAR45/135/405, DIN entrante1/4/16.
La diluizione mantiene un ambiente definito e contabilizza tutti gli scambi;
non è un sistema marino localizzato. Nessuna aggiunta continua di N nel caso
chiuso di assemblaggio. C e N iniziali coincidono nei confronti senza predatori.

La fisiologia comune usa l'ipotesi già trasferita alle comunità dai fit Ot/Tw,
con perdite e riciclo dichiarati nel registro. Non si interpreta tale fit come
identificazione univoca o media tassonomica universale. Si verificheranno separatamente
alcuni contrasti anche con i parametri fisiologici a priori se risultano decisivi.

Output: tutti gli stati, riserve C/N, divisione e grazing per classe, uptake,
fissazione e respirazione, ingestione effettiva C/N, produzione dei consumatori,
potenziali di crescita limitata da C/N, budget e inventari, hash sorgenti e dati.
Medie finali su365 giorni nelle comunità e365 nel piano minimo, controllando
se persiste deriva o variabilità; niente selezione del solo istante favorevole.
Stock ZN:PN e ZC:PC, rendimento Zprod/ingC e TTE Zprod/NPP restano distinti.

Evenness su14 intervalli fisici logaritmici comuni1–200µm, frazioni pico/nano/micro,
taglia geometrica e95° percentile. Le tre celle luce×N(45,1),(405,1),(405,16)
sono ripetute a h/2 e70/140classi: confrontare la direzione e la forma dei pattern,
non soltanto il flusso assoluto. Bilanci e positività sono condizioni necessarie.
La specifica forma del minimo del consumatore impone un vincolo di resa; sono
emergenti il regime in cui si attiva, la dieta incontrata e le retroazioni.

## Limiti da mantenere nella discussione

C:N dinamico e conservativo non include P, Fe, acidi grassi, tossicità o
selezione comportamentale per qualità. Taglie di riferimento fisse; riserve
non modificano il kernel. Niente interpretazione tassonomica delle classi.
Correlazioni e decomposizioni non equivalgono a causalità; i controlli
meccanicistici servono a discriminare spiegazioni nel toy world.
