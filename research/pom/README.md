# SPLAS come laboratorio di Pattern Oriented Modelling

La finalità dello studio è individuare un insieme relativamente semplice di regole
capace di generare congiuntamente pattern ecologici riconoscibili. La previsione di
un sito oceanografico e la riproduzione di un valore assoluto di esportazione non
sono il criterio principale di successo. Questo orientamento è stato esplicitato
dall'autore il 16 ottobre 2026.

La procedura riprende la logica sperimentale delle tesi De Pasquale (2018) e
Fianchini (2021). L'etichetta POM è una ricostruzione metodologica attuale: non viene
attribuita retrospettivamente agli autori come dichiarazione presente nelle tesi.

1. [Ricostruzione delle argomentazioni e verifica delle metriche](argument-review.md).
2. [Protocollo operativo del nuovo studio](STUDY_PROTOCOL.md).
3. [Diagnosi dei risultati già disponibili](CURRENT_ASSESSMENT.md).
4. [Catalogo esplicito dei pattern](pattern_catalogue.json).
5. [Schede ordinate dei prossimi esperimenti](next_experiments.json).

## Checkpoint

Il commit di partenza è `ae55985e0f05ae05f0930c1bdfc1cca6f04052b6`, identificato
anche dal tag locale `splas-cn-0.2-pom-baseline-20261016`. [BASELINE.json](BASELINE.json)
registra gli hash di 333 file: modello, test, configurazioni, risultati e rapporto
precedente. I 117 scenari di comunità restano 117: questo aggiornamento non esegue
nuove simulazioni ecologiche e non modifica le equazioni.

La nuova diagnosi è retrospettiva ed esplorativa; le schede dei futuri esperimenti
sono proposte operative, non esperimenti già effettuati. Il vecchio rapporto
quantitativo rimane un documento storico del checkpoint. Questa cartella ne
aggiorna l'interpretazione finalistica e corregge i confronti non omogenei.

```bash
python scripts/analyze_pom.py --check-only
python scripts/analyze_pom.py
```

La prima istruzione verifica il checkpoint; la seconda riproduce la diagnosi e le
figure senza integrare il modello. Gli archivi numerici voluminosi non sono in Git;
il PDF dello studio precedente contiene il proprio archivio incorporato. Il
manifest permette di riconoscere le copie corrispondenti. Il nuovo PDF presenta
protocollo, catalogo, schede e diagnosi; il dettaglio eseguibile resta nel repository.

