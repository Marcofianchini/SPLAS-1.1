# Raffinamento dopo il controllo primario

Il primo controllo mostra che full a PAR300/N0.5 ha stock quasi stazionari
diversi a h0.025 e0.0125: C:N14.41 contro21.05 e taglia1.77 contro2.52µm.
Questa differenza non è spiegata da una diversa fase delle oscillazioni.
Il fallimento del criterio operativo di robustezza viene conservato.

Si aggiungono prima dell'integrazione di questi controlli:

1. Tutto il piano full di9 ambienti a h0.00625, con le stesse condizioni iniziali,
   durata3650d e finestre. I criteri P1–P4 non cambiano. È un controllo successivo
   dell'intero piano, senza selezionare le sole condizioni favorevoli.
2. I3 angoli originari a h0.003125 per distinguere sensibilità del transiente
   da persistenza della differenza alle scale più fini.
3. A PAR300/N0.5, continuazioni di730d a h0.00625 da due stati esatti al giorno3650:
   il riferimento h0.025 e il suo controllo h0.0125. Ciascuna continuazione ha
   un controllo identico e una reintroduzione rara delle classi potenziali.
   La reintroduzione redistribuisce10−5 del PN residente uniformemente fra le
   classi, con quote di riferimento, prelevando prima la stessa frazione da
   tutti i pool P. C e N totali restano identici: l'eventuale differenza C delle
   nuove riserve è compensata in DIC. Lo stock Z, DIN, DON e detrito rimangono
   identici; la variazione di DIC è registrata. È un intervento biologico raro,
   non un apporto di nutriente né una modifica delle regole di crescita.

Le4 continuazioni non dimostrano da sole multistabilità: entrambi gli stati
potrebbero essere ancora transitori e le classi rare possono richiedere più
tempo per crescere. Si confrontano direzione del cambiamento, crescita delle
classi reintrodotte, stock e spettro. L'assenza di risposta non prova stabilità
globale. Non si riparte da un inventario aggiustato tacitamente.

Totale aggiunto:16 run; gli esiti primari e i loro fallimenti non vengono
sostituiti. Il modello e l'integratore restano congelati. Questa verifica è
interna e successiva, non una nuova conferma empirica indipendente.
