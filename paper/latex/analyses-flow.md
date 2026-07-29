# Flusso delle analisi

Questo file riassume cosa abbiamo fatto, in quale ordine e perché. È una guida
informale.

## La domanda

Vogliamo capire se l'uso delle scarpe chiodate con piastra in carbonio è
associato alla prevalenza di infortuni e se l'associazione cambia tra:

1. periodo di preparazione;
2. periodo competitivo.

È uno studio trasversale ma che raccoglie dati retrospettivi (possibile recall bias). Possiamo descrivere associazioni, non dimostrare che il carbonio causi gli infortuni o che l'esposizione li preceda.

## 1. Costruzione del campione

- 355 atleti avevano consenso, eleggibilità e questionario completo.
- Un atleta è stato escluso perché mancava il periodo dell'infortunio.
- Due atleti sono stati esclusi perché avevano un periodo con frequenza totale
  pari a zero: in quel caso la carbon share avrebbe denominatore zero.
- Campione finale: 352 atleti.
- Ogni atleta contribuisce con due righe, preparazione e competizione:
  704 osservazioni atleta-periodo.

Le due righe dello stesso atleta non sono indipendenti. Per questo le analisi
principali sono multilivello, con un'intercetta casuale condivisa per atleta.

## 2. Outcome ed esposizioni

L'outcome è binario e specifico per periodo:

- almeno un infortunio nel periodo di preparazione: sì/no;
- almeno un infortunio nel periodo competitivo: sì/no.

Per cinque domini di allenamento gli atleti hanno riportato separatamente la
frequenza settimanale con carbonio e senza carbonio. Le risposte vanno da 0 a
5, dove 5 significa “5 o più”.

Da queste risposte ricaviamo tre quantità.

### Adoption

`Adoption = 1` se nel periodo è riportata almeno una frequenza con carbonio;
altrimenti `0`.

### Carbon share

```text
carbon share =
  somma delle frequenze con carbonio
  -------------------------------------------------
  somma con carbonio + somma senza carbonio
```

Esempio: score carbonio 4 e score non-carbonio 6 danno share 4/(4+6) = 40%.

La share misura la quota relativa dell'allenamento riportato svolta con
carbonio, non il numero assoluto di sedute.

### Absolute carbon-frequency score

È il solo numeratore: la somma delle cinque frequenze con carbonio.

Non è un conteggio esatto di sedute uniche, perché:

- le risposte “5 o più” sono codificate come 5;
- una seduta potrebbe contribuire a più di un dominio.

Nel paper la chiamiamo quindi *absolute carbon-frequency score*, non “numero
esatto di sedute”.

## 3. Covariate e aggiustamento

Tutti i modelli aggiustano per lo stesso insieme di covariate:

- età, sesso competitivo, altezza, peso;
- anni di esperienza;
- infortunio nelle due stagioni precedenti;
- frequenza totale, specifica per periodo, nei cinque domini;
- frequenza di palestra;
- disciplina.

Le frequenze totali includono carbonio e non-carbonio. Quindi gli effetti di
share o score carbonio sono interpretati a parità delle frequenze totali
misurate e delle altre covariate.

Le covariate continue sono standardizzate. Le discipline sono parzialmente
pooled, così le categorie piccole non producono stime completamente separate e
instabili.

## 4. Modello multilivello primario

Il modello primario è una regressione logistica Bayesiana multilivello:

- intercetta diversa per preparazione e competizione;
- intercetta casuale per atleta, condivisa fra i due periodi;
- effetto dell'adoption diverso per periodo;
- effetto della carbon share diverso per periodo, solo tra gli utilizzatori;
- stesso insieme di covariate di aggiustamento.

Tra gli utilizzatori la share è modellata con un **processo Gaussiano
approssimato nello spazio di Hilbert (HSGP, kernel Matérn 5/2)**. Un modello
lineare imporrebbe che ogni aumento, per esempio di 10 punti percentuali di
share, produca lo stesso cambiamento nei log-odds di infortunio in ogni parte
della scala. Non avevamo una ragione clinica forte per imporre questa forma:
l'associazione potrebbe aumentare, attenuarsi o raggiungere un plateau.

La prima versione di questo modello (vedi cronologia git) usava una *natural
spline* a 3 gradi di libertà per la stessa ragione: lasciare che fossero i
dati a descrivere la curvatura, senza imporre linearità. Il limite di quella
scelta era che i 3 gradi di libertà erano comunque fissati da noi a priori, non
stimati; con 2 o 4 gradi di libertà la forma e la grandezza del contrasto
sarebbero potute cambiare, e non avevamo quantificato quella fragilità. L'HSGP
risolve questo punto direttamente: la quantità di curvatura che la curva può
assumere è governata da due iperparametri stimati dal modello, separatamente
per preparazione e competizione--un'ampiezza (`gp_alpha`, quanto la curva può
allontanarsi da una retta) e una lunghezza di scala (`gp_rho`, quanto
rapidamente può cambiare pendenza). Entrambi hanno prior debolmente
informative: `gp_alpha ~ Half-Normal(0, 0.7)`, la stessa scala già usata per
gli altri coefficienti di esposizione; `gp_rho ~ InvGamma(2,94, 0,80)`,
calibrata perché il suo 1° e 99° percentile corrispondano a un decimo e al
doppio dell'intervallo osservato della share centrata tra gli utilizzatori--
abbastanza ampia da non escludere a priori curve piuttosto mosse, ma senza
concentrare probabilità su lunghezze di scala troppo corte perché i dati le
possano risolvere. La curva è rappresentata con 40 funzioni di base; abbiamo
verificato che questo numero fosse sufficiente controllando che la densità
spettrale alla frequenza più alta inclusa fosse trascurabile rispetto a quella
alla frequenza più bassa, per l'intera distribuzione a posteriori di
`gp_alpha` e `gp_rho`--se non lo fosse stata, avremmo dovuto aumentare il
numero di funzioni di base.

Il risultato più importante del passaggio da spline a HSGP è che la direzione
e l'ordine di grandezza dell'associazione competitiva restano gli stessi, ma
l'incertezza è più onesta. Il contrasto P75−P25 in competizione passa da 14,3
punti percentuali (CrI 95% 2,7--25,8) con la spline a 3 gradi di libertà a
13,2 punti (CrI 95% −0,2--29,1) con l'HSGP: la stima centrale è quasi
identica, ma l'intervallo di credibilità ora include di poco lo zero, perché
il modello non sta più assumendo una quantità fissa di flessibilità--la sta
stimando, con l'incertezza che questo comporta. Confrontato con la
specificazione lineare (6,0 punti, CrI −0,3--12,6) e con prior HSGP più ampie
(14,6 punti, CrI −0,1--30,7), tutte e tre le specificazioni concordano sulla
direzione positiva e collocano almeno il 96,9% della probabilità a posteriori
su un'associazione positiva, pur toccando tutte lo zero al limite
dell'intervallo. Trattiamo questo come un segnale reale e clinicamente
rilevante--la probabilità direzionale resta alta in tutte e tre le
specificazioni--la cui grandezza precisa i dati attuali non permettono ancora
di fissare con certezza.

La curva HSGP non crea gruppi e non identifica automaticamente soglie
cliniche. Come già valeva per i nodi della spline, i valori 25%, 50%, e 75%
restano strumenti di lettura della curva continua, non punti in cui
affermiamo che il rischio cambi improvvisamente. Di conseguenza non
interpretiamo il 50% come change-point, e non presentiamo 13,2 punti
percentuali come un effetto stabile e preciso.

## 5. I confronti primari: percentili empirici

Per riassumere la curva continua abbiamo usato P25, P50 e P75 della share
osservata tra gli utilizzatori, separatamente per periodo.

| Periodo | P25 | P50 | P75 |
|---|---:|---:|---:|
| Preparazione | 28,6% | 42,1% | 57,1% |
| Competizione | 36,0% | 50,0% | 71,4% |

Questi sono percentili della distribuzione osservata. Non sono valori fissi del
25%, 50% e 75%.

I due contrasti primari, al fine dell'interpretazione clinica, sono:

1. **Adoption:** utilizzatore a P50 contro non-utilizzatore.
2. **Dose:** utilizzatore a P75 contro utilizzatore a P25.

Abbiamo usato P25 e P75 perché riassumono una variazione ampia ma ancora
centrale e ben supportata dai dati, evitando di basare il risultato sulle code
estreme. Non sono soglie biologiche o prescrittive.

Il contrasto P75-P25 population-standardized è:

- preparazione: 0,0 punti percentuali, CrI 95% da -12,0 a 10,2;
- competizione: +13,2 punti, CrI 95% da -0,2 a 29,1.

## 6. Stime di "rischio" per un profilo fisso (mediano)

Abbiamo due modi diversi, complementari, di trasformare il modello in rischi
assoluti.

### Stime population-standardized

Manteniamo la distribuzione osservata di covariate e discipline e calcoliamo una
media su tutto il campione. L'effetto casuale dell'atleta viene integrato.

Queste stime rispondono, in modo semplificato, alla domanda:

> quale prevalenza media prevedrebbe il modello nella popolazione osservata
> sotto ciascuno scenario di esposizione?

Sono le stime usate per i contrasti primari P75-P25 e P50-non-user.

### Curva per un profilo fisso (mediano)

Per rendere il risultato più leggibile clinicamente abbiamo anche fissato tutte
le caratteristiche a un unico profilo di riferimento e posto l'effetto casuale
dell'atleta a zero.

Questa curva risponde alla domanda:

> per lo stesso profilo, come cambia il rischio previsto al cambiare della
> share?

Non è una media di popolazione e non è un calcolatore clinico validato per tutti. Il
modello può produrre una curva analoga per qualsiasi profilo specificato, ma
servirebbe validazione esterna prima di usarlo per decisioni individuali.

## 7. Perché nella curva compaiono 25%, 50% e 75%

Nella Figura 7, 25%, 50% e 75% sono valori esatti di carbon share:

- 25%: un quarto dello score totale è carbonio;
- 50%: metà è carbonio;
- 75%: tre quarti sono carbonio.

Sono stati scelti dopo il fit come tre punti di lettura semplici:

- sono valori semplici da immaginare;
- sono equidistanti;
- 50% è la mediana pooled tra gli utilizzatori e il centro della curva HSGP;
- 25% e 75% sono simmetrici attorno al 50%;
- tutti e tre sono nel supporto osservato di entrambi i periodi.

Non sono percentili, gruppi, cut-off, nodi o soglie cliniche. La
curva è continua e viene stimata una volta sola. I passaggi 25→50 e 50→75 sono
semplicemente letture adiacenti della stessa curva, non due effetti stimati da
modelli diversi.

Per il profilo fisso, le mediane posteriori sono:

| Periodo | Share 25% | Share 50% | Share 75% |
|---|---:|---:|---:|
| Preparazione | 45,0% | 47,5% | 45,0% |
| Competizione | 46,5% | 54,3% | 64,3% |

In competizione la mediana della curva aumenta di 7,8 punti fra 25% e 50% e di
10,0 punti fra 50% e 75%. Questi due incrementi descrivono la forma della mediana
della curva; non sono due contrasti indipendenti e 50% non è un change-point.

## 8. Controlli del modello primario

Abbiamo verificato:

- convergenza delle catene e assenza di divergenze;
- posterior predictive checks;
- PSIS-LOO raggruppato per atleta;
- confronto con un modello identico ma senza termini di esposizione;
- sensibilità a una share lineare;
- sensibilità a prior più ampie;
- poststratificazione per età e sesso;
- stime descrittive condizionate per sesso ed età.

La direzione dell'associazione competitiva resta positiva, ma la sua grandezza
dipende dalla forma funzionale e dalle prior. Per questo non trattiamo +13,2
punti come una costante precisa o come una soglia. La poststratificazione e le
stime per sottogruppo elencate sopra sono state ricalcolate con la curva HSGP
(modello `weighted_model_hsgp.stan`, stessa base HSGP del modello primario):
il contrasto competitivo pesato è +13,7 punti (CrI 95% da -0,2 a 30,0),
sostanzialmente invariato rispetto al +13,2 non pesato, e resta coerente
(13,1-14,1 punti, $P_{\mathrm{dir}}$ 96,8% in ogni strato) condizionando su
sesso o fascia d'età--nessuna evidenza descrittiva che l'associazione sia
concentrata in un singolo sottogruppo.

## 9. Perché abbiamo aggiunto i due modelli share vs absolute

L'analisi primaria mostrava un'associazione con la share in competizione. La
domanda successiva era:

> è davvero la quota relativa di carbonio a essere associata all'infortunio,
> oppure basta la quantità assoluta riportata di uso del carbonio?

Abbiamo quindi costruito due modelli multilivello lineari appaiati:

1. **modello relative-use:** esposizione = carbon share;
2. **modello absolute-use:** esposizione = absolute carbon-frequency score.

I due modelli hanno:

- stesso campione e stesso outcome;
- stessa adoption;
- stessi effetti per periodo;
- stesse covariate e frequenze totali;
- stessa struttura multilivello e stesse prior;
- esposizione lineare standardizzata tra gli utilizzatori.

Li abbiamo stimati separatamente. Share, score carbonio e frequenza totale sono
matematicamente collegati; inserirli simultaneamente renderebbe i coefficienti
difficili da identificare e da spiegare.

La standardizzazione crea una scala confrontabile:

- +1 SD di share = +24,1 punti percentuali, qui da 37,3% a 61,4%;
- +1 SD di score assoluto = +3,36 unità dello score carbonio, qui da 2,83
  a 6,19.

Per ciascun modello confrontiamo lo stesso profilo da -0,5 SD a +0,5 SD: la
differenza totale è quindi +1 SD. Le 3,36 unità sono unità dello score
ordinale costruito dal questionario, non 3,36 sedute uniche.

Questa è un'analisi secondaria aggiunta per chiarire l'interpretazione del
modello primario; non va presentata retroattivamente come ipotesi
prospetticamente pre-specificata.

## 10. Risultato del confronto share vs absolute

Differenza di rischio prevista per +24,1 punti percentuali di share o +3,36
unità dello score assoluto tra gli utilizzatori:

| Periodo | Carbon share | Score assoluto carbonio |
|---|---:|---:|
| Preparazione | +3,6 pp (-3,9; 11,2) | +1,5 pp (-6,9; 9,4) |
| Competizione | +7,9 pp (1,4; 14,4) | +8,3 pp (1,6; 15,9) |

In preparazione entrambi gli intervalli includono zero.

In competizione entrambi i modelli mostrano un'associazione positiva e le stime
sono quasi uguali. Quindi i dati non sostengono la frase:

> “è la share, non la frequenza assoluta”.

La conclusione corretta è:

> tra gli utilizzatori, una maggiore esposizione riportata al carbonio nel
> periodo competitivo è associata a maggiore prevalenza di infortunio, sia
> quando l'esposizione è espressa come share relativa sia quando è espressa
> come score di frequenza assoluta.

I dati non permettono di stabilire quale delle due scale sia il meccanismo
clinicamente più rilevante.

## 11. Come si concatenano i passaggi

```text
Dati grezzi
  ↓
Campione bilanciato: 352 atleti × 2 periodi
  ↓
Descrizione di adoption, share, injury e burden
  ↓
Modello multilivello primario:
adoption + HSGP della share + covariate
  ↓
Contrasti primari population-standardized:
P50 vs non-user e P75 vs P25
  ↓
Curva continua per un profilo fisso:
letture illustrative a share 25%, 50%, 75%
  ↓
Controlli, sensibilità e poststratificazione
  ↓
Domanda interpretativa successiva:
share relativa o score assoluto?
  ↓
Due modelli lineari appaiati:
+24,1 punti percentuali di share vs +3,36 unità di score
  ↓
Risultato:
associazione competitiva simile sulle due scale
```

## 12. Cosa possiamo e non possiamo dire

Possiamo dire che:

- l'adoption da sola non mostra un'associazione chiara;
- tra gli utilizzatori, una maggiore esposizione competitiva è associata a
  maggiore prevalenza di infortunio;
- il segnale è presente sia sulla scala relativa sia su quella assoluta;
- la preparazione non mostra un'associazione chiara.

Non possiamo dire che:

- il carbonio causa gli infortuni;
- esista una soglia clinica a 25%, 50% o 75%;
- il 50% sia un change-point;
- la share sia superiore allo score assoluto;
- lo score assoluto sia il numero esatto di sedute;
- le curve siano già un calcolatore di rischio individuale validato.

Il passo successivo ideale è uno studio prospettico con ordine temporale noto e
misure oggettive di sedute uniche, durata, intensità e scarpa usata.
