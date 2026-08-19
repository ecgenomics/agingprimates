# Parametri e decisioni dell'analisi CAAStools `approach.pss`

## Obiettivo

Eseguire CAAStools sugli stessi allineamenti proteici delle altre analisi di
`05.gen.phen`, confrontando specie con longevity quotient alto (foreground,
`1`) e basso (background, `0`). La run attiva contiene due ipotesi PSS in
formato 4-vs-4.

Il fenotipo usato è `LQ_mammal` non trasformato. In ogni coppia PSS, la specie
con LQ maggiore è assegnata al foreground e quella con LQ minore al background.

## Decisione attiva: analisi 4-vs-4

Le precedenti configurazioni 5-vs-5 sono risultate restrittive. La nuova run
mantiene soltanto i primi due criteri PSS e usa esattamente quattro FG e quattro
BG per ciascun confronto:

- `configs/01.max-delta-lq.4v4.caas.cfg`;
- `configs/02.max-pss-score.4v4.caas.cfg`.

Per la massima ΔLQ vengono mantenute le coppie di rango 1, 2, 3 e 5. La coppia
di rango 4 (`Leontopithecus_rosalia` vs `Mico_humeralifer`) è esclusa perché
riutilizza `Mico_humeralifer`, già BG nel rango 1, e produrrebbe soltanto tre BG
distinti. In questo modo ogni specie FG conserva la propria controparte BG
della coppia originale.

Per il massimo `FinalScore` vengono mantenute direttamente le coppie di rango
1–4.

Il manifest numerico della selezione attiva è
`inputs/pss_pair_manifest.4v4.tsv`. Le configurazioni precedenti e i relativi
manifest restano nel repository come provenienza, ma il glob attivo è
`configs/*.4v4.caas.cfg`. In particolare,
`configs/03.absolute-lq-extremes.caas.cfg` è esplicitamente escluso.

### Ripresa della run esistente

`sbatch submit_pipeline_slurm.sh -resume` riusa il `RUN_ID` registrato in
`.last_run_id` e la cache in `work/`. I nomi nuovi dei config 4-vs-4 fanno sì
che Nextflow generi task distinti senza confonderli con le vecchie analisi.
I due file complessivi in `results/RUN_ID/merged/` vengono rigenerati usando
soltanto i due config attivi. Le vecchie cartelle per-config già pubblicate non
vengono cancellate automaticamente e vanno considerate risultati storici.

## Dati di origine

- Tabella fenotipica:
  `03.phenotype.shift.detection.lq/input/longevity_quotient.tsv`.
- Numero di specie nella tabella fenotipica: 162.
- Colonna fenotipica: `LQ_mammal`.
- Risultati PSS:
  `03.phenotype.shift.detection.lq/results/LQ.score_results.tsv`.
- Colonna dello score aggregato PSS: `FinalScore`.
- Numero totale di coppie PSS: 7.875.
- Top 1% PSS: 79 coppie, ottenute come `ceil(7875 × 0,01)` dopo ordinamento
  decrescente per `FinalScore`.
- Allineamenti:
  `05.gen.phen/lq.table2.nextflow/inputs/alignments/*.phy`.
- Copia di CAAStools:
  `05.gen.phen/lq.table2.nextflow/bin/caastools`.

## Ipotesi 1: massima differenza di LQ nel top 1% PSS

Config storico: `configs/01.max-delta-lq.caas.cfg`.

Procedura:

1. selezione delle 79 coppie con `FinalScore` più alto;
2. calcolo di `abs(TraitValue1 - TraitValue2)`;
3. selezione delle cinque coppie con differenza assoluta maggiore;
4. assegnazione della specie con LQ maggiore a FG e di quella con LQ minore a
   BG;
5. deduplicazione delle specie ripetute nel config.

| Rango | FG (`1`) | LQ | BG (`0`) | LQ | ΔLQ | FinalScore |
| ---: | --- | ---: | --- | ---: | ---: | ---: |
| 1 | `Cebus_olivaceus` | 1,612012 | `Mico_humeralifer` | 0,740275 | 0,871738 | 2,095528 |
| 2 | `Hylobates_muelleri` | 1,830577 | `Hylobates_klossii` | 1,111756 | 0,718821 | 16,152878 |
| 3 | `Cheirogaleus_medius` | 1,355345 | `Cheirogaleus_major` | 0,638673 | 0,716672 | 3,306366 |
| 4 | `Leontopithecus_rosalia` | 1,440173 | `Mico_humeralifer` | 0,740275 | 0,699899 | 2,846785 |
| 5 | `Lemur_catta` | 1,344976 | `Prolemur_simus` | 0,647946 | 0,697030 | 2,608588 |

`Mico_humeralifer` compare in due coppie e viene scritto una sola volta. Il
config finale contiene quindi 5 FG e 4 BG.

La provenienza numerica completa è in `inputs/pss_pair_manifest.tsv`.

## Ipotesi 2: massimo FinalScore PSS

Config storico: `configs/02.max-pss-score.caas.cfg`.

Procedura:

1. ordinamento di tutte le 7.875 coppie per `FinalScore` decrescente;
2. selezione delle prime cinque coppie;
3. assegnazione della specie con LQ maggiore a FG e di quella con LQ minore a
   BG.

| Rango | FG (`1`) | LQ | BG (`0`) | LQ | ΔLQ | FinalScore |
| ---: | --- | ---: | --- | ---: | ---: | ---: |
| 1 | `Leontopithecus_rosalia` | 1,440173 | `Leontopithecus_chrysomelas` | 0,988980 | 0,451194 | 31,164498 |
| 2 | `Saguinus_oedipus` | 1,298915 | `Saguinus_geoffroyi` | 0,982658 | 0,316257 | 29,184455 |
| 3 | `Nycticebus_coucang` | 1,172397 | `Nycticebus_bengalensis` | 1,022506 | 0,149891 | 17,996283 |
| 4 | `Callithrix_jacchus` | 1,184560 | `Callithrix_kuhlii` | 1,005881 | 0,178679 | 17,543259 |
| 5 | `Hylobates_muelleri` | 1,830577 | `Hylobates_klossii` | 1,111756 | 0,718821 | 16,152878 |

Il config contiene 5 FG e 5 BG. La provenienza numerica completa è in
`inputs/pss_pair_manifest.tsv`.

## Ipotesi 3: estremi assoluti di LQ

Config storico, escluso dalla run attiva:
`configs/03.absolute-lq-extremes.caas.cfg`.

Gruppo BG: le cinque specie con `LQ_mammal` più basso nell'intera tabella
fenotipica.

Gruppo FG: le specie vengono ordinate per `LQ_mammal` decrescente e vengono
selezionate le prime cinque appartenenti a cinque generi distinti. Per ogni
genere viene quindi mantenuta soltanto la specie con LQ maggiore.

| Gruppo | Rango nel gruppo | Rango nell'estremo di origine | Specie | LQ |
| --- | ---: | ---: | --- | ---: |
| FG | 1 | 1° più alto | `Hylobates_muelleri` | 1,830577 |
| FG | 2 | 2° più alto | `Cebus_capucinus` | 1,825399 |
| FG | 3 | 5° più alto | `Sapajus_apella` | 1,564604 |
| FG | 4 | 7° più alto | `Aotus_lemurinus` | 1,447364 |
| FG | 5 | 8° più alto | `Leontopithecus_rosalia` | 1,440173 |
| BG | 1 | 1° più basso | `Lepilemur_mustelinus` | 0,504586 |
| BG | 2 | 2° più basso | `Presbytis_melalophos` | 0,587634 |
| BG | 3 | 3° più basso | `Nasalis_larvatus` | 0,624430 |
| BG | 4 | 4° più basso | `Propithecus_diadema` | 0,624914 |
| BG | 5 | 5° più basso | `Cheirogaleus_major` | 0,638673 |

### Decisione: una specie per genere nel FG

È stato deciso di applicare al gruppo FG un vincolo generale di una sola specie
per genere, mantenendo la specie con LQ maggiore. La graduatoria originale
viene quindi attraversata in ordine decrescente:

1. `Hylobates_muelleri` viene mantenuta per il genere `Hylobates`;
2. `Cebus_capucinus` viene mantenuta per il genere `Cebus`;
3. `Cebus_olivaceus` viene esclusa perché `Cebus` è già rappresentato;
4. `Hylobates_lar` viene esclusa perché `Hylobates` è già rappresentato;
5. `Sapajus_apella` viene mantenuta;
6. `Hylobates_agilis` viene esclusa perché `Hylobates` è già rappresentato;
7. `Aotus_lemurinus` viene mantenuta;
8. `Leontopithecus_rosalia` viene mantenuta e completa il gruppo di cinque FG.

Il gruppo BG non contiene generi ripetuti e non richiede sostituzioni.

### Disponibilità genomica

`Cebus_capucinus`, `Lepilemur_mustelinus` e `Presbytis_melalophos` hanno
`genome_available=False`, zero campioni genomici e non compaiono negli
allineamenti locali. Sono mantenute nel config per rappresentare esattamente la
selezione fenotipica richiesta, ma con gli allineamenti disponibili i gruppi
effettivi possono contenere al massimo 4 FG e 3 BG.

`Aotus_lemurinus` ha quattro campioni genomici ed è marcata come disponibile
nella tabella fenotipica, ma non compare nel sottoinsieme locale di allineamenti
attualmente presente. Nel controllo locale il FG osservabile può quindi
ridursi ulteriormente a 3 specie; sul cluster dipenderà dagli allineamenti
effettivamente disponibili per ciascun gene.

La tabella completa con rango, LQ e disponibilità è in
`inputs/lq_extremes_manifest.tsv`.

## Parametri CAAStools

- Modalità: `discovery`.
- Formato allineamento: `phylip-relaxed`.
- Pattern analizzati: `1,2,3`.
- `max_bg_gaps`: `NO`.
- `max_fg_gaps`: `NO`.
- `max_gaps`: `NO`.
- `max_gaps_per_position`: `0.5`.
- `max_bg_miss`: `NO`.
- `max_fg_miss`: `NO`.
- `max_miss`: `NO`.
- Filtro dei risultati: `Pvalue <= 0.05`, inclusivo.
- Il `Pvalue` è quello ipergeometrico grezzo prodotto da CAAStools.
- Non viene applicata alcuna correzione FDR.

## Workflow Nextflow e risorse SLURM

- Nextflow DSL: 2.
- Executor: SLURM.
- Partizione predefinita: `std-cpu`.
- Numero massimo di job in coda: 100.
- Frequenza massima di invio: un batch ogni 10 secondi.
- Directory di lavoro: `approach.pss/work`.
- Identificatore della run: `run-YYMMDDhhmmss`.
- `-resume` riusa l'identificatore conservato in `.last_run_id`.

Processo CAAStools per gene/config:

- 1 CPU;
- 2 GB RAM;
- tempo massimo 10 minuti;
- nessun retry;
- `errorStrategy = ignore`: timeout ed errori non terminano l'intera run.

Processi di merge:

- 1 CPU;
- 4 GB RAM;
- tempo massimo 2 ore.

Driver Nextflow:

- 1 CPU;
- 16 GB RAM;
- tempo massimo 7 giorni;
- JVM: `-Xms512m -Xmx8g` salvo override di `NXF_OPTS`.

## Ambiente Conda

- Nome predefinito: `genphen-caas-pss`.
- OpenJDK 17.
- Nextflow.
- Python 3.11.
- `setuptools < 82`, necessario per l'importazione di `pkg_resources` da parte
  di CAAStools.
- Biopython.
- NumPy.
- SciPy.

Il percorso predefinito dello script di inizializzazione Conda del cluster è
`/homes/aplic/noarch/software/Miniconda3/23.9.0-0/etc/profile.d/conda.sh`.
Può essere sostituito tramite `GENPHEN_CONDA_SH`; il nome dell'ambiente può
essere sostituito tramite `GENPHEN_CONDA_ENV`.

## Merge e risultati

Per ogni gene/config, un risultato vuoto viene convertito in una tabella con il
solo header. I job falliti o terminati per timeout non producono un file e sono
omessi dal merge senza interrompere il workflow.

Il workflow produce:

1. un merge per ciascuno dei due config 4-vs-4 attivi;
2. un merge complessivo con un solo header;
3. un file complessivo filtrato per `Pvalue <= 0.05`.

Output finali in `results/RUN_ID/merged/`:

- `caas.all-configs.all-genes.tsv`;
- `caas.all-configs.all-genes.significant.tsv`.

## Controlli automatici

I test verificano:

- ricostruzione delle due selezioni PSS dai risultati originali;
- assegnazione FG/BG in base al valore LQ;
- deduplicazione delle specie ripetute nei config PSS;
- selezione degli estremi LQ con il vincolo di una specie per genere nel FG;
- corrispondenza fra config e manifest;
- presenza di esattamente 4 FG e 4 BG in ciascun config attivo;
- esclusione dei config storici dal glob attivo;
- merge con un solo header;
- filtro inclusivo a `Pvalue <= 0.05`.

I controlli locali non eseguono CAAStools sull'intero dataset e non inviano job
al cluster.
