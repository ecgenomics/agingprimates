# FQ: CAAStools pooled discovery

Copia del workflow LQ locale
`../../05.gen.phen/caas/lq.table2.bootstrap.nextflow`, con lo stesso codice
CAAStools, gli stessi due processi Nextflow, parametri analitici, risorse Slurm,
ambiente Conda e comportamento di ripresa/gestione degli errori.

## Grouping già definito

Il pool deriva esclusivamente da
`../03.phenotype.shift.detection.fq/results/FQ.shared_tail_ranking.tsv`,
copiato in `inputs/` per tracciabilità. Sono le specie comuni al top e bottom
1% delle 4.560 coppie PSS (46 coppie per coda), classificate con le soglie FQ
già usate nel progetto: FG > 1.3 e BG < 0.8. Le cinque specie intermedie sono
escluse. Non sono usate le selezioni alternative per clade o per coppia.

- FG (5): Leontocebus_fuscicolis, Galago_senegalensis,
  Trachypithecus_phayrei, Callimico_goeldii, Cercopithecus_diana.
- BG (7): Cebus_olivaceus, Cercopithecus_cephus, Ateles_belzebuth,
  Ateles_paniscus, Ateles_geoffroyi, Pongo_pygmaeus, Saimiri_sciureus.

`inputs/fertility.full-pools.caas.cfg` contiene `species<TAB>1/0`, senza
intestazione. I nomi delle specie sono conservati esattamente come nel grouping.

## Identità con LQ

`bin/caastools/` e `scripts/prepare_pooled_hypotheses.py` sono copie byte per
byte del riferimento LQ; sono escluse solo le cache Python. Anche gli script
di lancio/creazione ambiente e `environment.yml` sono copiati da LQ.
`LQ_COPY_AUDIT.tsv` registra SHA-256 originali e locali.

Le sole modifiche ai file operativi sono:

- `longevity.` diventa `fertility.` nei nomi del pool e dei metadati;
- il percorso relativo degli allineamenti raggiunge la stessa raccolta LQ
  dalla nuova directory;
- il messaggio dello script Slurm indica la nuova directory;
- il commento sul massimo numero di confronti indica 175 anziché 525.

Restano invariati: 100 ipotesi 4 FG contro 4 BG, seed 260811, selezione
deterministica SHA-256, formato `phylip-relaxed`, filtro posizionale 0.05,
pattern `1,2,3`, almeno 3 specie osservate per gruppo, limite di gap per
posizione 0.5 e gli altri limiti di gap/specie mancanti impostati a `NO`.
Con questo pool sono possibili `choose(5,4) * choose(7,4) = 175` confronti;
il default resta 100, come in LQ.

## Lancio sul cluster

In `conf/cluster.config` verificare il percorso della raccolta completa degli
allineamenti e le impostazioni del cluster. Il percorso predefinito riutilizza
la raccolta LQ locale, che qui contiene soltanto cinque allineamenti di esempio.
Il codice CAAStools è incluso localmente in questa cartella.

```bash
cd 04.caas.pooled.nextflow
bash create_conda_environment.sh  # solo se occorre predisporre phyloq
sbatch submit_pipeline_slurm.sh
```

Ripresa della stessa analisi:

```bash
sbatch submit_pipeline_slurm.sh -resume
```

Come in LQ, il driver e i task usano l'ambiente `phyloq`, la partizione
`std-cpu` e l'inizializzazione Conda Correfoc. La preparazione richiede
1 CPU/1 GB/10 minuti; ogni gene 1 CPU/2 GB/30 minuti; concorrenza 100 task.
Un gene fallito viene ignorato e non produce risultati; un errore nella
preparazione interrompe il workflow. Controllare quindi anche il log Nextflow.

## Output

```text
results/RUN_ID/metadata/fertility.pooled.hypotheses.tsv
results/RUN_ID/metadata/fertility.pooled.hypotheses.metadata.tsv
results/RUN_ID/caas-pooled/GENE.pooled.caas.tsv
results/RUN_ID/caas-pooled-events/GENE.pooled.caas.events.tsv
```

Il workflow genera una sola tabella condivisa di ipotesi per run e lancia
`ct pooled-discovery` una volta per allineamento, come in LQ. Il pooling
conserva i pool completi come denominatori degli eventi. Questo workflow
replica il lancio pooled; non aggiunge analisi successive o nuovi filtri.
