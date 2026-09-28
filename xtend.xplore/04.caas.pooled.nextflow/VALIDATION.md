# Verifica locale — 2026-09-28

- Grouping verificato direttamente sui 4.560 punteggi PSS: 46 coppie per
  coda, 17 specie nell'intersezione; gli Status salvati corrispondono alle
  soglie esistenti, producendo 5 FG e 7 BG.
- Verificati SHA-256 dei 30 file copiati dal workflow LQ: 26 identici;
  i quattro file adattati sono `main.nf`, `nextflow.config`,
  `conf/cluster.config` e `submit_pipeline_slurm.sh`. Nessuna modifica al
  codice CAAStools o allo script Python di preparazione.
- Sintassi dei tre script shell e caricamento della configurazione Slurm
  tramite `nextflow -c conf/cluster.config config -flat`: OK.
- Test reale Nextflow 26.04.6 con executor locale su
  `C4BPA.Homo_sapiens.filter2.phy` della raccolta LQ: entrambi i processi
  completati con exit code 0. Nessun job Slurm sottomesso.
- Ambiente temporaneo `/tmp/fq-caas-validation-env`: Python di sistema,
  librerie scientifiche già installate e DendroPy 5.1.0. La configurazione
  Conda del workflow è rimasta identica a LQ.
- L'allineamento C4BPA contiene 230 specie, incluse tutte le 12 del pool.
- Metadati: 175 confronti possibili, 100 selezionati, 4 FG/4 BG per
  ipotesi, seed 260811. SHA-256 della tabella di ipotesi:
  `bd7bf4faab6fa8b55fed57d72552171120199f419582ec12ec2ef1da6c42867c`.
- Output del test: 334 righe nella tabella pooled e 19 eventi più
  intestazione nella tabella eventi.

Comando eseguito dalla directory del workflow (percorso degli allineamenti
abbreviato qui con il relativo equivalente):

```bash
nextflow -log /tmp/fq-caas-smoke.log run main.nf \
  --run_id validation \
  --alignments '../../05.gen.phen/caas/lq.table2.nextflow/inputs/alignments/C4BPA.Homo_sapiens.filter2.phy' \
  --python_command /tmp/fq-caas-validation-env/bin/python \
  --results_root /tmp/fq-caas-validation-results \
  -work-dir /tmp/fq-caas-validation-work -ansi-log false
```

Risultati del test in `/tmp/fq-caas-validation-results/validation`.
Questa verifica riguarda un singolo gene e il funzionamento del workflow;
non costituisce il lancio della raccolta completa sul cluster.
