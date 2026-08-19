# PSS- and LQ-guided CAAStools workflow (`approach.pss`)

SLURM-ready Nextflow workflow that currently tests two 4-vs-4 hypotheses
derived from the pairwise phenotype shift score (PSS). Higher-LQ species are
foreground (`1`) and lower-LQ species are background (`0`). The older 5-vs-5
and absolute-LQ configurations remain in `configs/` only as analytical
provenance and are excluded from the active Nextflow input glob.

## Hypotheses and CAAStools configurations

The PSS table contains 7,875 pairs; the top 1% therefore contains 79 pairs.
Selections were calculated from
`../../03.phenotype.shift.detection.lq/results/LQ.score_results.tsv`; the LQ
values originate from
`../../03.phenotype.shift.detection.lq/input/longevity_quotient.tsv`.

### 1. Largest phenotype differences within the PSS top 1%

Active configuration: `configs/01.max-delta-lq.4v4.caas.cfg`

Four complete pairs were selected to obtain exactly 4 FG and 4 BG species.
Among the five top-1% PSS pairs with the largest absolute difference in
`LQ_mammal`, rank 4 was skipped because it repeats `Mico_humeralifer` as BG;
rank 5 was retained instead.

| FG (higher LQ) | BG (lower LQ) | Absolute LQ difference |
| --- | --- | ---: |
| `Cebus_olivaceus` | `Mico_humeralifer` | 0.871738 |
| `Hylobates_muelleri` | `Hylobates_klossii` | 0.718821 |
| `Cheirogaleus_medius` | `Cheirogaleus_major` | 0.716672 |
| `Lemur_catta` | `Prolemur_simus` | 0.697030 |

### 2. Highest PSS FinalScore

Active configuration: `configs/02.max-pss-score.4v4.caas.cfg`

The four pairs with the highest PSS `FinalScore` were selected. This
configuration contains 4 FG and 4 BG species.

| FG (higher LQ) | BG (lower LQ) | PSS FinalScore |
| --- | --- | ---: |
| `Leontopithecus_rosalia` | `Leontopithecus_chrysomelas` | 31.164498 |
| `Saguinus_oedipus` | `Saguinus_geoffroyi` | 29.184455 |
| `Nycticebus_coucang` | `Nycticebus_bengalensis` | 17.996283 |
| `Callithrix_jacchus` | `Callithrix_kuhlii` | 17.543259 |

The full numerical provenance, including both LQ values, absolute differences,
PSS scores and criterion-specific ranks, is recorded in
`inputs/pss_pair_manifest.4v4.tsv`. The former five-pair selection remains in
`inputs/pss_pair_manifest.tsv` for provenance.

### Historical analysis excluded from the active run

The former absolute-LQ configuration,
`configs/03.absolute-lq-extremes.caas.cfg`, is retained but is not matched by
the active `configs/*.4v4.caas.cfg` glob and will not be submitted.

| FG: five highest LQ | LQ | BG: five lowest LQ | LQ |
| --- | ---: | --- | ---: |
| `Hylobates_muelleri` | 1.830577 | `Lepilemur_mustelinus` | 0.504586 |
| `Cebus_capucinus` | 1.825399 | `Presbytis_melalophos` | 0.587634 |
| `Sapajus_apella` | 1.564604 | `Nasalis_larvatus` | 0.624430 |
| `Aotus_lemurinus` | 1.447364 | `Propithecus_diadema` | 0.624914 |
| `Leontopithecus_rosalia` | 1.440173 | `Cheirogaleus_major` | 0.638673 |

The FG group is obtained by scanning species in descending LQ order and
retaining only the first species from each genus. This removes
`Cebus_olivaceus`, `Hylobates_lar` and `Hylobates_agilis` because their genera
are already represented by species with higher LQ. The selected global ranks
are 1, 2, 5, 7 and 8.

`Cebus_capucinus`, `Lepilemur_mustelinus` and `Presbytis_melalophos` have
`genome_available=False` in the phenotype table and are absent from the local
alignment collection. They remain in the requested 5-vs-5 configuration for
analytical provenance, but the effective groups in the available alignments
can be at most 4 FG and 3 BG. `Aotus_lemurinus` has genomic samples but is also
absent from the local alignment subset, which can further reduce the observed
FG to 3 species locally. These availability data are recorded in
`inputs/lq_extremes_manifest.tsv`.

All selection rules, CAAStools parameters, cluster resources and decision
records are collected in `ANALYSIS_PARAMETERS.md`.

## Inputs and analytical settings

- Alignments: the same relaxed-PHYLIP collection used by the other runs,
  configured as `../lq.table2.nextflow/inputs/alignments/*.phy`.
- CAAStools: the same versioned copy at
  `../lq.table2.nextflow/bin/caastools`.
- Configurations: only the two `configs/*.4v4.caas.cfg` files.
- Discovery settings: patterns 1, 2 and 3; no explicit per-group gap or
  missing-species maximum; maximum total gap ratio 0.5.
- Per-job timeout: SLURM stops CAAStools after 10 minutes. Timed-out or failed
  jobs are ignored by Nextflow, so the remaining jobs and merges continue.
- Significance output: raw CAAStools hypergeometric `Pvalue <= 0.05`; no FDR
  correction is applied.

Paths, resources and thresholds can be changed in `conf/cluster.config`.

## Conda environment

On the cluster login node, create or update the dedicated environment once:

```bash
bash create_conda_environment.sh
```

This creates `genphen-caas-pss` with Nextflow, OpenJDK 17, Python 3.11,
Biopython, NumPy and SciPy, and verifies both CAAStools and Nextflow. Override
the cluster Conda installation or environment name with `GENPHEN_CONDA_SH` and
`GENPHEN_CONDA_ENV` if necessary.

## Cluster launch

Submit from this directory:

```bash
sbatch submit_pipeline_slurm.sh
```

Resume the latest run, preserving its timestamped run ID, with:

```bash
sbatch submit_pipeline_slurm.sh -resume
```

After pulling these changes on the cluster, that command reuses the previous
Nextflow work directory and run ID. Because the two 4-vs-4 configs have new
filenames and contents, their discovery tasks are new; the historical configs
are no longer supplied as workflow inputs.

The Nextflow driver and every worker activate the dedicated Conda environment.
Individual CAAStools jobs request 1 CPU, 2 GB and 10 minutes. A timeout or
other worker failure is ignored rather than terminating the workflow. Merge
jobs request 1 CPU, 4 GB and 2 hours.

## Outputs

Each successful gene/config task emits a headered result, including a
header-only result when no CAAS is found. A timed-out or failed job has no
result and is omitted from the merge, but does not terminate the workflow.
After all discovery tasks have completed or been ignored, Nextflow:

1. merges all genes separately for each hypothesis;
2. merges the two active 4-vs-4 hypothesis tables into one table;
3. retains exactly one header;
4. writes a second table containing rows with raw `Pvalue <= 0.05`.

Final files are written to `results/RUN_ID/merged/`:

- `caas.all-configs.all-genes.tsv`;
- `caas.all-configs.all-genes.significant.tsv`.

The per-hypothesis merges are under `results/RUN_ID/merged/by-config/`, while
individual gene results are under
`results/RUN_ID/caas/by-config/CONFIG_ID/`.

When `-resume` reuses an existing `RUN_ID`, the two generic final tables above
are overwritten with the new two-config merge. Previously published
per-config folders are not automatically deleted, so old 5-vs-5 or
absolute-LQ files may remain alongside the new folders; they are historical
files and are not included in the regenerated generic merge.

## Local checks

The lightweight checks do not run CAAStools or submit cluster jobs:

```bash
python3 -m unittest discover -s tests -v
bash -n create_conda_environment.sh run_pipeline.sh submit_pipeline_slurm.sh
```
