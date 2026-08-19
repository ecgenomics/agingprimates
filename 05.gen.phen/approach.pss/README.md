# PSS- and LQ-guided CAAStools workflow (`approach.pss`)

SLURM-ready Nextflow workflow that tests three alternative hypotheses: two
derived from the pairwise phenotype shift score (PSS) and one from the absolute
extremes of the untransformed longevity quotient (`LQ_mammal`). Higher-LQ
species are foreground (`1`) and lower-LQ species are background (`0`).

## Hypotheses and CAAStools configurations

The PSS table contains 7,875 pairs; the top 1% therefore contains 79 pairs.
Selections were calculated from
`../../03.phenotype.shift.detection.lq/results/LQ.score_results.tsv`; the LQ
values originate from
`../../03.phenotype.shift.detection.lq/input/longevity_quotient.tsv`.

### 1. Largest phenotype differences within the PSS top 1%

Configuration: `configs/01.max-delta-lq.caas.cfg`

The five top-1% PSS pairs with the largest absolute difference in `LQ_mammal`
were selected. `Mico_humeralifer` is the lower-LQ member of two pairs, so the
deduplicated CAAStools configuration contains 5 FG and 4 BG species.

| FG (higher LQ) | BG (lower LQ) | Absolute LQ difference |
| --- | --- | ---: |
| `Cebus_olivaceus` | `Mico_humeralifer` | 0.871738 |
| `Hylobates_muelleri` | `Hylobates_klossii` | 0.718821 |
| `Cheirogaleus_medius` | `Cheirogaleus_major` | 0.716672 |
| `Leontopithecus_rosalia` | `Mico_humeralifer` | 0.699899 |
| `Lemur_catta` | `Prolemur_simus` | 0.697030 |

### 2. Highest PSS FinalScore

Configuration: `configs/02.max-pss-score.caas.cfg`

The five pairs with the highest PSS `FinalScore` were selected. This
configuration contains 5 FG and 5 BG species.

| FG (higher LQ) | BG (lower LQ) | PSS FinalScore |
| --- | --- | ---: |
| `Leontopithecus_rosalia` | `Leontopithecus_chrysomelas` | 31.164498 |
| `Saguinus_oedipus` | `Saguinus_geoffroyi` | 29.184455 |
| `Nycticebus_coucang` | `Nycticebus_bengalensis` | 17.996283 |
| `Callithrix_jacchus` | `Callithrix_kuhlii` | 17.543259 |
| `Hylobates_muelleri` | `Hylobates_klossii` | 16.152878 |

The full numerical provenance, including both LQ values, absolute differences,
PSS scores and criterion-specific ranks, is recorded in
`inputs/pss_pair_manifest.tsv`.

### 3. Five highest versus five lowest absolute LQ values

Configuration: `configs/03.absolute-lq-extremes.caas.cfg`

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
- Configurations: all three `configs/*.caas.cfg` files.
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
2. merges the three hypothesis tables into one table;
3. retains exactly one header;
4. writes a second table containing rows with raw `Pvalue <= 0.05`.

Final files are written to `results/RUN_ID/merged/`:

- `caas.all-configs.all-genes.tsv`;
- `caas.all-configs.all-genes.significant.tsv`.

The per-hypothesis merges are under `results/RUN_ID/merged/by-config/`, while
individual gene results are under
`results/RUN_ID/caas/by-config/CONFIG_ID/`.

## Local checks

The lightweight checks do not run CAAStools or submit cluster jobs:

```bash
python3 -m unittest discover -s tests -v
bash -n create_conda_environment.sh run_pipeline.sh submit_pipeline_slurm.sh
```
