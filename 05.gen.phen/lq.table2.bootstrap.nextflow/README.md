# LQ table 2 CAAStools pooled-discovery Nextflow workflow

SLURM-ready workflow for species- and amino-acid-aware CAAStools pooled
discovery. The directory name is retained for continuity with the earlier
bootstrap implementation. RERconverge remains intentionally absent until the
complete gene-tree collection is available.

## Design

The analysis starts from one complete fixed-side longevity config containing
seven FG and six BG species. For four-species subsets there are:

```text
choose(7, 4) * choose(6, 4) = 35 * 15 = 525
```

possible unique four-vs-four comparisons. The workflow is configured to select
100 of them:

```groovy
pooled_fg_size = 4
pooled_bg_size = 4
pooled_comparisons = '100'
pooled_seed = 260811
```

Set `pooled_comparisons = 'max'` to use all 525. A selection smaller than the
maximum is determined by a seeded SHA-256 ranking, making it reproducible and
independent of Python's random-number implementation. It is not expected to
reproduce a historical selection created with R's `sample()` from the same
numeric seed.

The Nextflow graph has two stages:

1. `PREPARE_POOLED_HYPOTHESES` runs once, generates the shared headerless
   resampling-format table, calculates its SHA-256 checksum, and publishes the
   hypotheses plus metadata;
2. `CAASTOOLS_POOLED_DISCOVERY` runs once per alignment and reuses that exact
   shared table. CAAStools validates it against the complete fixed FG/BG pools
   before scanning the alignment.

The first stage prevents approximately 16,000 gene jobs from independently
materializing the same hypotheses. The second stage still generates only one
SLURM task per alignment, rather than one task per alignment/hypothesis pair.

CAAStools retains the complete seven-FG/six-BG pools as event denominators even
though individual hypotheses contain four species per side. Positive
hypotheses are merged into amino-acid-compatible events, species are counted
once per event, incompatible signatures remain separate, and the final table
reports event pattern, dominant amino acids, support fractions, conflicting
species, nominal event p-value, and positional p-value.

## Analytical filters

The optional positional hypergeometric prefilter remains configured as:

```groovy
caas_filter_significant = '0.05'
```

Set it to `'no'` to scan all otherwise eligible positions. This is distinct
from the nominal event p-value written to the event output.

Each four-vs-four hypothesis requires at least three observed species in each
group:

```groovy
caas_min_fg_observed = 3
caas_min_bg_observed = 3
```

Observed means present in the alignment and non-gapped at the position. A
hypothesis with fewer than three usable species on either side cannot be
positive. The broader gap and missing-species limits remain `NO`.

## Inputs

- `inputs/longevity.full-pools.caas.cfg`: complete fixed seven-FG/six-BG
  discovery pools, in `species<TAB>1/0` format;
- `../lq.table2.nextflow/inputs/alignments/*.phy`: existing alignment
  collection, reused read-only;
- `bin/caastools/`: bundled local CAAStools implementation;
- `scripts/prepare_pooled_hypotheses.py`: thin run-level wrapper around the
  CAAStools pooling function.

The earlier `longevity.100-comparisons.resampling.tsv` and four-vs-four
template remain historical inputs but are not used by this workflow.

## Configure

All user-editable analytical and cluster settings are in
`conf/cluster.config`, including:

- alignment and complete-pool paths;
- FG/BG subset sizes, comparison count, and selection seed;
- positional significance, minimum observed-species, gap, and missing-data
  filters;
- SLURM partition, account, QoS, constraint, memory, and time;
- Conda initialization and environment name.

The default partition is `std-cpu`. The one-time preparation task requests one
CPU, 1 GB RAM and ten minutes. Each gene task requests one CPU, 2 GB RAM and
30 minutes. At most 100 tasks are submitted concurrently. A failed or
time-limited gene is ignored by Nextflow so that the remaining alignments can
finish; the failed gene has no result and can be identified in the Nextflow
log. Failure of the shared hypothesis-preparation task still terminates the
workflow. The driver requests one CPU and 16 GB RAM through
`submit_pipeline_slurm.sh`.

## Conda environment

The workflow uses the `phyloq` environment. To create or update it:

```bash
bash create_conda_environment.sh
```

The pooled workflow requires Nextflow, Java, Python, Biopython, NumPy, SciPy,
and DendroPy. The environment file includes these dependencies.

## Launch

From this directory on the SLURM login node:

```bash
sbatch submit_pipeline_slurm.sh
```

This change introduces new Nextflow processes and output names, so the first
pooled-discovery analysis should be launched as a fresh run, without
`-resume`. Resume that new run later with:

```bash
sbatch submit_pipeline_slurm.sh -resume
```

The run identifier is stored in the Git-ignored `.last_run_id`. Work/cache,
results and logs are also Git-ignored. HTML reports, timelines, and large trace
tables are not generated.

## Outputs

The run-level preparation stage publishes:

```text
results/RUN_ID/metadata/longevity.pooled.hypotheses.tsv
results/RUN_ID/metadata/longevity.pooled.hypotheses.metadata.tsv
```

The metadata records complete-pool sizes, hypothesis sizes, maximum and
selected comparison counts, seed, selection method, and hypotheses SHA-256.

Each alignment produces a backwards-compatible position table:

```text
results/RUN_ID/caas-pooled/GENE.pooled.caas.tsv
```

and a headered species/amino-acid event table:

```text
results/RUN_ID/caas-pooled-events/GENE.pooled.caas.events.tsv
```

The historical table retains position, number of positive hypotheses, total
hypotheses, frequency, hypothesis identifiers, positional p-value, and config
path. The event table is the primary biological result; raw hypothesis counts
are retained for traceability rather than interpreted as independent evidence.
