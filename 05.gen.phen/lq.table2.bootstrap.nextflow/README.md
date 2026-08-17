# LQ table 2 CAAStools bootstrap Nextflow workflow

SLURM-ready CAAStools bootstrap workflow parallel to `../lq.table2.nextflow`.
RERconverge is intentionally absent until the complete gene-tree collection is
available.

## Design

The 100 unique four-FG versus four-BG longevity comparisons from
`../lq.table2.4species.comparisons/configurations/caas/` are represented as
100 rows in `inputs/longevity.100-comparisons.resampling.tsv`. Instead of
submitting one SLURM task for every gene/configuration pair, this workflow
submits one task per alignment and lets `ct bootstrap` evaluate all 100 trait
cycles inside that task.

For approximately 16,000 alignments, the design therefore requires roughly
16,000 SLURM tasks rather than 1.6 million tasks.

The bundled CAAStools copy supports the optional hypergeometric position
filter added for this analysis. Its default is configured as:

```groovy
caas_filter_significant = '0.05'
```

Set it to `'no'` in `conf/cluster.config` to scan all otherwise eligible
alignment positions.

## Inputs

- `inputs/longevity.100-comparisons.resampling.tsv`: the 100 longevity
  comparisons in CAAStools resampling format;
- `inputs/longevity.template.001.caas.cfg`: a four-FG/four-BG binary template
  used by CAAStools to initialize alignment slicing;
- `../lq.table2.nextflow/inputs/alignments/*.phy`: the existing alignment
  collection, reused read-only to avoid a second multi-gigabyte copy.

The resampling and template files can be regenerated deterministically from
the comparison configs:

```bash
python3 scripts/build_resampling.py
```

The generator validates the number of configs, their uniqueness, binary
format, group sizes, duplicate species and FG/BG overlap.

## Configure

All user-editable analytical and cluster settings are in
`conf/cluster.config`, including:

- alignment path;
- significance threshold and CAAStools gap/missing-data filters;
- SLURM partition, account, QoS and constraint;
- Conda initialization and environment name;
- process resources and optional pre-task commands.

The default partition is `std-cpu`. Each bootstrap task requests one CPU,
2 GB RAM and two hours. The driver requests one CPU and 16 GB RAM through the
`#SBATCH` header in `submit_pipeline_slurm.sh`.

## Conda environment

The workflow defaults to the same `phyloq` environment used by the discovery
workflow. To create or update it without removing any existing RERconverge
packages:

```bash
bash create_conda_environment.sh
```

The bootstrap workflow itself requires only Nextflow, Java, Python,
Biopython, NumPy and SciPy.

## Launch

From this directory on the SLURM login node:

```bash
sbatch submit_pipeline_slurm.sh
```

Resume the latest bootstrap run:

```bash
sbatch submit_pipeline_slurm.sh -resume
```

The run identifier is stored in the Git-ignored `.last_run_id`. Work/cache,
results and generated logs are also Git-ignored. Reports, timelines and large
trace tables are not generated.

## Outputs

Each alignment produces one table:

```text
results/RUN_ID/caas-bootstrap/GENE.bootstrap.caas.tsv
```

Each row contains the alignment position, number of resampling cycles with a
CAAS, total cycles, bootstrap value, positive cycle identifiers and the
template config path, following the native CAAStools bootstrap format.
