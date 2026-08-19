# Alternative CAAStools workflow (`approach.bis`)

SLURM-ready Nextflow workflow that runs CAAStools discovery for every protein
alignment against the four genus-reduced configurations in `configs/`.

## Inputs and analytical settings

- Alignments: the same relaxed-PHYLIP collection used by
  `../lq.table2.nextflow`, configured as
  `../lq.table2.nextflow/inputs/alignments/*.phy`.
- CAAStools: the same versioned copy at
  `../lq.table2.nextflow/bin/caastools`.
- Configurations: all four `configs/*.caas.cfg` files.
- Discovery settings: patterns 1, 2 and 3; no explicit per-group gap or
  missing-species maximum; maximum total gap ratio 0.5.
- Per-job timeout: SLURM stops CAAStools after 10 minutes. Timed-out or failed
  jobs are ignored by Nextflow, so the remaining jobs and merges continue.
- Significance: raw CAAStools hypergeometric `Pvalue <= 0.05`.

The paths and thresholds can be changed in `conf/cluster.config`.

## Conda environment

On the cluster login node, create or update the dedicated environment once:

```bash
bash create_conda_environment.sh
```

This creates `genphen-caas-bis` with Nextflow, OpenJDK 17, Python 3.11,
Biopython, NumPy and SciPy, and verifies both CAAStools and Nextflow. Override
the cluster Conda installation or environment name with `GENPHEN_CONDA_SH` and
`GENPHEN_CONDA_ENV` if necessary.

## Cluster launch

Submit from this directory:

```bash
sbatch submit_pipeline_slurm.sh
```

Resume the latest run with the same timestamped run ID:

```bash
sbatch submit_pipeline_slurm.sh -resume
```

The Nextflow driver and every worker activate the dedicated Conda environment.
Individual CAAStools jobs request 1 CPU, 2 GB and 10 minutes. A timeout or
other worker failure is ignored rather than terminating the workflow. Merge
jobs request 1 CPU, 4 GB and 2 hours.

## Merge and significance filter

Each successful gene/config task emits a headered result, including a
header-only result when no CAAS is found. A timed-out or failed job has no
result and is omitted from the merge, but does not terminate the workflow.
After every discovery task has completed or been ignored, Nextflow:

1. merges all genes separately for each configuration;
2. merges the four per-configuration tables into one table;
3. retains exactly one header;
4. writes a second table containing rows with `Pvalue <= 0.05`.

Final files are written to `results/RUN_ID/merged/`:

- `caas.all-configs.all-genes.tsv`;
- `caas.all-configs.all-genes.significant.tsv`.

The intermediate per-configuration merges are written under
`results/RUN_ID/merged/by-config/`, while individual gene results are under
`results/RUN_ID/caas/by-config/CONFIG_ID/`.
