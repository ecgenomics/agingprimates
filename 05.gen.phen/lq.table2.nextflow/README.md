# LQ table 2 Nextflow workflow

SLURM-ready workflow for the comparison configurations in
`../lq.table2.4species.comparisons`.

It follows the architecture of `traits.evolution.omar/phyloq/bootstrap` while
remaining self-contained: the CAAStools program and the RERconverge helper
scripts are copied into `bin/`.

## What it runs

- **CAAStools:** every alignment is combined with each of the 100 unique
  four-FG versus four-BG configurations.
- **RERconverge:** each of the 35 unique four-FG configurations is analysed
  against the complete fixed six-species BG set.

The two branches are independent and can be enabled or disabled in
`conf/cluster.config`.

## Configure the cluster run

1. Create or update the Conda environment once on the cluster login node:

   ```bash
   bash create_conda_environment.sh
   ```

   The script uses Correfoc's Miniconda installation and creates the `phyloq`
   environment from `environment.yml`. It also installs RERconverge from its
   official GitHub repository when it is missing.
2. Edit `conf/cluster.config`:
   - replace all `/ABSOLUTE/CLUSTER/PATH/...` entries;
   - set account, partition, QoS or constraint when required;
   - change `conda_init_script` or `conda_environment` if the cluster setup is
     different;
   - adjust resources and analytical parameters.
3. Edit the `#SBATCH` header in `submit_pipeline_slurm.sh` for the small
   Nextflow driver job. Optional account, partition, QoS and constraint lines
   are provided as inactive `##SBATCH` examples.
4. Submit the driver job. It initializes Miniconda and activates `phyloq`
   before invoking Nextflow; each analysis task repeats the same activation.

CAAStools requires Python with Biopython, SciPy and NumPy. RERconverge requires
R with `ape` and `RERconverge`. Gene-tree construction additionally requires
`phangorn`.

## Launch

From this directory on a SLURM login node:

```bash
sbatch submit_pipeline_slurm.sh
```

Resume the latest run:

```bash
sbatch submit_pipeline_slurm.sh -resume
```

For an interactive launch using the same Nextflow configuration:

```bash
bash run_pipeline.sh
```

Additional Nextflow arguments can be appended to either command.

## Outputs

The output root is controlled by `params.results_root` in
`conf/cluster.config`. Each launch receives a timestamped `run_id` and writes:

- `RUN_ID/caas/CONFIG_ID/`;
- `RUN_ID/rerconverge/CONFIG_ID/`.

Reports, traces and timelines are written to `logs/`. SLURM driver logs are
written as `genphen-driver-JOB_ID.out` and `genphen-driver-JOB_ID.err`.

## Bundled executables

- `bin/caastools/ct` and its Python modules;
- `bin/rerconverge/run_rerconverge.R`;
- `bin/rerconverge/build_gene_trees.R`.

The bundled files are snapshots copied from the reference workflow. Their
language and package dependencies are defined in `environment.yml`.
