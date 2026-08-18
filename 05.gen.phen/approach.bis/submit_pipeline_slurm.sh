#!/usr/bin/env bash

#SBATCH --job-name=genphen-caas-bis-driver
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=7-00:00:00
#SBATCH --partition=std-cpu
#SBATCH --output=logs/genphen-caas-bis-driver-%j.out
#SBATCH --error=logs/genphen-caas-bis-driver-%j.err

##SBATCH --account=YOUR_ACCOUNT
##SBATCH --qos=YOUR_QOS
##SBATCH --constraint=YOUR_CONSTRAINT

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "This is a SLURM batch script. Submit it with:" >&2
    echo "  sbatch submit_pipeline_slurm.sh" >&2
    exit 1
fi

pipeline_dir="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
if [[ ! -f "${pipeline_dir}/run_pipeline.sh" || ! -f "${pipeline_dir}/main.nf" ]]; then
    echo "Submit this script from the approach.bis directory." >&2
    echo "SLURM submission directory: ${pipeline_dir}" >&2
    exit 1
fi
cd "${pipeline_dir}"

export GENPHEN_NEXTFLOW_CONFIG="${GENPHEN_NEXTFLOW_CONFIG:-${pipeline_dir}/conf/cluster.config}"

conda_init_script="${GENPHEN_CONDA_SH:-/homes/aplic/noarch/software/Miniconda3/23.9.0-0/etc/profile.d/conda.sh}"
conda_environment="${GENPHEN_CONDA_ENV:-genphen-caas-bis}"
if [[ ! -r "${conda_init_script}" ]]; then
    echo "Cannot read Conda initialization script: ${conda_init_script}" >&2
    exit 1
fi
source "${conda_init_script}"
conda activate "${conda_environment}"

export NXF_OPTS="${NXF_OPTS:--Xms512m -Xmx8g}"

echo "GenPhen CAAStools BIS Nextflow driver"
echo "SLURM job: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Started: $(date --iso-8601=seconds)"
echo "Nextflow config: ${GENPHEN_NEXTFLOW_CONFIG}"
echo "Conda environment: ${CONDA_DEFAULT_ENV:-${conda_environment}}"
echo "Nextflow JVM options: ${NXF_OPTS}"

exec bash run_pipeline.sh "$@"
