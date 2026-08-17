#!/usr/bin/env bash

#SBATCH --job-name=genphen-driver
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00
#SBATCH --partition=std-cpu
#SBATCH --output=genphen-driver-%j.out
#SBATCH --error=genphen-driver-%j.err

# Uncomment and edit the directives required by the cluster.
##SBATCH --account=YOUR_ACCOUNT
##SBATCH --qos=YOUR_QOS
##SBATCH --constraint=YOUR_CONSTRAINT

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${script_dir}"

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "This is a SLURM batch script. Submit it with:" >&2
    echo "  sbatch submit_pipeline_slurm.sh" >&2
    exit 1
fi

export GENPHEN_NEXTFLOW_CONFIG="${GENPHEN_NEXTFLOW_CONFIG:-${script_dir}/conf/cluster.config}"

echo "GenPhen Nextflow driver"
echo "SLURM job: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Started: $(date --iso-8601=seconds)"
echo "Nextflow config: ${GENPHEN_NEXTFLOW_CONFIG}"

# If Nextflow/Java are provided as modules, load them here. For example:
# source /etc/profile.d/modules.sh
# module load nextflow

exec bash run_pipeline.sh "$@"
