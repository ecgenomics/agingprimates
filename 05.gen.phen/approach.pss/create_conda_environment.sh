#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
conda_init_script="${GENPHEN_CONDA_SH:-/homes/aplic/noarch/software/Miniconda3/23.9.0-0/etc/profile.d/conda.sh}"
conda_environment="${GENPHEN_CONDA_ENV:-genphen-caas-pss}"

if [[ ! -r "${conda_init_script}" ]]; then
    echo "Cannot read Conda initialization script: ${conda_init_script}" >&2
    exit 1
fi

source "${conda_init_script}"

if conda run -n "${conda_environment}" true >/dev/null 2>&1; then
    echo "Updating Conda environment: ${conda_environment}"
    conda env update \
        --name "${conda_environment}" \
        --file "${script_dir}/environment.yml" \
        --prune
else
    echo "Creating Conda environment: ${conda_environment}"
    conda env create \
        --name "${conda_environment}" \
        --file "${script_dir}/environment.yml"
fi

conda activate "${conda_environment}"

python3 -c 'import Bio, numpy, scipy; print("CAAStools Python dependencies: OK")'
python3 "${script_dir}/../lq.table2.nextflow/bin/caastools/ct" help >/dev/null
nextflow -version

echo "Conda environment ${conda_environment} is ready."
