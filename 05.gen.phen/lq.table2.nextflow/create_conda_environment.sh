#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
conda_init_script="${GENPHEN_CONDA_SH:-/homes/aplic/noarch/software/Miniconda3/23.9.0-0/etc/profile.d/conda.sh}"
conda_environment="${GENPHEN_CONDA_ENV:-phyloq}"

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

Rscript -e '
  if (!requireNamespace("RERconverge", quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) {
      stop("The remotes package is required to install RERconverge.")
    }
    remotes::install_github(
      "nclark-lab/RERconverge",
      ref = "master",
      dependencies = TRUE,
      upgrade = "never"
    )
  }
  suppressPackageStartupMessages(library(RERconverge))
  cat("RERconverge ", as.character(packageVersion("RERconverge")), " is ready.\n", sep = "")
'

python3 -c 'import Bio, numpy, scipy; print("Python dependencies: OK")'
Rscript -e 'suppressPackageStartupMessages({library(ape); library(phangorn); library(RERconverge)}); cat("R dependencies: OK\n")'
nextflow -version

echo "Conda environment ${conda_environment} is ready."
