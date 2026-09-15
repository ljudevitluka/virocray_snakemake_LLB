#!/usr/bin/env bash
set -euo pipefail

CONDA_EXE="${CONDA_EXE:-/proj/etools/conda/bin/conda}"

exec "$CONDA_EXE" run --no-capture-output -n snakemake \
    snakemake --use-conda "$@"
