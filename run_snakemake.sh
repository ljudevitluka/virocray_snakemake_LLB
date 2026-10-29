#!/usr/bin/env bash
set -euo pipefail

CONDA_EXE="${CONDA_EXE:-/proj/etools/conda/bin/conda}"
# Memory budget for the jobs (GB). Snakemake only starts jobs whose
# resources mem_mb fit in it.
MAX_MEM_GB="${MAX_MEM_GB:-500}"
# Hard limit for the whole run (GB). It also counts cached database files
# (the MMseqs2 nr database alone is 414 GB), so it is higher than MAX_MEM_GB.
HARD_MEM_GB="${HARD_MEM_GB:-900}"

cmd=("$CONDA_EXE" run --no-capture-output -n snakemake
     snakemake --use-conda --conda-frontend conda
     --resources mem_mb=$((MAX_MEM_GB * 1000)) "$@")

if systemd-run --user --scope -q true 2>/dev/null; then
    exec systemd-run --user --scope -q \
        -p MemoryMax="${HARD_MEM_GB}G" -p MemorySwapMax=0 "${cmd[@]}"
else
    echo "WARNING: systemd-run --user is not available; running without a hard memory limit" >&2
    exec "${cmd[@]}"
fi
