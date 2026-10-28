#!/usr/bin/env bash
set -euo pipefail

CONDA_EXE="${CONDA_EXE:-/proj/etools/conda/bin/conda}"
# Memory budget for the whole run (GB). Snakemake only starts jobs whose
# resources mem_mb fit in it, and the run is placed in a cgroup with a hard
# MemoryMax, so a job that exceeds it is killed alone instead of the server
# running out of memory (which on 2026-10-27 killed all screen sessions).
MAX_MEM_GB="${MAX_MEM_GB:-500}"

cmd=("$CONDA_EXE" run --no-capture-output -n snakemake
     snakemake --use-conda --conda-frontend conda
     --resources mem_mb=$((MAX_MEM_GB * 1000)) "$@")

if systemd-run --user --scope -q true 2>/dev/null; then
    exec systemd-run --user --scope -q \
        -p MemoryMax="${MAX_MEM_GB}G" -p MemorySwapMax=0 "${cmd[@]}"
else
    echo "WARNING: systemd-run --user is not available; running without a hard memory limit" >&2
    exec "${cmd[@]}"
fi
