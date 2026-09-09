#!/usr/bin/env python3

import argparse
import pandas as pd
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# -------------------------
# Argument parser
# -------------------------
parser = argparse.ArgumentParser(
    description="Merge per-sample CoverM TSV outputs into three matrices: counts, TPM, RPKM."
)
parser.add_argument(
    "--inputs", nargs="+", required=True,
    help="List of per-sample CoverM TSV files (one per sample)."
)
parser.add_argument(
    "--counts", required=True,
    help="Output TSV file for counts matrix."
)
parser.add_argument(
    "--tpm", required=True,
    help="Output TSV file for TPM matrix."
)
parser.add_argument(
    "--rpkm", required=True,
    help="Output TSV file for RPKM matrix."
)
args = parser.parse_args()

# -------------------------
# Load and merge
# -------------------------
dfs_counts = []
dfs_tpm = []
dfs_rpkm = []

logging.info(f"Processing {len(args.inputs)} CoverM TSV files...")

for f in args.inputs:
    sample_id = f.split("/")[-1].split("_coverm.tsv")[0].replace("05_", "")
    logging.info(f"Reading {f} for sample {sample_id}...")

    df = pd.read_csv(f, sep="\t")
    if "Contig" not in df.columns:
        raise ValueError(f"File {f} does not contain 'Contig' column")

    df = df.set_index("Contig")

    # Identify correct columns by searching
    count_col = [c for c in df.columns if "Read Count" in c][0]
    tpm_col   = [c for c in df.columns if "TPM" in c][0]
    rpkm_col  = [c for c in df.columns if "RPKM" in c][0]

    dfs_counts.append(df[count_col].rename(sample_id))
    dfs_tpm.append(df[tpm_col].rename(sample_id))
    dfs_rpkm.append(df[rpkm_col].rename(sample_id))

# -------------------------
# Concatenate and write
# -------------------------
logging.info("Merging counts...")
pd.concat(dfs_counts, axis=1).to_csv(args.counts, sep="\t")

logging.info("Merging TPM...")
pd.concat(dfs_tpm, axis=1).to_csv(args.tpm, sep="\t")

logging.info("Merging RPKM...")
pd.concat(dfs_rpkm, axis=1).to_csv(args.rpkm, sep="\t")

logging.info("Done! Outputs written:")
logging.info(f" - Counts: {args.counts}")
logging.info(f" - TPM:    {args.tpm}")
logging.info(f" - RPKM:   {args.rpkm}")
