#!/usr/bin/env python3

import argparse
import pandas as pd
import logging

# -----------------------------------------
# Set up logging
# -----------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# -----------------------------------------
# Parse command line arguments
# -----------------------------------------
parser = argparse.ArgumentParser(description="Split clusters into singletons vs non-singletons.")
parser.add_argument("--input", required=True, help="Input TSV file with contig->cluster and taxonomy info.")
parser.add_argument("--wo_singletons", required=True, help="Output TSV with clusters that have >1 contig.")
parser.add_argument("--singletons", required=True, help="Output TSV with singleton clusters (exactly 1 contig).")
args = parser.parse_args()

# -----------------------------------------
# Load input table
# -----------------------------------------
logging.info(f"Loading {args.input} ...")
df = pd.read_csv(args.input, sep="\t")

if "cluster_id" not in df.columns:
    raise ValueError("Input file must contain a 'cluster_id' column.")

logging.info(f"Loaded {len(df)} rows from {args.input}")

# -----------------------------------------
# Count contigs per cluster
# -----------------------------------------
cluster_counts = df["cluster_id"].value_counts()
logging.info(f"Found {len(cluster_counts)} unique clusters")

# Map counts back to dataframe
df["cluster_size"] = df["cluster_id"].map(cluster_counts)

# -----------------------------------------
# Split into singletons vs non-singletons
# -----------------------------------------
non_singletons_df = df[df["cluster_size"] > 1].drop(columns=["cluster_size"])
singletons_df = df[df["cluster_size"] == 1].drop(columns=["cluster_size"])

logging.info(f"Non-singletons: {len(non_singletons_df)} rows from {sum(cluster_counts>1)} clusters")
logging.info(f"Singletons: {len(singletons_df)} rows from {sum(cluster_counts==1)} clusters")

# -----------------------------------------
# Write outputs
# -----------------------------------------
non_singletons_df.to_csv(args.wo_singletons, sep="\t", index=False)
singletons_df.to_csv(args.singletons, sep="\t", index=False)

logging.info(f"Written non-singletons to {args.wo_singletons}")
logging.info(f"Written singletons to {args.singletons}")
