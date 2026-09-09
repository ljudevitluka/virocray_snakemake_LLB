#!/usr/bin/env python3

import argparse
import pandas as pd
from Bio import SeqIO
import logging

# -----------------------------------------
# Setup logging
# -----------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# -----------------------------------------
# Parse command line arguments
# -----------------------------------------
parser = argparse.ArgumentParser(description="Extract clusters containing viral contigs and add sequences.")
parser.add_argument("--input", required=True, help="Input TSV with merged clustering + taxonomy.")
parser.add_argument("--fasta", required=True, help="FASTA file with contig sequences.")
parser.add_argument("--output", required=True, help="Output TSV with sequences for viral clusters.")
args = parser.parse_args()

# -----------------------------------------
# Load merged TSV
# -----------------------------------------
logging.info("Loading merged clustering + taxonomy file...")
df = pd.read_csv(args.input, sep="\t")
logging.info(f"Loaded {len(df)} total contigs.")

# -----------------------------------------
# Find clusters that contain at least one viral contig
# -----------------------------------------
viral_contigs = df[df["lineage"].str.startswith("-_Viruses;", na=False)]
viral_clusters = viral_contigs["cluster_id"].unique()
logging.info(f"Identified {len(viral_clusters)} clusters containing at least one viral contig.")

# Filter to all contigs in these clusters
filtered_df = df[df["cluster_id"].isin(viral_clusters)]
logging.info(f"Filtered dataset has {len(filtered_df)} contigs across these clusters.")

# -----------------------------------------
# Load sequences from FASTA
# -----------------------------------------
logging.info("Loading sequences from FASTA...")
seq_dict = {record.id: str(record.seq) for record in SeqIO.parse(args.fasta, "fasta")}
logging.info(f"Loaded {len(seq_dict)} total sequences.")

# Add sequence column
filtered_df["sequence"] = filtered_df["contig_id"].map(seq_dict)
missing_seqs = filtered_df["sequence"].isna().sum()
if missing_seqs > 0:
    logging.warning(f"{missing_seqs} sequences could not be found in FASTA.")

# -----------------------------------------
# Write output
# -----------------------------------------
filtered_df.to_csv(args.output, sep="\t", index=False)
logging.info(f"Filtered table with sequences written to {args.output}")
logging.info("Extraction complete.")
