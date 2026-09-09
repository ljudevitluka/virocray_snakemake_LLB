#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description="Extract contigs labeled as unclassified from merged FASTA.")
parser.add_argument("--fasta", required=True, help="Input merged FASTA with renamed contigs")
parser.add_argument("--taxonomy", required=True, help="Merged mmseqs2 taxonomy TSV")
parser.add_argument("--output", required=True, help="Output FASTA for unclassified contigs")
args = parser.parse_args()

# Step 1: Load taxonomy file
df = pd.read_csv(args.taxonomy, sep="\t", header=None, dtype=str)

# Step 2: Filter rows where rank == "no rank" AND name == "unclassified" (case-insensitive)
df_unclassified = df[
    (df[2].str.lower() == "no rank") &
    (df[3].str.lower() == "unclassified")
].copy()

df_unclassified.columns = ["contig", "tax_id", "rank", "name", "retained", "assigned", "label_match", "support"] + ([ "lineage" ] if df.shape[1] == 9 else [])

print(f"Total unclassified contigs: {len(df_unclassified)}")

# Step 3: Extract matching contigs from FASTA
unclassified_ids = set(df_unclassified["contig"])

count_written = 0
with open(args.output, "w") as out_fasta:
    for record in SeqIO.parse(args.fasta, "fasta"):
        if record.id in unclassified_ids:
            SeqIO.write(record, out_fasta, "fasta")
            count_written += 1

print(f"Total unclassified contigs written: {count_written}")
