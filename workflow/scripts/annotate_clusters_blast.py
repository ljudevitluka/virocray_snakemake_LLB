#!/usr/bin/env python3

import argparse
import pandas as pd
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s:%(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

parser = argparse.ArgumentParser(description="Extract intra-cluster BLAST pairs with sample IDs.")
parser.add_argument("--blast", required=True, help="Filtered BLAST pairs (TSV).")
parser.add_argument("--clusters", required=True, help="Contig -> cluster mapping (TSV).")
parser.add_argument("--out", required=True, help="Output TSV with intra-cluster BLAST pairs.")
args = parser.parse_args()

# Load blast results
blast_df = pd.read_csv(args.blast, sep="\t")
logging.info(f"Loaded {len(blast_df)} BLAST pairs")

# Load clusters
clusters_df = pd.read_csv(args.clusters, sep="\t")
cluster_map = dict(zip(clusters_df["contig_id"], clusters_df["cluster_id"]))
logging.info(f"Loaded {len(cluster_map)} contig->cluster mappings")

# Annotate clusters
blast_df["q_cluster"] = blast_df["qseqid"].map(cluster_map)
blast_df["s_cluster"] = blast_df["sseqid"].map(cluster_map)

# Keep only intra-cluster pairs
intra_df = blast_df[blast_df["q_cluster"] == blast_df["s_cluster"]].copy()
logging.info(f"Found {len(intra_df)} intra-cluster pairs")

# Extract sample IDs (everything before "_NODE")
intra_df["q_id"] = intra_df["qseqid"].str.extract(r"^(.*?_E)_NODE")[0]
intra_df["s_id"] = intra_df["sseqid"].str.extract(r"^(.*?_E)_NODE")[0]

# Reorder and rename columns
final_df = intra_df[[
    "qseqid", "q_id",
    "sseqid", "s_id",
    "pident", "length", "qstart", "qend", "sstart", "send",
    "q_cluster"
]].rename(columns={"q_cluster": "cluster_id"})

# Save
final_df.to_csv(args.out, sep="\t", index=False)
logging.info(f"Intra-cluster BLAST pairs written to {args.out}")
