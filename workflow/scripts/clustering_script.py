#!/usr/bin/env python3

import argparse
import pandas as pd
from Bio import SeqIO
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
parser = argparse.ArgumentParser(description="Filter BLAST hits and extract clusters.")
parser.add_argument("--fasta", required=True)
parser.add_argument("--blast", required=True)
parser.add_argument("--edges", required=True)
parser.add_argument("--filtered", required=True)
parser.add_argument("--clusters", required=True, help="Output cluster list TSV (contig -> cluster ID).")
parser.add_argument("--min_identity", type=float, default=90.0)
parser.add_argument("--min_length", type=int, default=200)
args = parser.parse_args()

# -----------------------------------------
# Load BLAST table
# -----------------------------------------
logging.info("Loading BLAST file...")
df = pd.read_csv(args.blast, sep="\t", header=None)
df.columns = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send"]
logging.info(f"Loaded {len(df)} total BLAST hits.")

# -----------------------------------------
# Remove self hits
# -----------------------------------------
df = df[df["qseqid"] != df["sseqid"]]
logging.info(f"{len(df)} hits remaining after removing self-hits.")

# -----------------------------------------
# Deduplicate hits
# -----------------------------------------
df["pair_key"] = df.apply(lambda row: tuple(sorted([row["qseqid"], row["sseqid"]])), axis=1)
df = df.drop_duplicates(subset="pair_key").drop(columns=["pair_key"])
logging.info(f"{len(df)} hits remaining after deduplication.")

# -----------------------------------------
# Apply thresholds
# -----------------------------------------
filtered_df = df[(df["pident"] >= args.min_identity) & (df["length"] > args.min_length)]
logging.info(f"{len(filtered_df)} hits remaining after filtering at "
             f"{args.min_identity}% identity and {args.min_length} bp.")

# Write filtered BLAST table
filtered_df.to_csv(args.filtered, sep="\t", index=False)
logging.info(f"Filtered BLAST table written to {args.filtered}")

# -----------------------------------------
# Write edge list
# -----------------------------------------
with open(args.edges, "w") as out:
    for _, row in filtered_df.iterrows():
        out.write(f"{row['qseqid']}\t{row['sseqid']}\n")
logging.info(f"Edge list written to {args.edges}")

# -----------------------------------------
# Extract clusters (connected components)
# -----------------------------------------
import networkx as nx

G = nx.Graph()
all_contigs = [record.id for record in SeqIO.parse(args.fasta, "fasta")]
G.add_nodes_from(all_contigs)

with open(args.edges) as f:
    for line in f:
        q, s = line.strip().split("\t")
        G.add_edge(q, s)
logging.info(f"Graph constructed with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")

clusters = list(nx.connected_components(G))
logging.info(f"Identified {len(clusters)} clusters.")

# Write cluster list
with open(args.clusters, "w") as out:
    out.write("contig_id\tcluster_id\n")
    for i, cluster in enumerate(clusters, 1):
        cluster_id = f"cluster_{i}"
        for node in cluster:
            out.write(f"{node}\t{cluster_id}\n")
logging.info(f"Cluster list written to {args.clusters}")
