#!/usr/bin/env python3

import argparse
import logging
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s:%(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

parser = argparse.ArgumentParser(
    description="Select longest representative contig per cluster and write FASTA + TSV."
)
parser.add_argument("--fasta", required=True,
                    help="FASTA file with all viral contigs (e.g. viral_contigs.fasta).")
parser.add_argument("--clusters", required=True,
                    help="TSV file with columns: contig_id, cluster_id.")
parser.add_argument("--reps-fasta", required=True,
                    help="Output FASTA file with representative contigs.")
parser.add_argument("--reps-tsv", required=True,
                    help="Output TSV file with cluster_id, rep_contig, rep_contig_length.")
args = parser.parse_args()

# Load clusters table
clusters_df = pd.read_csv(args.clusters, sep="\t")
if not {"contig_id", "cluster_id"}.issubset(clusters_df.columns):
    raise ValueError("Clusters file must contain columns: contig_id, cluster_id")

logging.info(f"Loaded {len(clusters_df)} contig->cluster assignments")

# Load all viral contigs into a dict
seq_dict = {}
for rec in SeqIO.parse(args.fasta, "fasta"):
    seq_dict[rec.id] = rec
logging.info(f"Loaded {len(seq_dict)} contigs from FASTA")

# Check which contigs from clusters are missing in FASTA
missing = set(clusters_df["contig_id"]) - set(seq_dict.keys())
if missing:
    logging.warning(f"{len(missing)} contigs from clusters file not found in FASTA")

# Add length column (only for contigs present in FASTA)
clusters_df["length"] = clusters_df["contig_id"].map(
    lambda cid: len(seq_dict[cid].seq) if cid in seq_dict else 0
)

# Remove contigs with length 0 (not present in FASTA)
before = len(clusters_df)
clusters_df = clusters_df[clusters_df["length"] > 0].copy()
after = len(clusters_df)
if after < before:
    logging.warning(f"Removed {before - after} contigs with length 0 (missing in FASTA)")

# For each cluster, pick longest contig as representative
reps_df = (
    clusters_df
    .sort_values("length", ascending=False)
    .groupby("cluster_id", as_index=False)
    .first()
)

reps_df = reps_df[["cluster_id", "contig_id", "length"]]
reps_df = reps_df.rename(columns={
    "contig_id": "rep_contig",
    "length": "rep_contig_length"
})

logging.info(f"Selected {len(reps_df)} representative contigs (one per cluster)")

# Write TSV
reps_df.to_csv(args.reps_tsv, sep="\t", index=False)
logging.info(f"Wrote representatives table to {args.reps_tsv}")

# Write FASTA: header = original contig ID + "_cluster_X"
with open(args.reps_fasta, "w") as out_fa:
    for _, row in reps_df.iterrows():
        cluster_id = row["cluster_id"]
        contig_id = row["rep_contig"]
        rec = seq_dict[contig_id]
        new_id = f"{contig_id}_{cluster_id}"
        # Create a new SeqRecord to avoid mutating the original
        new_rec = SeqRecord(
            rec.seq,
            id=new_id,
            name=new_id,
            description=""
        )
        SeqIO.write(new_rec, out_fa, "fasta")

logging.info(f"Wrote representative contigs FASTA to {args.reps_fasta}")
