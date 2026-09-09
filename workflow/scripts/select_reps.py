#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def parse_args():
    p = argparse.ArgumentParser(description="Select longest representative contig per cluster and write FASTA + TSV + stats")
    p.add_argument("--fasta", required=True, help="FASTA file with all viral contigs")
    p.add_argument("--clusters", required=True, help="TSV file with columns: contig_id, cluster_id")
    p.add_argument("--reps-fasta", required=True, help="Output FASTA file with representative contigs")
    p.add_argument("--reps-tsv", required=True, help="Output TSV file with representative metadata")
    p.add_argument("--stats", required=True, help="Output TSV/TXT file with cluster statistics")
    return p.parse_args()

def log(msg):
    print(f"[INFO] {msg}", flush=True)

def main():
    args = parse_args()
    log(f"Loading clusters table from {args.clusters}")

    clusters_df = pd.read_csv(args.clusters, sep="\t")
    if not {"contig_id", "cluster_id"}.issubset(clusters_df.columns):
        raise ValueError("Clusters file must contain columns: contig_id, cluster_id")
    log(f"Loaded {len(clusters_df)} contig cluster assignments")

    cluster_sizes = clusters_df.groupby("cluster_id").size().rename("cluster_size").reset_index()

    log(f"Loading contigs from FASTA {args.fasta}")
    seq_dict = {rec.id: rec for rec in SeqIO.parse(args.fasta, "fasta")}
    log(f"Loaded {len(seq_dict)} contigs from FASTA")

    missing = set(clusters_df["contig_id"]) - set(seq_dict.keys())
    if missing:
        log(f"WARNING: {len(missing)} contigs from clusters file not found in FASTA")

    clusters_df["length"] = clusters_df["contig_id"].map(
        lambda cid: len(seq_dict[cid].seq) if cid in seq_dict else 0
    )
    before = len(clusters_df)
    clusters_df = clusters_df[clusters_df["length"] > 0].copy()
    after = len(clusters_df)
    if after < before:
        log(f"WARNING: Removed {before - after} contigs with length 0")

    reps_df = clusters_df.sort_values("length", ascending=False).groupby("cluster_id", as_index=False).first()
    reps_df = reps_df[["cluster_id", "contig_id", "length"]].rename(columns={
        "contig_id": "rep_contig",
        "length": "rep_contig_length"
    })

    reps_df = reps_df.merge(cluster_sizes, on="cluster_id", how="left")
    reps_df["singleton"] = (reps_df["cluster_size"] == 1).astype(int)
    reps_df = reps_df.drop(columns=["cluster_size"])

    log(f"Selected {len(reps_df)} representative contigs (one per cluster)")

    log(f"Writing representatives TSV to {args.reps_tsv}")
    reps_df.to_csv(args.reps_tsv, sep="\t", index=False)

    log(f"Writing representatives FASTA to {args.reps_fasta}")
    with open(args.reps_fasta, "w") as out_fa:
        for _, row in reps_df.iterrows():
            rec = seq_dict[row["rep_contig"]]
            new_id = f"{row['rep_contig']}_{row['cluster_id']}"
            new_rec = SeqRecord(rec.seq, id=new_id, name=new_id, description="")
            SeqIO.write(new_rec, out_fa, "fasta")

    total_clusters = len(reps_df)
    n_singletons = int((reps_df["singleton"] == 1).sum())
    n_non_singletons = total_clusters - n_singletons

    log(f"Writing cluster statistics to {args.stats}")
    with open(args.stats, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"total_clusters\t{total_clusters}\n")
        fh.write(f"singletons\t{n_singletons}\n")
        fh.write(f"non_singletons\t{n_non_singletons}\n")

    log("Done.")

if __name__ == "__main__":
    main()
