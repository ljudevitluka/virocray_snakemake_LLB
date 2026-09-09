#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(description="Merge cluster assignments with merged MMseqs2 taxonomy data.")
    p.add_argument("--clusters", required=True, help="TSV file with contig_id and cluster_id.")
    p.add_argument("--taxonomy", required=True, help="Merged MMseqs2 taxonomy file.")
    p.add_argument("--output", required=True, help="Output merged TSV file.")
    return p.parse_args()

def log(msg):
    print(f"[INFO] {msg}", flush=True)

def main():
    args = parse_args()

    log("Loading cluster file...")
    clusters_df = pd.read_csv(args.clusters, sep="\t")
    if not {"contig_id", "cluster_id"}.issubset(clusters_df.columns):
        raise ValueError("clusters file must contain columns: contig_id, cluster_id")
    log(f"Loaded {len(clusters_df)} cluster assignments.")

    cluster_sizes = clusters_df.groupby("cluster_id")["contig_id"].size().rename("cluster_size")
    clusters_df = clusters_df.merge(cluster_sizes, on="cluster_id", how="left")
    clusters_df["singleton"] = (clusters_df["cluster_size"] == 1).astype(int)

    log(f"Reading taxonomy file: {args.taxonomy}")
    taxonomy_df = pd.read_csv(args.taxonomy, sep="\t", header=None, dtype=str)

    if taxonomy_df.shape[1] < 9:
        raise ValueError(
            f"Taxonomy file must contain at least 9 columns, but {taxonomy_df.shape[1]} were found."
        )

    if taxonomy_df.shape[1] > 9:
        log(f"Taxonomy file has {taxonomy_df.shape[1]} columns; keeping only the first 9 MMseqs2 columns.")
        taxonomy_df = taxonomy_df.iloc[:, :9].copy()

    taxonomy_df.columns = [
        "contig_id", "taxonomic_identifier", "taxonomic_rank", "taxonomic_name",
        "no_of_retained_fragments", "no_of_assigned_fragments",
        "in_agreement_with_contig_label", "support", "lineage"
    ]
    log(f"Total taxonomy entries loaded: {len(taxonomy_df)}")

    merged_df = clusters_df.merge(taxonomy_df, on="contig_id", how="left")
    log(f"Merged table has {len(merged_df)} rows.")

    merged_df = merged_df.drop(columns=["cluster_size"])
    merged_df.to_csv(args.output, sep="\t", index=False)
    log(f"Merged table written to {args.output}")
    log("Merge complete.")

if __name__ == "__main__":
    main()