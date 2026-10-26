#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from datetime import datetime

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Merge representative TPM, MMseqs2 taxonomy and blastx results into one table"
    )
    p.add_argument("--tpm", required=True, help="TPM table with rep_contig and one column per sample")
    p.add_argument("--taxonomy", required=True, help="clusters_with_taxonomy.tsv (contig_id, cluster_id, MMseqs2 taxonomy)")
    p.add_argument("--blastx", required=True, help="blastx table from fill_missing_blastx.py (qseqid = <contig>_cluster_<n>)")
    p.add_argument("--output", required=True)
    return p.parse_args()


def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)


def main():
    args = parse_args()

    tpm = pd.read_csv(args.tpm, sep="\t")
    taxonomy = pd.read_csv(args.taxonomy, sep="\t")
    blastx = pd.read_csv(args.blastx, sep="\t", dtype=str)
    log(f"Representatives in TPM table: {len(tpm)}")

    tpm = tpm.rename(columns={"rep_contig": "contig_id"})

    # Cluster size is counted over all contigs, before taxonomy is reduced to representatives
    cluster_size = taxonomy.groupby("cluster_id").size().rename("cluster_size")
    tpm.insert(tpm.columns.get_loc("cluster_id") + 1, "cluster_size", tpm["cluster_id"].map(cluster_size))

    # MMseqs2 taxonomy of the representative contig itself
    taxonomy = taxonomy.drop(columns=["cluster_id"])
    merged = tpm.merge(taxonomy, on="contig_id", how="left", validate="one_to_one")

    # Representative FASTA headers (and so blastx qseqid) are <contig_id>_cluster_<n>
    blastx["contig_id"] = blastx["qseqid"].str.rsplit("_cluster_", n=1).str[0]
    blastx = blastx.drop(columns=["qseqid"])
    blastx = blastx.rename(columns={c: f"blastx_{c}" for c in blastx.columns
                                    if c not in ("contig_id", "blastx_status")})
    merged = merged.merge(blastx, on="contig_id", how="left", validate="one_to_one")

    log(f"With MMseqs2 taxonomy: {merged['taxonomic_name'].notna().sum()}")
    log(f"With blastx hit: {(merged['blastx_status'] == 'hit').sum()}")
    if merged["blastx_status"].isna().any():
        log(f"[WARN] Missing from blastx table: {merged['blastx_status'].isna().sum()}")

    merged.to_csv(args.output, sep="\t", index=False)
    log(f"Wrote {len(merged)} rows to {args.output}")


if __name__ == "__main__":
    main()
