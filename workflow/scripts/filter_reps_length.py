#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
from datetime import datetime
import pandas as pd
from Bio import SeqIO


def parse_args():
    p = argparse.ArgumentParser(description="Filter representative contigs by minimum length")
    p.add_argument("--reps-fasta", required=True)
    p.add_argument("--reps-tsv", required=True)
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--stats", required=True)
    p.add_argument("--min-len", type=int, required=True)
    return p.parse_args()


def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)


def main():
    args = parse_args()
    log(f"Minimum length cutoff: {args.min_len} nt")

    reps_df = pd.read_csv(args.reps_tsv, sep="\t")

    # ----------------------------
    # REMOVE UNNECESSARY COLUMN
    # ----------------------------
    if "contig_id" in reps_df.columns:
        reps_df = reps_df.drop(columns=["contig_id"])

    id_col = "rep_contig"
    len_col = "rep_contig_length"
    singleton_col = "singleton"
    has_singleton = singleton_col in reps_df.columns

    missing_cols = {id_col, len_col} - set(reps_df.columns)
    if missing_cols:
        sys.exit(f"ERROR: Missing required column(s) in TSV: {', '.join(missing_cols)}")

    keep_df = reps_df[reps_df[len_col] >= args.min_len].copy()
    keep_ids = set(keep_df[id_col].astype(str))

    log(f"Representatives kept after length filter: {len(keep_df)}")

    keep_df.to_csv(args.out_tsv, sep="\t", index=False)

    with open(args.out_fasta, "w") as out_fa:
        for rec in SeqIO.parse(args.reps_fasta, "fasta"):
            rec_id_base = rec.id.split("_cluster_")[0]
            if rec_id_base in keep_ids:
                SeqIO.write(rec, out_fa, "fasta")

    total_reps = len(keep_df)
    n_singletons = int((keep_df[singleton_col] == 1).sum()) if has_singleton else 0
    n_non_singletons = total_reps - n_singletons

    with open(args.stats, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"total_representatives\t{total_reps}\n")
        fh.write(f"singletons\t{n_singletons}\n")
        fh.write(f"non_singletons\t{n_non_singletons}\n")

    log("Done.")


if __name__ == "__main__":
    main()