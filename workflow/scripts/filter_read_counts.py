#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from datetime import datetime
import argparse


def parse_args():
    p = argparse.ArgumentParser(description="Filter vOTUs by min reads and breadth")
    p.add_argument("--samples", required=True)
    p.add_argument("--breadth", required=True)
    p.add_argument("--out-filtered", required=True)
    p.add_argument("--out-stats", required=True)
    p.add_argument("--min-reads", type=int, default=2)
    p.add_argument("--min-breadth", type=float, default=0.0)
    return p.parse_args()


def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)


def get_readcount_cols(df):
    return [c for c in df.columns if c not in [
        "cluster_id",
        "rep_contig",
        "rep_contig_length",
        "singleton"
    ]]


def ensure_int(df, cols):
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    return df


def apply_read_breadth_filter(counts_df, breadth_df, min_reads, min_breadth):
    out = counts_df.copy()

    for sample in out.columns:
        if sample not in breadth_df.columns:
            continue

        breadth_vals = out.index.map(breadth_df[sample]).fillna(0)

        low_reads = out[sample] < min_reads
        low_breadth = breadth_vals < min_breadth

        out.loc[low_reads | low_breadth, sample] = 0

    return out


def count_vOTU_detections(df):
    return (df > 0).sum().sum()


def main():
    args = parse_args()

    log("Reading samples table")
    samples = pd.read_csv(args.samples, sep="\t").set_index("rep_contig")

    sample_cols = get_readcount_cols(samples)
    samples = ensure_int(samples, sample_cols)

    log("Reading breadth table")
    breadth = pd.read_csv(args.breadth, sep="\t").set_index("rep_contig")
    breadth = breadth.apply(pd.to_numeric, errors="coerce").fillna(0)

    log(f"Applying min_reads={args.min_reads} and min_breadth={args.min_breadth}")
    filtered_counts = apply_read_breadth_filter(
        samples[sample_cols],
        breadth,
        args.min_reads,
        args.min_breadth
    )

    filtered_out = samples.copy()
    filtered_out[sample_cols] = filtered_counts
    filtered_out.to_csv(args.out_filtered, sep="\t")

    stats = pd.DataFrame({
        "step": ["after_read_breadth"],
        "vOTU_detections": [(filtered_counts > 0).sum().sum()]
    })

    stats.to_csv(args.out_stats, sep="\t", index=False)

    log("Done.")


if __name__ == "__main__":
    main()