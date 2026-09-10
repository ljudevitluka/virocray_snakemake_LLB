#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from datetime import datetime
import argparse


def parse_args():
    p = argparse.ArgumentParser(description="Filter vOTUs by reads and CoverM coverage")
    p.add_argument("--samples", required=True)
    p.add_argument("--covered-fraction", required=True)
    p.add_argument("--covered-bases", required=True)
    p.add_argument("--sample-names", nargs="+", required=True)
    p.add_argument("--out-filtered", required=True)
    p.add_argument("--out-stats", required=True)
    p.add_argument("--min-reads", type=int, default=2)
    p.add_argument("--min-covered-fraction", type=float, default=0.0)
    p.add_argument("--min-covered-bases", type=float, default=0.0)
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


def load_coverage_matrix(file, sample_names):
    coverage = pd.read_csv(file, sep="\t")
    if "rep_contig" not in coverage.columns:
        raise ValueError(f"Coverage matrix {file} must contain a rep_contig column")

    missing_samples = [sample for sample in sample_names if sample not in coverage.columns]
    if missing_samples:
        raise ValueError(
            f"Coverage matrix {file} is missing sample columns: {missing_samples}"
        )

    coverage = coverage.set_index("rep_contig")
    coverage.index = coverage.index.astype(str).str.strip()
    return coverage[sample_names].apply(pd.to_numeric, errors="coerce").fillna(0)


def apply_read_coverage_filter(
    counts_df,
    covered_fraction,
    covered_bases,
    min_reads,
    min_covered_fraction,
    min_covered_bases,
):
    out = counts_df.copy()

    for sample in out.columns:
        covered_fraction_values = out.index.map(covered_fraction[sample]).fillna(0)
        covered_bases_values = out.index.map(covered_bases[sample]).fillna(0)

        low_reads = out[sample] < min_reads
        low_covered_fraction = covered_fraction_values < min_covered_fraction
        low_covered_bases = covered_bases_values < min_covered_bases

        out.loc[
            low_reads | low_covered_fraction | low_covered_bases,
            sample,
        ] = 0

    return out


def count_vOTU_detections(df):
    return (df > 0).sum().sum()


def main():
    args = parse_args()

    log("Reading samples table")
    samples = pd.read_csv(args.samples, sep="\t").set_index("rep_contig")

    sample_cols = get_readcount_cols(samples)
    samples = ensure_int(samples, sample_cols)

    log("Reading merged coverage metrics")
    covered_fraction = load_coverage_matrix(args.covered_fraction, args.sample_names)
    covered_bases = load_coverage_matrix(args.covered_bases, args.sample_names)

    log(
        f"Applying min_reads={args.min_reads}, "
        f"min_covered_fraction={args.min_covered_fraction}, "
        f"min_covered_bases={args.min_covered_bases}"
    )
    filtered_counts = apply_read_coverage_filter(
        samples[sample_cols],
        covered_fraction,
        covered_bases,
        args.min_reads,
        args.min_covered_fraction,
        args.min_covered_bases,
    )

    filtered_out = samples.copy()
    filtered_out[sample_cols] = filtered_counts
    filtered_out.to_csv(args.out_filtered, sep="\t")

    stats = pd.DataFrame({
        "step": ["after_read_coverage"],
        "vOTU_detections": [(filtered_counts > 0).sum().sum()]
    })

    stats.to_csv(args.out_stats, sep="\t", index=False)

    log("Done.")


if __name__ == "__main__":
    main()