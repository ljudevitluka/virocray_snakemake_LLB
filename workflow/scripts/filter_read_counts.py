#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from datetime import datetime
import argparse


def parse_args():
    p = argparse.ArgumentParser(description="Filter vOTUs by reads, breadth, and CoverM coverage")
    p.add_argument("--samples", required=True)
    p.add_argument("--breadth", required=True)
    p.add_argument("--coverm-files", nargs="+", required=True)
    p.add_argument("--sample-names", nargs="+", required=True)
    p.add_argument("--out-filtered", required=True)
    p.add_argument("--out-stats", required=True)
    p.add_argument("--min-reads", type=int, default=2)
    p.add_argument("--min-breadth", type=float, default=0.0)
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


def normalise_column_name(column):
    return "".join(str(column).lower().split()).replace("_", "")


def load_coverage_metrics(files, sample_names):
    if len(files) != len(sample_names):
        raise ValueError("CoverM files and sample names must have same length")

    metrics = {}
    for file, sample in zip(files, sample_names):
        df = pd.read_csv(file, sep="\t")
        columns = {normalise_column_name(c): c for c in df.columns}
        contig_col = columns.get("contig")
        fraction_col = columns.get("coveredfraction")
        bases_col = columns.get("coveredbases")
        if not contig_col or not fraction_col or not bases_col:
            raise ValueError(
                f"CoverM file {file} must contain Contig, Covered Fraction, and Covered Bases"
            )

        sample_metrics = df[[contig_col, fraction_col, bases_col]].copy()
        sample_metrics.columns = ["rep_contig", "covered_fraction", "covered_bases"]
        sample_metrics["rep_contig"] = (
            sample_metrics["rep_contig"].astype(str).str.strip()
        )
        sample_metrics["rep_contig"] = sample_metrics["rep_contig"].str.replace(
            r"_cluster_\d+$", "", regex=True
        )
        sample_metrics = sample_metrics.set_index("rep_contig")
        metrics[sample] = sample_metrics.apply(pd.to_numeric, errors="coerce").fillna(0)

    return metrics


def apply_read_breadth_filter(
    counts_df,
    breadth_df,
    coverage_metrics,
    min_reads,
    min_breadth,
    min_covered_fraction,
    min_covered_bases,
):
    out = counts_df.copy()

    for sample in out.columns:
        breadth_vals = (
            out.index.map(breadth_df[sample]).fillna(0)
            if sample in breadth_df.columns
            else 0
        )
        coverage = coverage_metrics.get(sample)
        covered_fraction = (
            out.index.map(coverage["covered_fraction"]).fillna(0)
            if coverage is not None
            else 0
        )
        covered_bases = (
            out.index.map(coverage["covered_bases"]).fillna(0)
            if coverage is not None
            else 0
        )

        low_reads = out[sample] < min_reads
        low_breadth = breadth_vals < min_breadth
        low_covered_fraction = covered_fraction < min_covered_fraction
        low_covered_bases = covered_bases < min_covered_bases

        out.loc[
            low_reads | low_breadth | low_covered_fraction | low_covered_bases,
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

    log("Reading breadth table")
    breadth = pd.read_csv(args.breadth, sep="\t").set_index("rep_contig")
    breadth = breadth.apply(pd.to_numeric, errors="coerce").fillna(0)

    log("Reading CoverM coverage metrics")
    coverage_metrics = load_coverage_metrics(args.coverm_files, args.sample_names)

    log(
        f"Applying min_reads={args.min_reads}, min_breadth={args.min_breadth}, "
        f"min_covered_fraction={args.min_covered_fraction}, "
        f"min_covered_bases={args.min_covered_bases}"
    )
    filtered_counts = apply_read_breadth_filter(
        samples[sample_cols],
        breadth,
        coverage_metrics,
        args.min_reads,
        args.min_breadth,
        args.min_covered_fraction,
        args.min_covered_bases,
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