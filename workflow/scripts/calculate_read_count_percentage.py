#!/usr/bin/env python3

import argparse
import gzip

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate per-contig read-count percentages"
    )
    parser.add_argument("--counts", required=True)
    parser.add_argument("--r1-files", nargs="+", required=True)
    parser.add_argument("--r2-files", nargs="+", required=True)
    parser.add_argument("--sample-names", nargs="+", required=True)
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def count_fastq_reads(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        records = sum(1 for _ in handle)

    if records % 4 != 0:
        raise ValueError(f"FASTQ file does not contain complete records: {path}")
    return records // 4


def main():
    args = parse_args()

    if not (
        len(args.r1_files)
        == len(args.r2_files)
        == len(args.sample_names)
    ):
        raise ValueError("R1 files, R2 files, and sample names must have the same length")

    counts = pd.read_csv(args.counts, sep="\t")
    if "rep_contig" not in counts.columns:
        raise ValueError("Counts table must contain a rep_contig column")

    missing_samples = [sample for sample in args.sample_names if sample not in counts.columns]
    if missing_samples:
        raise ValueError(f"Counts table is missing sample columns: {missing_samples}")

    counts = counts.set_index("rep_contig")
    counts[args.sample_names] = counts[args.sample_names].apply(
        pd.to_numeric, errors="coerce"
    ).fillna(0)

    total_reads = {}
    for sample, r1, r2 in zip(args.sample_names, args.r1_files, args.r2_files):
        total_reads[sample] = count_fastq_reads(r1) + count_fastq_reads(r2)

    percentages = counts[args.sample_names].copy().astype(float)
    for sample in args.sample_names:
        denominator = total_reads[sample]
        if denominator == 0:
            raise ValueError(f"Trimmed FASTQ files contain zero reads for {sample}")
        percentages[sample] = percentages[sample] / denominator * 100.0

    percentages.reset_index().to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()