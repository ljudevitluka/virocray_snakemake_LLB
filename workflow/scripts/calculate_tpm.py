#!/usr/bin/env python3

import argparse
import os

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate TPM from filtered per-contig read counts"
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def get_sample_columns(dataframe):
    metadata = {"rep_contig", "cluster_id", "rep_contig_length", "singleton"}
    return [column for column in dataframe.columns if column not in metadata]


def main():
    args = parse_args()
    dataframe = pd.read_csv(args.input, sep="\t")

    required = {"rep_contig", "rep_contig_length"}
    missing = required.difference(dataframe.columns)
    if missing:
        raise ValueError(f"Input is missing required columns: {sorted(missing)}")

    sample_columns = get_sample_columns(dataframe)
    if not sample_columns:
        raise ValueError("No sample columns detected")

    lengths_kb = pd.to_numeric(dataframe["rep_contig_length"], errors="coerce") / 1000
    if lengths_kb.isna().any() or (lengths_kb <= 0).any():
        raise ValueError("rep_contig_length must contain positive numeric values")

    counts = dataframe[sample_columns].apply(pd.to_numeric, errors="coerce").fillna(0)
    reads_per_kb = counts.div(lengths_kb, axis=0)
    scaling_factors = reads_per_kb.sum(axis=0)
    if (scaling_factors <= 0).any():
        invalid = scaling_factors[scaling_factors <= 0].index.tolist()
        raise ValueError(f"Cannot calculate TPM for samples with no retained signal: {invalid}")

    tpm = reads_per_kb.div(scaling_factors, axis=1) * 1_000_000
    output = dataframe[[column for column in dataframe.columns if column not in sample_columns]].copy()
    output[sample_columns] = tpm

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()