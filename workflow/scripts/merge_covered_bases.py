#!/usr/bin/env python3

import argparse

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge per-sample CoverM covered-bases values"
    )
    parser.add_argument("--viral-reps", required=True)
    parser.add_argument("--mode", choices=["samples"], required=True)
    parser.add_argument("--files", nargs="+", required=True)
    parser.add_argument("--sample-names", nargs="+", required=True)
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def normalise_column_name(column):
    return "".join(str(column).lower().split()).replace("_", "")


def load_covered_bases(file, sample_name):
    dataframe = pd.read_csv(file, sep="\t")
    dataframe.columns = dataframe.columns.str.strip()
    columns = {normalise_column_name(column): column for column in dataframe.columns}

    contig_column = columns.get("contig")
    bases_column = columns.get("coveredbases")
    if not contig_column or not bases_column:
        raise ValueError(
            f"CoverM file {file} must contain Contig and Covered Bases columns"
        )

    output = dataframe[[contig_column, bases_column]].copy()
    output.columns = ["rep_contig", sample_name]
    output["rep_contig"] = output["rep_contig"].astype(str).str.strip()
    output["rep_contig"] = output["rep_contig"].str.replace(
        r"_cluster_\d+$", "", regex=True
    )
    output[sample_name] = pd.to_numeric(output[sample_name], errors="coerce").fillna(0)

    return output.set_index("rep_contig")


def merge(files, sample_names):
    if len(files) != len(sample_names):
        raise ValueError("CoverM files and sample names must have the same length")

    merged = load_covered_bases(files[0], sample_names[0])
    for file, sample_name in zip(files[1:], sample_names[1:]):
        merged = merged.join(
            load_covered_bases(file, sample_name),
            how="outer",
        )

    return merged


def main():
    args = parse_args()

    representatives = pd.read_csv(args.viral_reps, sep="\t")
    representatives["rep_contig"] = representatives["rep_contig"].astype(str).str.strip()
    base = representatives.set_index("rep_contig")

    merged = merge(args.files, args.sample_names)
    output = base.join(merged, how="left").fillna(0)
    output.reset_index().to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()