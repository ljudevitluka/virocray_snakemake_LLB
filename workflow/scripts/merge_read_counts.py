#!/usr/bin/env python3

import argparse
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument("--viral-reps", required=True)
    p.add_argument("--mode", required=True)

    # samples mode
    p.add_argument("--files", nargs="*")
    p.add_argument("--sample-names", nargs="*")

    # controls mode
    p.add_argument("--nki-files", nargs="*")
    p.add_argument("--nki-names", nargs="*")
    p.add_argument("--carry-files", nargs="*")
    p.add_argument("--carry-names", nargs="*")

    p.add_argument("--out")
    p.add_argument("--out-nki")
    p.add_argument("--out-carry")

    return p.parse_args()


def load_coverm(file, sample_name):
    df = pd.read_csv(file, sep="\t")
    df.columns = df.columns.str.strip()

    read_cols = [c for c in df.columns if "Read" in c and "Count" in c]
    if not read_cols:
        raise ValueError(f"No Read Count column in file: {file}")

    read_col = read_cols[0]

    out = df[["Contig", read_col]].copy()
    out.columns = ["rep_contig", sample_name]

    out["rep_contig"] = out["rep_contig"].astype(str).str.strip()
    out["rep_contig"] = out["rep_contig"].str.replace(
        r"_cluster_\d+$", "", regex=True
    )

    return out.set_index("rep_contig")


def merge(files, names):
    if len(files) != len(names):
        raise ValueError("Files and sample names must have same length")

    dfs = [load_coverm(f, n) for f, n in zip(files, names)]

    out = dfs[0]
    for df in dfs[1:]:
        out = out.join(df, how="outer")

    return out


def main():
    args = parse_args()

    reps = pd.read_csv(args.viral_reps, sep="\t")
    reps["rep_contig"] = reps["rep_contig"].astype(str).str.strip()
    base = reps.set_index("rep_contig")

    if args.mode == "samples":
        merged = merge(args.files, args.sample_names)
        out = base.join(merged, how="left").fillna(0)
        out.reset_index().to_csv(args.out, sep="\t", index=False)

    elif args.mode == "controls":
        nki = merge(args.nki_files, args.nki_names)
        carry = merge(args.carry_files, args.carry_names)

        nki_out = base.join(nki, how="left").fillna(0).reset_index()
        carry_out = base.join(carry, how="left").fillna(0).reset_index()

        nki_out.to_csv(args.out_nki, sep="\t", index=False)
        carry_out.to_csv(args.out_carry, sep="\t", index=False)

    else:
        raise ValueError(f"Unknown mode: {args.mode}")


if __name__ == "__main__":
    main()