#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import sys

def parse_args():
    p = argparse.ArgumentParser(
        description="Remove self-hits, canonicalize BLAST pairs, and sort for ANI/AF calculation"
    )
    p.add_argument("--input", required=True, help="Input BLAST TSV")
    p.add_argument("--output", required=True, help="Output sorted TSV.gz")
    return p.parse_args()

def log(msg):
    print(f"[INFO] {msg}", flush=True)

def main():
    args = parse_args()
    log(f"Reading input BLAST TSV from {args.input}")

    records = []

    with open(args.input, "r") as fin:
        for i, line in enumerate(fin, 1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                log(f"[WARN] Skipping malformed line {i}")
                continue

            q, s = parts[0], parts[1]

            if q == s:
                continue

            a, b = sorted([q, s])

            records.append((a, b, parts))

    log(f"Total valid pairs collected: {len(records)}")
    log("Sorting canonical pairs")
    records.sort(key=lambda x: (x[0], x[1]))

    log(f"Writing sorted output to {args.output}")
    with gzip.open(args.output, "wt") as fout:
        for a, b, parts in records:
            fout.write("\t".join([a, b] + parts) + "\n")

    log("Done.")

if __name__ == "__main__":
    main()
