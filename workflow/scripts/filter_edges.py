#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
from datetime import datetime
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser(description="Compute ANI/AF/maxHSP per canonical contig pair and filter to edges")
    p.add_argument("--blast_sorted_gz", required=True, help="Sorted BLAST TSV.gz with canonical columns A,B prepended")
    p.add_argument("--fasta", required=True, help="FASTA with contigs (for lengths)")
    p.add_argument("--out_pairs", required=True, help="Output filtered pairs TSV")
    p.add_argument("--out_edges", required=True, help="Output edges TSV (A\\tB)")
    p.add_argument("--min_ani", type=float, default=90.0, help="Minimum ANI in percent")
    p.add_argument("--filter_mode", choices=["HSP", "AF"], required=True, help="Use 'HSP' or 'AF' filter")
    p.add_argument("--min_hsp_len", type=int, default=150, help="Minimum HSP length if filter_mode=HSP")
    p.add_argument("--min_af", type=float, default=0.0, help="Minimum aligned fraction if filter_mode=AF")
    return p.parse_args()

def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)

def read_fasta_lengths(fasta_path):
    lengths = {}
    cur, cur_len = None, 0
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if cur is not None:
                    lengths[cur] = cur_len
                cur = line[1:].split()[0]
                cur_len = 0
            else:
                cur_len += len(line.strip())
        if cur is not None:
            lengths[cur] = cur_len
    return lengths

def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return merged

def intervals_total_len(intervals):
    return sum(e - s + 1 for s, e in intervals)

def compute_ani(hsps):
    num = sum(h["len"] * h["pid"] for h in hsps)
    den = sum(h["len"] for h in hsps)
    return num / den if den else 0.0

def compute_af(hsps, lenA, lenB):
    q_int = merge_intervals([h["qcoords"] for h in hsps])
    s_int = merge_intervals([h["scoords"] for h in hsps])
    afA = intervals_total_len(q_int) / lenA if lenA > 0 else 0.0
    afB = intervals_total_len(s_int) / lenB if lenB > 0 else 0.0
    return afA, afB

def flush_pair(A, B, hsps, lengths, args, out_pairs_f, out_edges_f):
    if not hsps:
        return
    lenA = lengths.get(A, 0)
    lenB = lengths.get(B, 0)
    ani = compute_ani(hsps) * 100.0
    max_hsp = max(h["len"] for h in hsps)
    afA, afB = compute_af(hsps, lenA, lenB)
    min_af = min(afA, afB)

    if ani < args.min_ani:
        return
    if args.filter_mode == "HSP" and max_hsp < args.min_hsp_len:
        return
    if args.filter_mode == "AF" and min_af < args.min_af:
        return

    out_pairs_f.write(
        "\t".join([A, B, f"{ani:.4f}", str(max_hsp), f"{afA:.6f}", f"{afB:.6f}", f"{min_af:.6f}", str(len(hsps))]) + "\n"
    )
    out_edges_f.write(f"{A}\t{B}\n")

def main():
    args = parse_args()
    log("Reading contig lengths")
    lengths = read_fasta_lengths(args.fasta)

    if not (0.0 <= args.min_af <= 1.0):
        raise ValueError("--min_af must be between 0 and 1")

    log(f"Processing BLAST results: {args.blast_sorted_gz}")
    with open(args.out_pairs, "w") as out_pairs_f, open(args.out_edges, "w") as out_edges_f:
        out_pairs_f.write("\t".join(["contig_A","contig_B","ANI_percent","max_hsp_len","AF_A","AF_B","AF_min","num_hsps"]) + "\n")
        current, hsps = None, []

        with gzip.open(args.blast_sorted_gz, "rt") as f:
            for line in f:
                line = line.rstrip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 10:
                    continue

                A, B = parts[0], parts[1]
                pident = float(parts[4]) / 100.0
                alen = int(float(parts[5]))
                qstart, qend = sorted([int(parts[6]), int(parts[7])])
                sstart, send = sorted([int(parts[8]), int(parts[9])])
                pair = (A, B)

                if current is None:
                    current = pair
                if pair != current:
                    flush_pair(current[0], current[1], hsps, lengths, args, out_pairs_f, out_edges_f)
                    current, hsps = pair, []

                hsps.append({"pid": pident, "len": alen, "qcoords": [qstart, qend], "scoords": [sstart, send]})

        if current is not None:
            flush_pair(current[0], current[1], hsps, lengths, args, out_pairs_f, out_edges_f)

    log(f"Finished writing pairs to {args.out_pairs} and edges to {args.out_edges}")

if __name__ == "__main__":
    main()
