#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
visualize_clustering.py
-----------------------
Generate clustering QC plots after the Leiden clustering step.

Produces a single multi-panel PDF / PNG figure with three panels:
  1. Clustered vs. non-clustered contigs (stacked bar)
  2. Histogram of the number of edges per contig (node degree distribution)
  3. Histogram of cluster sizes (number of contigs per cluster)

Inputs (all produced by upstream rules):
  --fasta      : merged_viral_contigs.fasta          (merge_contigs rule)
  --edges      : edges_filtered.tsv                  (filter_edges rule)
  --clusters   : clusters.tsv                        (clustering rule)
  --stats      : clusters_stats.tsv                  (clustering rule)

Output:
  --outdir     : directory where the figure(s) are written
  --prefix     : filename prefix (default: clustering_qc)
  --format     : pdf | png | svg  (default: pdf)
"""

import argparse
import os
from collections import Counter
from datetime import datetime

import pandas as pd
from Bio import SeqIO
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np


# ── helpers ──────────────────────────────────────────────────────────────────

def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)


def parse_args():
    p = argparse.ArgumentParser(description="Leiden clustering visualisation")
    p.add_argument("--fasta",    required=True, help="FASTA: merged_viral_contigs.fasta")
    p.add_argument("--edges",    required=True, help="TSV: edges_filtered.tsv")
    p.add_argument("--clusters", required=True, help="TSV: clusters.tsv (contig_id, cluster_id)")
    p.add_argument("--stats",    required=True, help="TSV: clusters_stats.tsv")
    p.add_argument("--outdir",   required=True, help="Output directory")
    p.add_argument("--prefix",   default="clustering_qc", help="Output filename prefix")
    p.add_argument("--format",   default="pdf", choices=["pdf", "png", "svg"],
                   help="Output figure format")
    return p.parse_args()


# ── colour palette (consistent with project style) ───────────────────────────
BLUE   = "#2E86AB"   # clustered / multi-member clusters
ORANGE = "#E07A5F"   # non-clustered / singletons
GREY   = "#8D8D92"


# ── data loading ─────────────────────────────────────────────────────────────

def load_fasta_lengths(fasta_path):
    log(f"Reading FASTA: {fasta_path}")
    lengths = {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}
    log(f"  {len(lengths):,} contigs loaded")
    return lengths


def load_edges(edges_path):
    log(f"Reading edges: {edges_path}")
    df = pd.read_csv(edges_path, sep="\t", header=None, names=["a", "b"])
    log(f"  {len(df):,} edges loaded")
    return df


def load_clusters(clusters_path):
    log(f"Reading clusters: {clusters_path}")
    df = pd.read_csv(clusters_path, sep="\t")
    if not {"contig_id", "cluster_id"}.issubset(df.columns):
        raise ValueError("clusters TSV must contain columns: contig_id, cluster_id")
    log(f"  {len(df):,} contig–cluster rows loaded")
    return df


def load_stats(stats_path):
    log(f"Reading stats: {stats_path}")
    df = pd.read_csv(stats_path, sep="\t")
    return df


# ── derived statistics ────────────────────────────────────────────────────────

def node_degree(edges_df, all_contigs):
    """Return Series: contig -> number of edges (degree). Contigs with 0 edges included."""
    deg = Counter()
    for _, row in edges_df.iterrows():
        deg[row["a"]] += 1
        deg[row["b"]] += 1
    # ensure every contig is represented
    for c in all_contigs:
        if c not in deg:
            deg[c] = 0
    return pd.Series(deg)


def cluster_sizes(clusters_df):
    """Return Series of cluster sizes (number of contigs per cluster)."""
    return clusters_df.groupby("cluster_id").size()


# ── plotting helpers ──────────────────────────────────────────────────────────

def _despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _add_count_label(ax, patches, fontsize=8):
    """Annotate bar tops with count."""
    for p in patches:
        h = p.get_height()
        if h == 0:
            continue
        ax.text(
            p.get_x() + p.get_width() / 2,
            h + ax.get_ylim()[1] * 0.01,
            f"{int(h):,}",
            ha="center", va="bottom", fontsize=fontsize
        )


# ── panel 1 : clustered vs non-clustered ─────────────────────────────────────

def panel_clustered_vs_not(ax, clusters_df, all_contigs):
    """
    Stacked bar: total contigs = clustered (in multi-member clusters)
                               + singletons (alone in their cluster)
    """
    sizes = cluster_sizes(clusters_df)
    singleton_clusters = set(sizes[sizes == 1].index)
    singleton_contigs  = set(clusters_df[clusters_df["cluster_id"].isin(singleton_clusters)]["contig_id"])
    multi_contigs      = set(clusters_df["contig_id"]) - singleton_contigs
    unclustered        = set(all_contigs) - set(clusters_df["contig_id"])  # should be 0

    n_multi      = len(multi_contigs)
    n_singletons = len(singleton_contigs)
    n_unclust    = len(unclustered)

    x      = [0]
    width  = 0.45

    b1 = ax.bar(x, [n_multi],      width=width, label=f"Multi-member clusters ({n_multi:,})",
                color=BLUE, edgecolor="white", linewidth=0.5)
    b2 = ax.bar(x, [n_singletons], width=width, label=f"Singletons ({n_singletons:,})",
                color=ORANGE, edgecolor="white", linewidth=0.5, bottom=[n_multi])
    if n_unclust:
        b3 = ax.bar(x, [n_unclust], width=width, label=f"Not in clusters ({n_unclust:,})",
                    color=GREY, edgecolor="white", linewidth=0.5,
                    bottom=[n_multi + n_singletons])

    total = n_multi + n_singletons + n_unclust
    ax.set_ylim(0, total * 1.18)
    ax.set_xticks([])
    ax.set_ylabel("Number of contigs", fontsize=10)
    ax.set_title("Clustered vs. singleton contigs", fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, frameon=False)
    _despine(ax)

    # annotate stacked segments
    for bar_list, bottom in [(b1, 0), (b2, n_multi)]:
        for p in bar_list:
            h = p.get_height()
            if h == 0:
                continue
            ax.text(p.get_x() + p.get_width() / 2,
                    bottom + h / 2,
                    f"{int(h):,}",
                    ha="center", va="center", fontsize=8,
                    color="white", fontweight="bold")


# ── panel 2 : node degree (edges per contig) histogram ───────────────────────

def panel_degree_histogram(ax, edges_df, all_contigs):
    """
    Histogram of node degree (number of edges per contig).
    Separate colour for degree-0 nodes (no edges at all).
    """
    deg = node_degree(edges_df, all_contigs)
    zero_deg  = (deg == 0).sum()
    nonzero   = deg[deg > 0]

    if len(nonzero) == 0:
        ax.text(0.5, 0.5, "No edges found", ha="center", va="center",
                transform=ax.transAxes, fontsize=10)
        ax.set_title("Node degree distribution", fontsize=11, fontweight="bold")
        return

    # choose bins: log-spaced makes long tails visible
    max_deg = int(nonzero.max())
    if max_deg <= 50:
        bins = range(1, max_deg + 2)
    else:
        bins = np.unique(np.floor(np.geomspace(1, max_deg + 1, 40)).astype(int))

    ax.hist(nonzero, bins=bins, color=BLUE, edgecolor="white", linewidth=0.4,
            label=f"Contigs with ≥1 edge (n={len(nonzero):,})")

    # annotate zero-degree count in a text box
    ax.axvline(x=0, color="grey", linewidth=0)  # invisible; just anchors the box
    ax.text(0.97, 0.95,
            f"No edges: {zero_deg:,} contigs",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color=ORANGE,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=ORANGE, lw=0.8))

    ax.set_xlabel("Number of edges (degree)", fontsize=10)
    ax.set_ylabel("Number of contigs", fontsize=10)
    ax.set_title("Node degree distribution", fontsize=11, fontweight="bold")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.legend(fontsize=8, frameon=False)
    _despine(ax)


# ── panel 3 : cluster size distribution ──────────────────────────────────────

def panel_cluster_size_histogram(ax, clusters_df):
    """
    Histogram of cluster sizes (number of contigs per cluster).
    Singletons shown separately in ORANGE; multi-member in BLUE.
    """
    sizes = cluster_sizes(clusters_df)
    singles = sizes[sizes == 1]
    multi   = sizes[sizes > 1]

    if len(multi) == 0:
        ax.text(0.5, 0.5, "No multi-member clusters", ha="center", va="center",
                transform=ax.transAxes, fontsize=10)
        ax.set_title("Cluster size distribution", fontsize=11, fontweight="bold")
        return

    max_size = int(multi.max())
    if max_size <= 50:
        bins = range(2, max_size + 2)
    else:
        bins = np.unique(np.floor(np.geomspace(2, max_size + 1, 40)).astype(int))

    ax.hist(multi, bins=bins, color=BLUE, edgecolor="white", linewidth=0.4,
            label=f"Multi-member (n={len(multi):,} clusters)")

    ax.text(0.97, 0.95,
            f"Singletons: {len(singles):,} clusters",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color=ORANGE,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=ORANGE, lw=0.8))

    ax.set_xlabel("Cluster size (contigs per cluster)", fontsize=10)
    ax.set_ylabel("Number of clusters", fontsize=10)
    ax.set_title("Cluster size distribution\n(multi-member clusters only)", fontsize=11, fontweight="bold")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.legend(fontsize=8, frameon=False)
    _despine(ax)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    lengths_dict = load_fasta_lengths(args.fasta)
    all_contigs  = list(lengths_dict.keys())
    edges_df     = load_edges(args.edges)
    clusters_df  = load_clusters(args.clusters)
    stats_df     = load_stats(args.stats)

    # ── print summary to stdout ───────────────────────────────────────────
    log("=== Clustering summary ===")
    for col in stats_df.columns:
        log(f"  {col}: {stats_df[col].iloc[0]}")
    log(f"  Total contigs in FASTA : {len(all_contigs):,}")
    log(f"  Total edges            : {len(edges_df):,}")

    # ── figure ────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    fig.suptitle("ViroCray – Leiden clustering QC", fontsize=13, fontweight="bold", y=1.01)

    panel_clustered_vs_not(axes[0], clusters_df, all_contigs)
    panel_degree_histogram(axes[1], edges_df, all_contigs)
    panel_cluster_size_histogram(axes[2], clusters_df)

    plt.tight_layout()

    out_path = os.path.join(args.outdir, f"{args.prefix}.{args.format}")
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    log(f"Figure saved to {out_path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
