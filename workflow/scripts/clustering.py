#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from Bio import SeqIO
import pandas as pd
import igraph as ig
import leidenalg as la
from datetime import datetime

def parse_args():
    p = argparse.ArgumentParser(description="Run Leiden clustering from an edge list and report stats.")
    p.add_argument("--fasta", required=True, help="FASTA with contigs (node list).")
    p.add_argument("--edges", required=True, help="Edges TSV with two columns: contig_A, contig_B.")
    p.add_argument("--clusters", required=True, help="Output clusters TSV (contig_id, cluster_id).")
    p.add_argument("--stats", required=True, help="Output stats TSV.")
    p.add_argument("--resolution", type=float, default=1.0, help="Leiden resolution parameter.")
    return p.parse_args()

def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)

def main():
    args = parse_args()

    log("Reading FASTA and preparing nodes")
    contigs = [rec.id for rec in SeqIO.parse(args.fasta, "fasta")]
    node_to_idx = {node: i for i, node in enumerate(contigs)}

    log("Reading edges")
    edges_df = pd.read_csv(args.edges, sep="\t", header=None, names=["a", "b"])
    edges = []
    skipped = 0
    for a, b in zip(edges_df["a"], edges_df["b"]):
        if a in node_to_idx and b in node_to_idx:
            edges.append((node_to_idx[a], node_to_idx[b]))
        else:
            skipped += 1

    log("Building graph")
    g = ig.Graph(n=len(contigs), edges=edges, directed=False)
    g.vs["name"] = contigs

    log("Calculating connected components")
    comps = g.components()
    n_comp = len(comps)
    comp_sizes = comps.sizes()
    comp_singletons = sum(1 for s in comp_sizes if s == 1)
    comp_multi = sum(1 for s in comp_sizes if s > 1)

    log("Running Leiden clustering")
    partition = la.find_partition(
        g,
        la.RBConfigurationVertexPartition,
        resolution_parameter=args.resolution
    )
    n_leiden = len(partition)
    leiden_sizes = [len(c) for c in partition]
    leiden_singletons = sum(1 for s in leiden_sizes if s == 1)
    leiden_multi = sum(1 for s in leiden_sizes if s > 1)

    log(f"Writing clusters to {args.clusters}")
    with open(args.clusters, "w") as out:
        out.write("contig_id\tcluster_id\n")
        for cid, community in enumerate(partition, start=1):
            cluster_name = f"cluster_{cid}"
            for vid in community:
                out.write(f"{g.vs[vid]['name']}\t{cluster_name}\n")

    log(f"Writing stats to {args.stats}")
    with open(args.stats, "w") as out:
        out.write("\t".join([
            "num_nodes",
            "num_edges",
            "edges_skipped_missing_nodes",
            "num_connected_components",
            "connected_component_singletons",
            "connected_component_multi",
            "num_leiden_clusters",
            "leiden_singletons",
            "leiden_multi",
            "leiden_resolution"
        ]) + "\n")

        out.write("\t".join(map(str, [
            g.vcount(),
            g.ecount(),
            skipped,
            n_comp,
            comp_singletons,
            comp_multi,
            n_leiden,
            leiden_singletons,
            leiden_multi,
            args.resolution
        ])) + "\n")

    log("Finished Leiden clustering script")

if __name__ == "__main__":
    main()
