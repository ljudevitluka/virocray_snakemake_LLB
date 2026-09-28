"""
visualize_clustering.smk
------------------------
Generates clustering QC figures after the Leiden clustering step.

Inputs  (produced by upstream rules):
  filter_edges rule  → edges_filtered.tsv
  clustering rule    → clusters.tsv, clusters_stats.tsv
  merge_contigs rule → merged_viral_contigs.fasta

Output:
  results/merged/clustering_qc/clustering_qc.pdf   (main figure)

Add the output path to rule all: to include it in the default target.
"""

rule visualize_clustering:
    input:
        fasta    = RESULTS_DIR + "/merged/merged_viral_contigs.fasta",
        edges    = RESULTS_DIR + "/merged/edges_filtered.tsv",
        clusters = RESULTS_DIR + "/merged/clusters.tsv",
        stats    = RESULTS_DIR + "/merged/clusters_stats.tsv"
    output:
        figure   = RESULTS_DIR + "/merged/clustering_qc/clustering_qc.pdf"
    params:
        outdir   = RESULTS_DIR + "/merged/clustering_qc",
        prefix   = "clustering_qc",
        fmt      = "pdf"
    log:
        logO = "logs/visualize_clustering/visualize_clustering.log",
        logE = "logs/visualize_clustering/visualize_clustering.err.log"
    conda:
        "../envs/visualize_clustering_env.yaml"
    shell:
        """
        python workflow/scripts/visualize_clustering.py \
            --fasta    {input.fasta}    \
            --edges    {input.edges}    \
            --clusters {input.clusters} \
            --stats    {input.stats}    \
            --outdir   {params.outdir}  \
            --prefix   {params.prefix}  \
            --format   {params.fmt}     \
            > {log.logO} 2> {log.logE}
        """
