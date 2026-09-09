rule clustering:
    input:
        fasta = RESULTS_DIR + "/merged/merged_viral_contigs.fasta",
        edges = RESULTS_DIR + "/merged/edges_filtered.tsv"
    output:
        clusters = RESULTS_DIR + "/merged/clusters.tsv",
        stats    = RESULTS_DIR + "/merged/clusters_stats.tsv"
    params:
        resolution = 1.0
    log:
        logO = "logs/clustering/clustering.log",
        logE = "logs/clustering/clustering.err.log"
    conda:
        "../envs/clustering_env.yaml"
    shell:
        """
        python workflow/scripts/clustering.py \
            --fasta {input.fasta} \
            --edges {input.edges} \
            --clusters {output.clusters} \
            --stats {output.stats} \
            --resolution {params.resolution} \
            > {log.logO} 2> {log.logE}
        """
