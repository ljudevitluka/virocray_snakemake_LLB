rule merge_taxonomy_clusters:
    input:
        clusters = RESULTS_DIR + "/merged/clusters.tsv",
        taxonomy = RESULTS_DIR + "/merged/merged_viral_taxonomy.tsv"
    output:
        merged = RESULTS_DIR + "/merged/clusters_with_taxonomy.tsv"
    log:
        logO = "logs/merge_taxonomy_clusters/merge_taxonomy_clusters.log",
        logE = "logs/merge_taxonomy_clusters/merge_taxonomy_clusters.err.log"
    conda:
        "../envs/clustering_env.yaml"
    shell:
        """
        python workflow/scripts/merging_taxonomy_and_clusters_script.py \
          --clusters {input.clusters} \
          --taxonomy {input.taxonomy} \
          --output {output.merged} \
          > {log.logO} 2> {log.logE}
        """
