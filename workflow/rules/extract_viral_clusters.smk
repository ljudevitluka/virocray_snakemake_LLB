rule extract_viral_clusters:
    input:
        merged = "results/merged/clusters_with_taxonomy.tsv",
        fasta = "results/merged/merged_contigs.fasta"
    output:
        viral_clusters = "results/merged/viral_clusters.tsv"
    log:
        logO = "logs/extract_viral_clusters/extract_viral_clusters.log",
        logE = "logs/extract_viral_clusters/extract_viral_clusters.err.log"
    conda:
        "../envs/clustering_env.yaml"
    shell:
        """
        python workflow/scripts/extracting_viral_clusters_script.py \
          --input {input.merged} \
          --fasta {input.fasta} \
          --output {output.viral_clusters} \
          > {log.logO} 2> {log.logE}
        """
