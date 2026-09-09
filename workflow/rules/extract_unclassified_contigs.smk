rule extract_unclassified_contigs:
    input:
        fasta = RESULTS_DIR + "/merged/merged_contigs.fasta",
        taxonomy = RESULTS_DIR + "/merged/merged_taxonomy.tsv"
    output:
        unclassified_fasta = RESULTS_DIR + "/merged/unclassified_contigs.fasta"
    log:
        logO = "logs/extract_unclassified_contigs/extract_unclassified_contigs.log",
        logE = "logs/extract_unclassified_contigs/extract_unclassified_contigs.err.log"
    conda:
        "../envs/clustering_env.yaml"
    shell:
        """
        python workflow/scripts/extract_unclassified_contigs.py \
            --fasta {input.fasta} \
            --taxonomy {input.taxonomy} \
            --output {output.unclassified_fasta} \
            > {log.logO} 2> {log.logE}
        """
