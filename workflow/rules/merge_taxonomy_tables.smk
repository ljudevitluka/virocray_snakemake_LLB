rule merge_taxonomy_tables:
    input:
        taxonomy=expand(
            RESULTS_DIR + "/{sample}/04_{sample}_viral_taxonomy.tsv",
            sample=sorted(samples["sample"].tolist())
        )
    output:
        merged=RESULTS_DIR + "/merged/merged_viral_taxonomy.tsv"
    log:
        "logs/merge_taxonomy_tables/merge_taxonomy_tables.log"
    shell:
        """
        mkdir -p {RESULTS_DIR}/merged
        cat {input.taxonomy} > {output.merged} 2> {log}
        """