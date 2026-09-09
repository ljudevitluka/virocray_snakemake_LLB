rule merge_contigs:
    input:
        expand(
            RESULTS_DIR + "/{sample}/04_{sample}_viral_contigs.fasta",
            sample=samples["sample"].tolist()
        )
    output:
        merged=RESULTS_DIR + "/merged/merged_viral_contigs.fasta"
    log:
        "logs/merge_viral_contigs/merge_viral_contigs.log"
    shell:
        """
        mkdir -p {RESULTS_DIR}/merged
        cat {input} > {output.merged}
        """