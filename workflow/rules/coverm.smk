rule coverm:
    input:
        bam = RESULTS_DIR + "/{sample}/07_{sample}_bbmap_reps.bam",
        bai = RESULTS_DIR + "/{sample}/07_{sample}_bbmap_reps.bam.bai"
    output:
        tsv     = RESULTS_DIR + "/{sample}/08_{sample}_coverm_filtered_reps.tsv",
        filtered_bam = RESULTS_DIR + "/{sample}/08_{sample}_coverm_filtered_reps.bam",
        check   = RESULTS_DIR + "/{sample}/08_{sample}_coverm_filtered_reps.done"
    log:
        logO = "logs/coverm/{sample}.log",
        logE = "logs/coverm/{sample}.err.log"
    conda:
        "../envs/mapping_env.yaml"
    threads: 12
    shell:
        """
        coverm filter \
            --bam-files {input.bam} \
            --threads {threads} \
            --min-read-percent-identity 95 \
            --min-read-aligned-percent 95 \
            --output-bam-files {output.filtered_bam} \
            > {log.logO} 2> {log.logE}

        samtools index {output.filtered_bam} 2>> {log.logE}

        coverm contig \
            --bam-files {output.filtered_bam} \
            --threads {threads} \
            --methods count tpm rpkm \
            --output-file {output.tsv} \
            >> {log.logO} 2>> {log.logE}

        touch {output.check}
        """
