rule subsample_reads:
    input:
        r1 = "resources/{sample_ID}/{sample_ID}_1.fq.gz",
        r2 = "resources/{sample_ID}/{sample_ID}_2.fq.gz"
    output:
        r1 = RESULTS_DIR + "/{sample_ID}/00_{sample_ID}_subsampled_1.fq.gz",
        r2 = RESULTS_DIR + "/{sample_ID}/00_{sample_ID}_subsampled_2.fq.gz",
        check = RESULTS_DIR + "/{sample_ID}/00_{sample_ID}_subsampled.done"
    params:
        n = 24678480
    log:
        logO = "logs/subsample/{sample_ID}.log",
        logE = "logs/subsample/{sample_ID}.err.log"
    conda:
        "../envs/seqtk_env.yaml"
    threads: 2
    shell:
        """
        mkdir -p {RESULTS_DIR}/{wildcards.sample_ID}

        seqtk sample -s100 {input.r1} {params.n} | gzip > {output.r1} 2>> {log.logE}
        seqtk sample -s100 {input.r2} {params.n} | gzip > {output.r2} 2>> {log.logE}

        touch {output.check}
        """
