rule fill_blastx_missing:
    input:
        fasta = RESULTS_DIR + "/merged/cluster_representatives_filtered_length.fasta",
        blast = RESULTS_DIR + "/merged/cluster_representatives_blastx_reps_nr.raw.tsv"

    output:
        tsv = RESULTS_DIR + "/merged/cluster_representatives_blastx_reps_nr.tsv"

    log:
        logO = "logs/fill_blastx_missing.log",
        logE = "logs/fill_blastx_missing.err.log"

    conda:
        "../envs/core_env.yaml"

    shell:
        """
        python workflow/scripts/fill_missing_blastx.py \
            {input.fasta} \
            {input.blast} \
            {output.tsv} \
            > {log.logO} 2> {log.logE}
        """