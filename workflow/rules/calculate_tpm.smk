rule calculate_tpm:
    input:
        counts = RESULTS_DIR + "/merged/reps_read_counts_samples_filtered.tsv"
    output:
        tpm = RESULTS_DIR + "/merged/reps_tpm_samples_filtered.tsv"
    log:
        logO = "logs/calculate_tpm/calculate_tpm.log",
        logE = "logs/calculate_tpm/calculate_tpm.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/calculate_tpm.py \
            --input {input.counts} \
            --output {output.tpm} \
            > {log.logO} 2> {log.logE}
        """