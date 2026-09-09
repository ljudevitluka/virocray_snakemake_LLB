rule calculate_rpkm:
    input:
        counts = RESULTS_DIR + "/merged/reps_read_counts_samples_filtered.tsv"
    output:
        rpkm_long = RESULTS_DIR + "/merged/rpkm_long.tsv",
        pa_long   = RESULTS_DIR + "/merged/presence_absence_long.tsv"
    log:
        logO = "logs/calculate_rpkm/calculate_rpkm.log",
        logE = "logs/calculate_rpkm/calculate_rpkm.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/calculate_rpkm.py \
            --input {input.counts} \
            --output_rpkm {output.rpkm_long} \
            --output_pa {output.pa_long} \
            > {log.logO} 2> {log.logE}
        """
