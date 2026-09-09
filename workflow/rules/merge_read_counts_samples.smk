rule merge_read_counts_samples:
    input:
        reps = RESULTS_DIR + "/merged/cluster_representatives_filtered_length.tsv",
        coverm = expand(
            RESULTS_DIR + "/{sample}/08_{sample}_coverm_filtered_reps.tsv",
            sample=samples["sample"].tolist()
        )
    output:
        samples = RESULTS_DIR + "/merged/reps_read_counts_samples.tsv"
    log:
        logO = "logs/merge_read_counts/merge_read_counts_samples.log",
        logE = "logs/merge_read_counts/merge_read_counts_samples.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/merge_read_counts.py \
            --viral-reps {input.reps} \
            --mode samples \
            --files {input.coverm} \
            --sample-names {sample_names} \
            --out {output.samples} \
            > {log.logO} 2> {log.logE}
        """