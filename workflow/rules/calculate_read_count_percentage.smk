rule calculate_read_count_percentage:
    input:
        counts = RESULTS_DIR + "/merged/reps_read_counts_samples.tsv",
        r1 = expand(
            RESULTS_DIR + "/{sample}/01_{sample}_trim_1_paired.fq.gz",
            sample=samples["sample"].tolist()
        ),
        r2 = expand(
            RESULTS_DIR + "/{sample}/01_{sample}_trim_2_paired.fq.gz",
            sample=samples["sample"].tolist()
        )
    output:
        percentages = RESULTS_DIR + "/merged/reps_read_count_percentage_samples.tsv"
    log:
        logO = "logs/read_count_percentage/read_count_percentage.log",
        logE = "logs/read_count_percentage/read_count_percentage.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/calculate_read_count_percentage.py \
            --counts {input.counts} \
            --r1-files {input.r1} \
            --r2-files {input.r2} \
            --sample-names {sample_names} \
            --out {output.percentages} \
            > {log.logO} 2> {log.logE}
        """