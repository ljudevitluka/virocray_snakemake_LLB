rule filter_read_counts:
    input:
        samples = RESULTS_DIR + "/merged/reps_read_counts_samples.tsv",
        breadth = RESULTS_DIR + "/merged/merged_breadth.tsv"
    output:
        filtered = RESULTS_DIR + "/merged/reps_read_counts_samples_filtered.tsv",
        stats    = RESULTS_DIR + "/merged/reps_read_counts_samples_filtered_stats.tsv"
    params:
        min_reads = 2,
        min_breadth = 0.25
    log:
        logO = "logs/filter_read_counts/filter_breadth.log",
        logE = "logs/filter_read_counts/filter_breadth.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/filter_read_counts.py \
            --samples {input.samples} \
            --breadth {input.breadth} \
            --out-filtered {output.filtered} \
            --out-stats {output.stats} \
            --min-reads {params.min_reads} \
            --min-breadth {params.min_breadth} \
            > {log.logO} 2> {log.logE}
        """