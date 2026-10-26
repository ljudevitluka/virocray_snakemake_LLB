rule filter_read_counts:
    input:
        samples = RESULTS_DIR + "/merged/reps_read_counts_samples.tsv",
        covered_fraction = RESULTS_DIR + "/merged/reps_covered_fraction_samples.tsv",
        covered_bases = RESULTS_DIR + "/merged/reps_covered_bases_samples.tsv"
    output:
        filtered = RESULTS_DIR + "/merged/reps_read_counts_samples_filtered.tsv",
        stats    = RESULTS_DIR + "/merged/reps_read_counts_samples_filtered_stats.tsv"
    params:
        min_reads = config["filter_read_counts"]["min_reads"],
        min_covered_fraction = config["filter_read_counts"]["min_covered_fraction"],
        min_covered_bases = config["filter_read_counts"]["min_covered_bases"]
    log:
        logO = "logs/filter_read_counts/filter_breadth.log",
        logE = "logs/filter_read_counts/filter_breadth.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/filter_read_counts.py \
            --samples {input.samples} \
            --covered-fraction {input.covered_fraction} \
            --covered-bases {input.covered_bases} \
            --sample-names {sample_names} \
            --out-filtered {output.filtered} \
            --out-stats {output.stats} \
            --min-reads {params.min_reads} \
            --min-covered-fraction {params.min_covered_fraction} \
            --min-covered-bases {params.min_covered_bases} \
            > {log.logO} 2> {log.logE}
        """