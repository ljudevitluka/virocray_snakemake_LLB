rule filter_read_counts:
    input:
        samples = RESULTS_DIR + "/merged/reps_read_counts_samples.tsv",
        breadth = RESULTS_DIR + "/merged/merged_breadth.tsv",
        coverm = expand(
            RESULTS_DIR + "/{sample}/08_{sample}_coverm_filtered_reps.tsv",
            sample=samples["sample"].tolist()
        )
    output:
        filtered = RESULTS_DIR + "/merged/reps_read_counts_samples_filtered.tsv",
        stats    = RESULTS_DIR + "/merged/reps_read_counts_samples_filtered_stats.tsv"
    params:
        min_reads = 2,
        min_breadth = 0.25,
        min_covered_fraction = 0.25,
        min_covered_bases = 400
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
            --coverm-files {input.coverm} \
            --sample-names {sample_names} \
            --out-filtered {output.filtered} \
            --out-stats {output.stats} \
            --min-reads {params.min_reads} \
            --min-breadth {params.min_breadth} \
            --min-covered-fraction {params.min_covered_fraction} \
            --min-covered-bases {params.min_covered_bases} \
            > {log.logO} 2> {log.logE}
        """