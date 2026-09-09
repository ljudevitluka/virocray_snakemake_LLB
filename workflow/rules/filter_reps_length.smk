rule filter_reps_length:
    input:
        reps_fasta = RESULTS_DIR + "/merged/cluster_representatives.fasta",
        reps_tsv   = RESULTS_DIR + "/merged/cluster_representatives_list.tsv"
    output:
        reps_fasta = RESULTS_DIR + "/merged/cluster_representatives_filtered_length.fasta",
        reps_tsv   = RESULTS_DIR + "/merged/cluster_representatives_filtered_length.tsv",
        stats      = RESULTS_DIR + "/merged/cluster_representatives_filtered_length_stats.tsv"
    params:
        min_len = 500
    log:
        logO = "logs/filter_reps_length/filter_reps_length.log",
        logE = "logs/filter_reps_length/filter_reps_length.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/filter_reps_length.py \
          --reps-fasta {input.reps_fasta} \
          --reps-tsv {input.reps_tsv} \
          --out-fasta {output.reps_fasta} \
          --out-tsv {output.reps_tsv} \
          --stats {output.stats} \
          --min-len {params.min_len} \
          > {log.logO} 2> {log.logE}
        """
