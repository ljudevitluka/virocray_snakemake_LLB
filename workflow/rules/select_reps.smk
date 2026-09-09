rule select_reps:
    input:
        fasta    = RESULTS_DIR + "/merged/merged_viral_contigs.fasta",
        clusters = RESULTS_DIR + "/merged/clusters.tsv"
    output:
        reps_fasta = RESULTS_DIR + "/merged/cluster_representatives.fasta",
        reps_tsv   = RESULTS_DIR + "/merged/cluster_representatives_list.tsv",
        stats      = RESULTS_DIR + "/merged/cluster_representatives_stats.tsv"
    log:
        logO = "logs/select_reps/select_reps.log",
        logE = "logs/select_reps/select_reps.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/select_reps.py \
          --fasta {input.fasta} \
          --clusters {input.clusters} \
          --reps-fasta {output.reps_fasta} \
          --reps-tsv {output.reps_tsv} \
          --stats {output.stats} \
          > {log.logO} 2> {log.logE}
        """
