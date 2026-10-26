rule merge_final_table:
    input:
        tpm      = RESULTS_DIR + "/merged/reps_tpm_samples_filtered.tsv",
        taxonomy = RESULTS_DIR + "/merged/clusters_with_taxonomy.tsv",
        blastx   = RESULTS_DIR + "/merged/cluster_representatives_blastx_reps_nr.tsv",
        annot    = RESULTS_DIR + "/merged/functional_annotation/reps_functional_annotation.tsv"
    output:
        table = RESULTS_DIR + "/merged/final_table_tpm_taxonomy_blastx.tsv"
    log:
        logO = "logs/merge_final_table/merge_final_table.log",
        logE = "logs/merge_final_table/merge_final_table.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/merge_final_table.py \
            --tpm {input.tpm} \
            --taxonomy {input.taxonomy} \
            --blastx {input.blastx} \
            --annotation {input.annot} \
            --output {output.table} \
            > {log.logO} 2> {log.logE}
        """
