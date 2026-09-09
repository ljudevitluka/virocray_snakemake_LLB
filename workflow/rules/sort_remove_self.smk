rule sort_remove_self:
    input:
        blast = RESULTS_DIR + "/merged/blastn_pairs.tsv"
    output:
        sorted = RESULTS_DIR + "/merged/blastn_pairs_sorted_no_self.tsv.gz"
    log:
        logO = "logs/sort_remove_self/sort_remove_self.log",
        logE = "logs/sort_remove_self/sort_remove_self.err.log"
    conda:
        "../envs/core_env.yaml"
    threads: 1
    shell:
        """
        python workflow/scripts/sort_remove_self.py \
            --input {input.blast} \
            --output {output.sorted} \
            > {log.logO} 2> {log.logE}
        """
