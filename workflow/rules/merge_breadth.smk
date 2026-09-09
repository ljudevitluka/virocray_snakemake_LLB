rule merge_breadth:
    input:
        breadth_files = expand(
            RESULTS_DIR + "/{sample}/09_{sample}_breadth.tsv",
            sample=samples["sample"].tolist()
        )
    output:
        breadth_matrix = RESULTS_DIR + "/merged/merged_breadth.tsv"
    log:
        logO = "logs/merge_breadth/merge_breadth.log",
        logE = "logs/merge_breadth/merge_breadth.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/merge_breadth.py \
            --breadth-files {input.breadth_files} \
            --out {output.breadth_matrix} \
            > {log.logO} 2> {log.logE}
        """