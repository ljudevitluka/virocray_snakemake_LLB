rule calculate_breadth:
    input:
        bam = RESULTS_DIR + "/{sample}/08_{sample}_coverm_filtered_reps.bam"
    output:
        breadth = RESULTS_DIR + "/{sample}/09_{sample}_breadth.tsv"
    log:
        logO = "logs/calculate_breadth/{sample}.log",
        logE = "logs/calculate_breadth/{sample}.err.log"
    conda:
        "../envs/calculate_breadth_env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.breadth})

        cov_file=$(mktemp $(dirname {output.breadth})/tmp_genomecov_XXXXXX.txt)

        bedtools genomecov -d -ibam {input.bam} > "$cov_file" 2>> {log.logE}

        python workflow/scripts/calculate_breadth.py "$cov_file" {output.breadth} >> {log.logO} 2>> {log.logE}

        rm -f "$cov_file"
        """  