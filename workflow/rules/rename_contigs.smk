rule rename_contigs:
    input:
        fasta = "results/{sample}/02_{sample}_spades_rnaviral_contigs.fasta"
    output:
        renamed = "results/{sample}/02_{sample}_spades_rnaviral_contigs_renamed.fasta"
    log:
        log = "logs/rename_contigs/{sample}.log"
    conda:
        "../envs/rename_contigs_env.yaml"
    shell:
        """
        python workflow/scripts/renaming_contigs_script.py \
            --input {input.fasta} \
            --output {output.renamed} \
            --prefix {wildcards.sample} \
            > {log.log} 2>&1
        """
