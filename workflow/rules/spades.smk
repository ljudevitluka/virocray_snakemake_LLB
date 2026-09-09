rule spades:
    input:
        r1="results/{sample}/01_{sample}_trim_1_paired.fq.gz",
        r2="results/{sample}/01_{sample}_trim_2_paired.fq.gz"
    output:
        d=temp(directory("results/{sample}/02_{sample}_spades_rnaviral_dir")),
        f="results/{sample}/02_{sample}_spades_rnaviral_contigs.fasta"
    log:
        logO="logs/spades/{sample}.log",
        logE="logs/spades/{sample}.err.log"
    conda:
        "../envs/spades_env.yaml"
    threads: 30
    shell:
        r"""
        rnaviralspades.py -t {threads} \
          -1 {input.r1} -2 {input.r2} \
          -o {output.d} \
          > {log.logO} 2> {log.logE}

        cp {output.d}/contigs.fasta {output.f}
        """