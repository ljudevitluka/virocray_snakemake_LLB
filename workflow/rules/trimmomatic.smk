rule trimmomatic:
    input:
        r1="resources/{sample}_1.fq.gz",
        r2="resources/{sample}_2.fq.gz",
        iclip="resources/TruSeq3-PE.fa",
    output:
        tp1="results/{sample}/01_{sample}_trim_1_paired.fq.gz",
        tp2="results/{sample}/01_{sample}_trim_2_paired.fq.gz",
        tup1="results/{sample}/01_{sample}_trim_1_unpaired.fq.gz",
        tup2="results/{sample}/01_{sample}_trim_2_unpaired.fq.gz",
        check="results/{sample}/01_{sample}_trim_pe.done",
    log:
        logO="logs/trim_pe/{sample}.log",
        logE="logs/trim_pe/{sample}.err.log",
    conda:
        "../envs/trimmomatic_env.yaml"
    threads: 4
    shell:
        """
        trimmomatic PE -threads {threads} -phred33 {input.r1} {input.r2} {output.tp1} {output.tup1} {output.tp2} {output.tup2} ILLUMINACLIP:{input.iclip}:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 > {log.logO} 2> {log.logE}
        touch {output.check}
        """