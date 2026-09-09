rule bbmap_to_reps:
    input:
        r1  = RESULTS_DIR + "/{sample}/01_{sample}_trim_1_paired.fq.gz",
        r2  = RESULTS_DIR + "/{sample}/01_{sample}_trim_2_paired.fq.gz",
        ref = RESULTS_DIR + "/merged/cluster_representatives.fasta"
    output:
        bam   = RESULTS_DIR + "/{sample}/07_{sample}_bbmap_reps.bam",
        bai   = RESULTS_DIR + "/{sample}/07_{sample}_bbmap_reps.bam.bai",
        check = RESULTS_DIR + "/{sample}/07_{sample}_bbmap_reps.done"
    log:
        logE = "logs/bbmap_to_reps/{sample}.err.log"
    conda:
        "../envs/mapping_env.yaml"
    threads: 20
    shell:
        r"""
        bbmap.sh \
            nodisk=t \
            in={input.r1} in2={input.r2} \
            ref={input.ref} \
            ambig=random \
            minid=0.90 \
            maxindel=3 \
            threads={threads} \
            out={output.bam}.sam \
            2> {log.logE}

        samtools view -b -@ {threads} {output.bam}.sam \
        | samtools sort -@ {threads} -o {output.bam}

        samtools index -@ {threads} {output.bam} 2>> {log.logE}

        rm {output.bam}.sam
        touch {output.check}
        """
