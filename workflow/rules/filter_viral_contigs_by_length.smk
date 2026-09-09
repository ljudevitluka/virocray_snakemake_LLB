rule filter_viral_contigs_by_length:
    input:
        f = RESULTS_DIR + "/{sample}/02_{sample}_spades_rnaviral_contigs_renamed.fasta"
    output:
        f = RESULTS_DIR + "/{sample}/02_{sample}_spades_rnaviral_contigs_renamed_min500.fasta"
    log:
        logO = "logs/filter_length/{sample}.log",
        logE = "logs/filter_length/{sample}.err.log"
    conda:
        "../envs/seqkit.yaml"
    threads: 1
    shell:
        r"""
        set -euo pipefail
        seqkit seq -m 500 {input.f} > {output.f} 2> {log.logE}
        """
