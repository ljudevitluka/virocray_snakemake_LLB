rule blastx_nr:
    input:
        f = RESULTS_DIR + "/{sample}/04_{sample}_viral_contigs.fasta"
    output:
        t = RESULTS_DIR + "/{sample}/05_{sample}_viral_contigs_blastx.tsv"
    log:
        logO = "logs/blastx/{sample}.blastx.log",
        logE = "logs/blastx/{sample}.blastx.err.log"
    conda:
        "../envs/blast_env.yaml"
    threads: 60
    shell:
        """
        blastx \
            -query {input.f} \
            -db /biodbs/blastdb/120.ncbi/nr \
            -out {output.t} \
            -outfmt "6 qseqid sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids sscinames" \
            -max_target_seqs 10 \
            -evalue 1e-5 \
            -num_threads {threads} \
            > {log.logO} 2> {log.logE}
        """
