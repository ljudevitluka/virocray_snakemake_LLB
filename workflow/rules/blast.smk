rule blast:
    input:
        f = RESULTS_DIR + "/merged/merged_viral_contigs.fasta"
    output:
        db = temp(RESULTS_DIR + "/merged/blastdb.nhr"),
        t = RESULTS_DIR + "/merged/blastn_pairs.tsv"
    log:
        logO = "logs/blast/blast.log",
        logE = "logs/blast/blast.err.log"
    conda:
        "../envs/blast_env.yaml"
    threads: 40
    shell:
        """
        makeblastdb -in {input.f} -dbtype nucl \
            -out {RESULTS_DIR}/merged/blastdb \
            > {log.logO} 2> {log.logE}
        
        blastn -query {input.f} \
            -db {RESULTS_DIR}/merged/blastdb \
            -out {output.t} \
            -outfmt "6 qseqid sseqid pident length qstart qend sstart send" \
            -num_threads {threads} \
            >> {log.logO} 2>> {log.logE}
        """