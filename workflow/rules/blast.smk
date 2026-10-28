rule blast:
    input:
        f = RESULTS_DIR + "/merged/merged_viral_contigs.fasta"
    output:
        db = temp(directory(RESULTS_DIR + "/merged/blastdb")),
        t = RESULTS_DIR + "/merged/blastn_pairs.tsv"
    log:
        logO = "logs/blast/blast.log",
        logE = "logs/blast/blast.err.log"
    conda:
        "../envs/blast_env.yaml"
    threads: 40
    resources:
        mem_mb = 20000
    shell:
        """
        mkdir -p {output.db}
        makeblastdb -in {input.f} -dbtype nucl \
            -out {output.db}/blastdb \
            > {log.logO} 2> {log.logE}
        
        blastn -task blastn \
            -word_size 28 \
            -query {input.f} \
            -db {output.db}/blastdb \
            -out {output.t} \
            -outfmt "6 qseqid sseqid pident length qstart qend sstart send" \
            -num_threads {threads} \
            >> {log.logO} 2>> {log.logE}
        """