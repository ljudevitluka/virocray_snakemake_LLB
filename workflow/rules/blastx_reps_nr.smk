rule blastx_reps_nr:
    input:
        fasta=RESULTS_DIR + "/merged/cluster_representatives_filtered_length.fasta"
    output:
        tsv=RESULTS_DIR + "/merged/cluster_representatives_blastx_reps_nr.raw.tsv"
    log:
        stdout="logs/blastx_reps_nr.stdout.log",
        stderr="logs/blastx_reps_nr.stderr.log"
    conda:
        "../envs/blast_env.yaml"
    threads: 60
    shell:
        """
        blastx \
            -query {input.fasta} \
            -db /biodbs/blastdb/120.ncbi/nr \
            -out {output.tsv} \
            -outfmt "6 qseqid sacc stitle sscinames staxids bitscore evalue pident length qcovs qstart qend sstart send slen" \
            -evalue 1e-5 \
            -max_target_seqs 1 \
            -num_threads {threads} \
            > {log.stdout} 2> {log.stderr}
        """