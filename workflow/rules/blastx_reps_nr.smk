rule blastx_reps_nr:
    input:
        fasta=RESULTS_DIR + "/merged/cluster_representatives.fasta"
    output:
        tsv=RESULTS_DIR + "/merged/cluster_representatives_blastx_reps_nr.raw.tsv"
    params:
        db=config["diamond_nr_db"],  # BLAST nr; DIAMOND reads it directly (no .dmnd needed)
        taxdump=config["diamond_taxdump"]
    log:
        stdout="logs/blastx_reps_nr.stdout.log",
        stderr="logs/blastx_reps_nr.stderr.log"
    conda:
        "../envs/diamond_env.yaml"
    threads: 60
    shell:
        # Column order matches fill_missing_blastx.py; DIAMOND has no sacc/qcovs,
        # so sseqid (= accession for nr) and qcovhsp (= qcovs with --max-hsps 1) are used instead.
        """
        diamond blastx \
            --query {input.fasta} \
            --db {params.db} \
            --taxdump {params.taxdump} \
            --out {output.tsv} \
            --outfmt 6 qseqid sseqid stitle sscinames staxids bitscore evalue pident length qcovhsp qstart qend sstart send slen \
            --evalue 1e-5 \
            --max-target-seqs 1 \
            --max-hsps 1 \
            --sensitive \
            --threads {threads} \
            > {log.stdout} 2> {log.stderr}
        """
