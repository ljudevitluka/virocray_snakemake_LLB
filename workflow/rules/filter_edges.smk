rule filter_edges:
    input:
        fasta = RESULTS_DIR + "/merged/merged_viral_contigs.fasta",
        blast = RESULTS_DIR + "/merged/blastn_pairs_sorted_no_self.tsv.gz"
    output:
        pairs = RESULTS_DIR + "/merged/blastn_pairs_sorted_no_self_filtered.tsv",
        edges = RESULTS_DIR + "/merged/edges_filtered.tsv"
    params:
        min_ani = 95.0,
        filter_mode = "AF",
        min_hsp_len = 150,
        min_af = 0.15
    log:
        logO="logs/filter_edges/filter_edges.log",
        logE="logs/filter_edges/filter_edges.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        python workflow/scripts/filter_edges.py \
            --blast_sorted_gz {input.blast} \
            --fasta {input.fasta} \
            --out_pairs {output.pairs} \
            --out_edges {output.edges} \
            --min_ani {params.min_ani} \
            --filter_mode {params.filter_mode} \
            --min_hsp_len {params.min_hsp_len} \
            --min_af {params.min_af} \
            > {log.logO} 2> {log.logE}
        """
