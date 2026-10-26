ANNOT_DIR = RESULTS_DIR + "/merged/functional_annotation"
ORF_CHUNKS = [f"{i:03d}" for i in range(1, config["annotation_chunks"] + 1)]


rule getorf_reps:
    input:
        fasta = RESULTS_DIR + "/merged/cluster_representatives.fasta"
    output:
        orfs = ANNOT_DIR + "/reps_orfs.faa"
    params:
        min_orf_nt = config["getorf"]["min_orf_nt"],
        table      = config["getorf"]["table"]
    log:
        logO = "logs/functional_annotation/getorf.log",
        logE = "logs/functional_annotation/getorf.err.log"
    conda:
        "../envs/emboss_env.yaml"
    shell:
        # -find 1: translate regions between stop codons in all six frames, so
        # partial genes at contig ends (no start codon) are kept.
        # InterProScan rejects '*', so any stop symbols are removed.
        r"""
        getorf -sequence {input.fasta} -outseq {output.orfs}.tmp \
            -find 1 -minsize {params.min_orf_nt} -table {params.table} \
            > {log.logO} 2> {log.logE}
        sed '/^>/! s/\*//g' {output.orfs}.tmp > {output.orfs}
        rm {output.orfs}.tmp
        """


rule split_orfs:
    input:
        orfs = ANNOT_DIR + "/reps_orfs.faa"
    output:
        chunks = expand(ANNOT_DIR + "/chunks/reps_orfs.part_{chunk}.faa", chunk=ORF_CHUNKS)
    params:
        outdir = ANNOT_DIR + "/chunks",
        n      = config["annotation_chunks"]
    log:
        logE = "logs/functional_annotation/split_orfs.err.log"
    conda:
        "../envs/emboss_env.yaml"
    shell:
        """
        seqkit split2 -p {params.n} -O {params.outdir} --force {input.orfs} 2> {log.logE}
        """


rule interproscan_chunk:
    input:
        chunk = ANNOT_DIR + "/chunks/reps_orfs.part_{chunk}.faa"
    output:
        tsv = ANNOT_DIR + "/interproscan/reps_orfs.part_{chunk}.tsv"
    params:
        ips  = config["interproscan"]["path"],
        appl = config["interproscan"]["applications"],
        tmp  = ANNOT_DIR + "/interproscan/tmp_{chunk}"
    log:
        logO = "logs/functional_annotation/interproscan_{chunk}.log",
        logE = "logs/functional_annotation/interproscan_{chunk}.err.log"
    threads: 8
    shell:
        # -dp: no lookup service, match calculation is done locally
        """
        {params.ips} -i {input.chunk} -f TSV -o {output.tsv} \
            -appl {params.appl} -dp -iprlookup \
            -cpu {threads} -T {params.tmp} \
            > {log.logO} 2> {log.logE}
        rm -rf {params.tmp}
        """


rule neordrp_hmmsearch_chunk:
    input:
        chunk = ANNOT_DIR + "/chunks/reps_orfs.part_{chunk}.faa"
    output:
        domtbl = ANNOT_DIR + "/neordrp/reps_orfs.part_{chunk}.domtbl"
    params:
        hmm    = config["neordrp"]["hmm"],
        evalue = config["neordrp"]["evalue"]
    log:
        logE = "logs/functional_annotation/neordrp_{chunk}.err.log"
    conda:
        "../envs/hmmer_env.yaml"
    threads: 4
    shell:
        """
        hmmsearch --cpu {threads} -E {params.evalue} --domE {params.evalue} \
            --domtblout {output.domtbl} --noali -o /dev/null \
            {params.hmm} {input.chunk} \
            2> {log.logE}
        """


rule summarize_functional_annotation:
    input:
        reps     = RESULTS_DIR + "/merged/cluster_representatives.fasta",
        orfs     = ANNOT_DIR + "/reps_orfs.faa",
        ips      = expand(ANNOT_DIR + "/interproscan/reps_orfs.part_{chunk}.tsv", chunk=ORF_CHUNKS),
        neordrp  = expand(ANNOT_DIR + "/neordrp/reps_orfs.part_{chunk}.domtbl", chunk=ORF_CHUNKS),
        rdrp_ids = config["rdrp"]["interpro_ids"]
    output:
        table    = ANNOT_DIR + "/reps_functional_annotation.tsv",
        ips      = ANNOT_DIR + "/reps_orfs_interproscan.tsv",
        rdrp_faa = ANNOT_DIR + "/reps_rdrp_proteins.faa"
    params:
        regex   = config["rdrp"]["signature_regex"],
        exclude = " ".join(config["rdrp"]["exclude_signatures"]),
        evalue  = config["neordrp"]["evalue"]
    log:
        logO = "logs/functional_annotation/summarize.log",
        logE = "logs/functional_annotation/summarize.err.log"
    conda:
        "../envs/core_env.yaml"
    shell:
        """
        cat {input.ips} > {output.ips}
        python workflow/scripts/summarize_functional_annotation.py \
            --reps {input.reps} \
            --orfs {input.orfs} \
            --interproscan {output.ips} \
            --neordrp {input.neordrp} \
            --rdrp-interpro-ids {input.rdrp_ids} \
            --rdrp-signature-regex "{params.regex}" \
            --exclude-signatures {params.exclude} \
            --neordrp-evalue {params.evalue} \
            --out-table {output.table} \
            --out-rdrp-faa {output.rdrp_faa} \
            > {log.logO} 2> {log.logE}
        """
