# Used when mmseqs2_mode: combined. All samples are searched in one MMseqs2 run,
# so the nr index (built once per database split) is built once instead of once
# per sample. Memory is set by the nr index (--split-memory-limit), not by the
# number of query contigs. The result is split back into the same per-sample
# 03_{sample}_mmseqs2_taxonomy.tsv files the per-sample mode produces.
MMSEQS_COMBINED_DIR = RESULTS_DIR + "/merged/mmseqs2_combined"


rule mmseqs2_combined:
    input:
        fasta = expand(
            RESULTS_DIR + "/{sample}/02_{sample}_spades_contigs_renamed_min500.fasta",
            sample=samples["sample"].tolist()
        )
    output:
        tsv    = MMSEQS_COMBINED_DIR + "/all_samples_mmseqs2_taxonomy.tsv",
        tmpdir = temp(directory(MMSEQS_COMBINED_DIR + "/tmp"))
    params:
        queryprefix  = MMSEQS_COMBINED_DIR + "/queryDB",
        resultprefix = MMSEQS_COMBINED_DIR + "/resultDB",
        db           = config["mmseqs2_db"],
        mem_limit    = config["mmseqs2_split_memory_limit"]
    log:
        logO = "logs/mmseqs2/combined.log",
        logE = "logs/mmseqs2/combined.err.log"
    conda:
        "../envs/mmseqs2_env.yaml"
    # One job for all samples, so it gets all cores
    threads: 60
    resources:
        mem_mb = 230000  # split-memory-limit + overhead
    shell:
        # mmseqs refuses to overwrite existing databases, so clear old DB files first
        r"""
        rm -f {params.queryprefix} {params.queryprefix}.* {params.queryprefix}_h* \
              {params.resultprefix} {params.resultprefix}.* {params.resultprefix}_*

        mmseqs createdb {input.fasta} {params.queryprefix} \
          > {log.logO} 2> {log.logE}

        mmseqs taxonomy \
          {params.queryprefix} \
          {params.db} \
          {params.resultprefix} \
          {output.tmpdir} \
          --threads {threads} --tax-lineage 1 \
          --split-memory-limit {params.mem_limit} \
          >> {log.logO} 2>> {log.logE}

        mmseqs createtsv {params.queryprefix} {params.resultprefix} {output.tsv} \
          >> {log.logO} 2>> {log.logE}
        """


rule split_mmseqs2_taxonomy:
    input:
        tsv = MMSEQS_COMBINED_DIR + "/all_samples_mmseqs2_taxonomy.tsv"
    output:
        tsv = RESULTS_DIR + "/{sample}/03_{sample}_mmseqs2_taxonomy.tsv"
    shell:
        # Contig names are <sample>_NODE_..., so the prefix also keeps UL1 apart from UL10
        r"""
        awk -F'\t' -v p="{wildcards.sample}_NODE_" 'index($1, p) == 1' {input.tsv} > {output.tsv}
        """
