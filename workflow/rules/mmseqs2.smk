rule mmseqs2:
    input:
        f = RESULTS_DIR + "/{sample}/02_{sample}_spades_contigs_renamed_min500.fasta"
    output:
        querydb_idx  = RESULTS_DIR + "/{sample}/03_{sample}_queryDB.index",
        resultdb_idx = RESULTS_DIR + "/{sample}/03_{sample}_resultDB.index",
        tmpdir       = temp(directory(RESULTS_DIR + "/{sample}/03_{sample}_tmp"))
    params:
        queryprefix  = RESULTS_DIR + "/{sample}/03_{sample}_queryDB",
        resultprefix = RESULTS_DIR + "/{sample}/03_{sample}_resultDB",
        db           = config["mmseqs2_db"]
    log:
        logO = "logs/mmseqs2/{sample}.db.log",
        logE = "logs/mmseqs2/{sample}.db.err.log"
    conda:
        "../envs/mmseqs2_env.yaml"
    threads: 30
    shell:
        # mmseqs refuses to overwrite existing databases, and Snakemake only removes
        # the declared .index files before a rerun, so clear the old DB files first
        r"""
        rm -f {params.queryprefix} {params.queryprefix}.* {params.queryprefix}_h* \
              {params.resultprefix} {params.resultprefix}.* {params.resultprefix}_*

        mmseqs createdb {input.f} {params.queryprefix} \
          > {log.logO} 2> {log.logE}

        mmseqs taxonomy \
          {params.queryprefix} \
          {params.db} \
          {params.resultprefix} \
          {output.tmpdir} \
          --threads {threads} --tax-lineage 1 \
          >> {log.logO} 2>> {log.logE}
        """
