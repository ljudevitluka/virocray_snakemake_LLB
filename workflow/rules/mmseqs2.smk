rule mmseqs2:
    input:
        f = RESULTS_DIR + "/{sample}/02_{sample}_spades_rnaviral_contigs_renamed_min500.fasta"
    output:
        querydb_idx  = RESULTS_DIR + "/{sample}/03_{sample}_queryDB.index",
        resultdb_idx = RESULTS_DIR + "/{sample}/03_{sample}_resultDB.index",
        tmpdir       = temp(directory(RESULTS_DIR + "/{sample}/03_{sample}_tmp"))
    params:
        queryprefix  = RESULTS_DIR + "/{sample}/03_{sample}_queryDB",
        resultprefix = RESULTS_DIR + "/{sample}/03_{sample}_resultDB"
    log:
        logO = "logs/mmseqs2/{sample}.db.log",
        logE = "logs/mmseqs2/{sample}.db.err.log"
    conda:
        "../envs/mmseqs2_env.yaml"
    threads: 30
    shell:
        r"""
        mmseqs createdb {input.f} {params.queryprefix} \
          > {log.logO} 2> {log.logE}

        mmseqs taxonomy \
          {params.queryprefix} \
          /biodbs/mmseqs2/nr_database/nr.fnaDB \
          {params.resultprefix} \
          {output.tmpdir} \
          --threads {threads} --tax-lineage 1 \
          >> {log.logO} 2>> {log.logE}
        """
