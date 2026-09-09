rule mmseqs2_tsv:
    input:
        querydb_idx  = RESULTS_DIR + "/{sample}/03_{sample}_queryDB.index",
        resultdb_idx = RESULTS_DIR + "/{sample}/03_{sample}_resultDB.index"
    params:
        qprefix = RESULTS_DIR + "/{sample}/03_{sample}_queryDB",
        rprefix = RESULTS_DIR + "/{sample}/03_{sample}_resultDB"
    output:
        tsv = RESULTS_DIR + "/{sample}/03_{sample}_mmseqs2_taxonomy.tsv"
    log:
        logO = "logs/mmseqs2/{sample}.tsv.log",
        logE = "logs/mmseqs2/{sample}.tsv.err.log"
    conda:
        "../envs/mmseqs2_env.yaml"
    threads: 4
    shell:
        r"""
        mmseqs createtsv \
          {params.qprefix} \
          {params.rprefix} \
          {output.tsv} \
          > {log.logO} 2> {log.logE}
        """
