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
        tsv    = MMSEQS_COMBINED_DIR + "/all_samples_mmseqs2_taxonomy.tsv"
    params:
        queryprefix  = MMSEQS_COMBINED_DIR + "/queryDB",
        resultprefix = MMSEQS_COMBINED_DIR + "/resultDB",
        # Not an output, so Snakemake keeps it after a crash and MMseqs2 resumes from it
        tmpdir       = MMSEQS_COMBINED_DIR + "/tmp",
        db           = config["mmseqs2_db"],
        mem_limit    = config["mmseqs2_combined_split_memory_limit"]
    log:
        logO = "logs/mmseqs2/combined.log",
        logE = "logs/mmseqs2/combined.err.log"
    conda:
        "../envs/mmseqs2_env.yaml"
    # One job for all samples, so it gets all cores
    threads: 60
    resources:
        # split-memory-limit + 30 GB overhead
        mem_mb = int(config["mmseqs2_combined_split_memory_limit"].rstrip("Gg")) * 1000 + 30000
    shell:
        # After a crash, MMseqs2 resumes from its tmp folder if the input contigs, threads
        # and memory limit are unchanged (MMseqs2 names the folder after these and after
        # queryDB, so queryDB is then kept as is; rebuilding it makes MMseqs2 start over).
        # If any of them changed, the old tmp folder is removed and the search starts fresh.
        # mmseqs refuses to overwrite existing databases, so clear old DB files first
        r"""
        sig=$( (stat -c '%n %s %Y' {input.fasta}; echo "threads={threads} mem={params.mem_limit}") | md5sum | cut -d' ' -f1)
        rm -f {params.resultprefix} {params.resultprefix}.* {params.resultprefix}_*

        if [ -d {params.tmpdir} ] && [ -e {params.queryprefix}.dbtype ] && [ "$(cat {params.tmpdir}/inputs.md5 2>/dev/null)" = "$sig" ]; then
            echo "Resuming from {params.tmpdir}" >> {log.logO}
        else
            rm -rf {params.tmpdir}
            rm -f {params.queryprefix} {params.queryprefix}.* {params.queryprefix}_h*
            mkdir -p {params.tmpdir}
            echo "$sig" > {params.tmpdir}/inputs.md5

            mmseqs createdb {input.fasta} {params.queryprefix} \
              >> {log.logO} 2>> {log.logE}
        fi

        mmseqs taxonomy \
          {params.queryprefix} \
          {params.db} \
          {params.resultprefix} \
          {params.tmpdir} \
          --threads {threads} --tax-lineage 1 \
          --split-memory-limit {params.mem_limit} \
          >> {log.logO} 2>> {log.logE}

        mmseqs createtsv {params.queryprefix} {params.resultprefix} {output.tsv} \
          >> {log.logO} 2>> {log.logE}

        rm -rf {params.tmpdir}
        """


rule split_mmseqs2_taxonomy:
    input:
        tsv = MMSEQS_COMBINED_DIR + "/all_samples_mmseqs2_taxonomy.tsv"
    output:
        tsv = RESULTS_DIR + "/{sample}/03_{sample}_mmseqs2_taxonomy.tsv"
    shell:
        r"""
        awk -F'\t' -v p="{wildcards.sample}_NODE_" 'index($1, p) == 1' {input.tsv} > {output.tsv}
        """
