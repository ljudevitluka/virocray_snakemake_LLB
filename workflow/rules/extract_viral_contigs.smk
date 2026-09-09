rule extract_viral_contigs:
    input:
        fasta="results/{sample}/02_{sample}_spades_rnaviral_contigs_renamed.fasta",
        taxonomy="results/{sample}/03_{sample}_mmseqs2_taxonomy.tsv",
    output:
        viral_fasta="results/{sample}/04_{sample}_viral_contigs.fasta",
        viral_tsv="results/{sample}/04_{sample}_viral_taxonomy.tsv",
    log:
        logO="logs/extract_viral_contigs/{sample}.log",
        logE="logs/extract_viral_contigs/{sample}.err.log",
    conda:
        "../envs/extract_viral_contigs_env.yaml"
    shell:
        r"""
        set -euo pipefail
        python workflow/scripts/extract_viral_contigs.py \
          --fasta {input.fasta} \
          --taxonomy {input.taxonomy} \
          --output_fasta {output.viral_fasta} \
          --output_tsv {output.viral_tsv} \
          > {log.logO} 2> {log.logE}
        """