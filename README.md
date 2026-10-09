# ViroCray Snakemake pipeline

## Running the pipeline

Run the pipeline from the project root with the included launcher:

```bash
cd /path/to/root/virocray_snakemake_LLB
./run_snakemake.sh --cores 60 -k -p
```

The launcher automatically enables Snakemake's per-rule Conda environments. Any additional Snakemake options can be passed after the script name. For example, to perform a dry run:

```bash
./run_snakemake.sh --cores 1 --dry-run
```

## Why use the launcher?

The pipeline uses Conda environments for individual Snakemake rules. Previously, a shell in `screen` could report that the `snakemake` environment was active while `python` still resolved to the base interpreter:

```text
/proj/etools/conda/bin/python
```

That interpreter did not contain packages required by the workflow, such as `pandas`, causing Python rules to fail with `ModuleNotFoundError`.

`run_snakemake.sh` avoids relying on shell activation and runs Snakemake explicitly through:

```bash
conda run -n snakemake
```

This makes the workflow use the intended Snakemake environment consistently in normal terminals and in `screen` sessions. Snakemake then activates the environments declared by each rule through `--use-conda`.

## Requirements

The configured Conda installation must be available at `/proj/etools/conda/bin/conda`, and the `snakemake` environment must exist there. The launcher also accepts a `CONDA_EXE` override:

```bash
CONDA_EXE=/path/to/conda ./run_snakemake.sh --cores 1 --dry-run
```

The workflow configuration and sample list are defined in `config/config.yaml` and the file referenced by its `samples` setting.

## Pipeline changes

1. The `filter_edges` rule now uses a minimum average nucleotide identity (`min_ani`) of `90.0`, reduced from `95.0` and AF (`min_af`) of `0.80`, reduced from `85`.

2. Clustering QC visualisation 

After the Leiden clustering step, the workflow also generates a clustering quality figure at:

```text
results/merged/clustering_qc/clustering_qc.pdf
```

The plot summarises three things:

- clustered vs singleton contigs,
- contig degree distribution (number of edges per contig),
- cluster-size distribution across multi-member clusters.


## Troubleshooting

Inspect the rule-specific error log when a job fails. Logs are stored under `logs/`, with a separate directory for each rule.

To check the workflow without running jobs:

```bash
./run_snakemake.sh --cores 1 --dry-run
```

