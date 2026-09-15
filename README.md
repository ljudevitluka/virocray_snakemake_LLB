# ViroCray Snakemake pipeline

## Running the pipeline

Run the pipeline from the project root with the included launcher:

```bash
cd /scratch/lukabostjancic/ViroCray/virocray_snakemake_LLB
./run_snakemake.sh --cores 60 -k -p
```

The launcher automatically enables Snakemake's per-rule Conda environments. Any additional Snakemake options can be passed after the script name. For example, to perform a dry run:

```bash
./run_snakemake.sh --cores 1 --dry-run
```

To run the workflow in a `screen` session:

```bash
screen -S virocray
cd /scratch/lukabostjancic/ViroCray/virocray_snakemake_LLB
./run_snakemake.sh --cores 60 -k -p
```

Detach from the session with `Ctrl-a` followed by `d`. Reconnect later with:

```bash
screen -r virocray
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

## Troubleshooting

Inspect the rule-specific error log when a job fails. Logs are stored under `logs/`, with a separate directory for each rule.

To check the workflow without running jobs:

```bash
./run_snakemake.sh --cores 1 --dry-run
```

A Conda warning about strict channel priority is separate from the Python import failure. For more reproducible Conda resolution, strict channel priority can be enabled with:

```bash
conda config --set channel_priority strict
```
