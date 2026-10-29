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

3. MMseqs2 taxonomy can run per sample or for all samples combined (see `mmseqs2_mode` below). This is a major speed up for the whole pipeline. 

4. Reads are mapped to the cluster representatives with BBMap `mappedonly=t`. This reduced the storage memory usage of the pipeline. 

5. Functional annotation of all cluster representatives: six-frame ORF calling with EMBOSS `getorf`, then InterProScan 5.78 and NeoRdRp2 `hmmsearch`. A contig gets `has_RdRp = yes` when either program detects an RdRp; E-values, scores and HMM coverage are reported as confidence values. Outputs:

```text
results/merged/functional_annotation/reps_functional_annotation.tsv
results/merged/functional_annotation/reps_orfs_interproscan.tsv
results/merged/functional_annotation/reps_rdrp_proteins.faa
```

6. Final table combining TPM, taxonomy, blastx and functional annotation:

```text
results/merged/final_table_tpm_taxonomy_blastx.tsv
```

7. Memory limit: jobs are limited to 500 GB RAM and the whole run to a hard 900 GB by `run_snakemake.sh` (see `MAX_MEM_GB` and `HARD_MEM_GB` below).

8. Hardcoded values were moved to `config/config.yaml`.

## Options

### Launcher (`run_snakemake.sh`)

| Option | Default | Description |
|---|---|---|
| `MAX_MEM_GB` | `500` | RAM budget for the jobs in GB. Snakemake only starts jobs that fit in it. |
| `HARD_MEM_GB` | `900` | Hard memory limit for the whole run in GB, including cached database files. The run is killed if it goes over. |
| `CONDA_EXE` | `/proj/etools/conda/bin/conda` | Conda executable used to run Snakemake. |

### `config/config.yaml`

| Option | Default | Description |
|---|---|---|
| `mmseqs2_db` | `.../path_to_database/nr.fnaDB` | MMseqs2 database for contig taxonomy |
| `mmseqs2_split_memory_limit` | `200G` | Memory per MMseqs2 job (lower = less memory/slower) |
| `mmseqs2_combined_split_memory_limit` | `400G` | Memory for the MMseqs2 job in `combined` mode |
| `mmseqs2_mode` | `combined` | `per_sample`: one MMseqs2 job per sample. `combined`: all samples in one job. Both write the same `03_{sample}_mmseqs2_taxonomy.tsv` files. In `combined` mode, a crashed MMseqs2 job resumes where it stopped when the run is restarted with the same inputs, threads and memory limit. |
| `spades_memory_gb` | `120` | SPAdes memory limit (`-m`) |
| `diamond_nr_db` | `/path_to_database/nr` | Database for DIAMOND blastx |
| `diamond_taxdump` | `/path_to_database/taxdump` | NCBI taxdump for DIAMOND |
| `min_contig_length` | `500` | Minimum contig length after assembly (nt) |
| `filter_edges.min_ani` | `90.0` | Minimum ANI between two contigs |
| `filter_edges.filter_mode` | `AF` | `AF` (aligned fraction) or `HSP` (longest HSP length) |
| `filter_edges.min_hsp_len` | `150` | Used when `filter_mode: HSP` |
| `filter_edges.min_af` | `0.80` | Used when `filter_mode: AF` |
| `clustering.resolution` | `1.0` | Leiden resolution |
| `bbmap.minid` / `bbmap.maxindel` | `0.90` / `3` | BBMap mapping settings |
| `coverm.min_read_percent_identity` / `coverm.min_read_aligned_percent` | `95` / `95` | CoverM read filters |
| `filter_read_counts.min_reads` | `2` | Minimum reads per representative |
| `filter_read_counts.min_covered_fraction` | `0.25` | Minimum covered fraction |
| `filter_read_counts.min_covered_bases` | `400` | Minimum covered bases |
| `getorf.min_orf_nt` / `getorf.table` | `300` / `1` | Minimum ORF length (nt) and genetic code |
| `annotation_chunks` | `20` | Number of chunks for InterProScan / hmmsearch |
| `interproscan.path` | `.../interproscan-5.78-109.0/interproscan.sh` | InterProScan installation |
| `interproscan.applications` | `Pfam,CDD,SUPERFAMILY,Gene3D,ProSiteProfiles,NCBIfam` | InterProScan member databases |
| `neordrp.hmm` / `neordrp.evalue` | `/path_to_database/NeoRdRp.2.1.hmm` / `1e-10` | NeoRdRp2 HMMs and E-value cutoff |
| `rdrp.interpro_ids` | `config/rdrp_interpro_ids.tsv` | InterPro entries counted as RdRp |
| `rdrp.signature_regex` | see config | Signature descriptions counted as RdRp |
| `rdrp.exclude_signatures` | `[PF05183]` | Signatures never counted as RdRp |

Any option can also be set for a single run without editing the file:

```bash
./run_snakemake.sh --cores 60 -k -p --config mmseqs2_mode=combined
```


## Troubleshooting

Inspect the rule-specific error log when a job fails. Logs are stored under `logs/`, with a separate directory for each rule.

To check the workflow without running jobs:

```bash
./run_snakemake.sh --cores 1 --dry-run
```

