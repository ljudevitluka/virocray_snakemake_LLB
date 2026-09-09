#!/usr/bin/env python3

import argparse
import pandas as pd
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(
        description="Extract contigs classified as Viruses from FASTA based on mmseqs2 taxonomy TSV."
    )
    parser.add_argument("--fasta", required=True, help="Input FASTA file with renamed contigs")
    parser.add_argument("--taxonomy", required=True, help="Taxonomy TSV file from mmseqs2")
    parser.add_argument("--output_fasta", required=True, help="Output FASTA file with only viral contigs")
    parser.add_argument("--output_tsv", required=True, help="Output TSV file with only viral contig rows")
    args = parser.parse_args()

    df = pd.read_csv(args.taxonomy, sep="\t", header=None, dtype=str)

    if df.shape[1] < 9:
        # Nothing to do (or you can raise)
        df_viral = df.iloc[0:0].copy()
        viral_ids = set()
    else:
        df = df.iloc[:, :9].copy()
        df.columns = ["contig", "tax_id", "rank", "name", "retained", "assigned", "label_match", "support", "lineage"]

        # keep only rows where lineage exists
        df = df[df["lineage"].notna()].copy()

        # domain from lineage; clean leading -_ tokens
        df["domain"] = (
            df["lineage"]
            .str.split(";", expand=False)
            .str[0]
            .str.replace(r"^[-_]+", "", regex=True)
        )

        df_viral = df[df["domain"] == "Viruses"].copy()
        viral_ids = set(df_viral["contig"].tolist())

    print(f"Total rows in taxonomy TSV: {len(df)}")
    print(f"Total contigs classified as 'Viruses': {len(viral_ids)}")

    # write filtered TSV
    df_viral.to_csv(args.output_tsv, sep="\t", header=False, index=False)

    # write viral FASTA
    count_written = 0
    with open(args.output_fasta, "w") as out_fasta:
        for record in SeqIO.parse(args.fasta, "fasta"):
            if record.id in viral_ids:
                SeqIO.write(record, out_fasta, "fasta")
                count_written += 1

    print(f"Total contigs written to output FASTA: {count_written}")

if __name__ == "__main__":
    main()
