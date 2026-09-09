#!/usr/bin/env python3

import sys
from Bio import SeqIO
import csv


fasta_file = sys.argv[1]
blast_file = sys.argv[2]
output_file = sys.argv[3]


blast_columns = [
    "qseqid",
    "sacc",
    "stitle",
    "sscinames",
    "staxids",
    "bitscore",
    "evalue",
    "pident",
    "length",
    "qcovs",
    "qstart",
    "qend",
    "sstart",
    "send",
    "slen"
]


# read fasta IDs
contigs = []

for record in SeqIO.parse(fasta_file, "fasta"):
    contigs.append(record.id)


# read BLAST hits
blast_hits = {}

with open(blast_file) as f:
    reader = csv.DictReader(
        f,
        fieldnames=blast_columns,
        delimiter="\t"
    )

    for row in reader:
        blast_hits[row["qseqid"]] = row


# write final table

with open(output_file, "w") as out:

    writer = csv.writer(out, delimiter="\t")

    header = [
        "qseqid",
        "blastx_status",
        "sacc",
        "stitle",
        "sscinames",
        "staxids",
        "bitscore",
        "evalue",
        "pident",
        "length",
        "qcovs",
        "qstart",
        "qend",
        "sstart",
        "send",
        "slen"
    ]

    writer.writerow(header)


    for contig in contigs:

        if contig in blast_hits:

            row = blast_hits[contig]

            writer.writerow([
                contig,
                "hit",
                row["sacc"],
                row["stitle"],
                row["sscinames"],
                row["staxids"],
                row["bitscore"],
                row["evalue"],
                row["pident"],
                row["length"],
                row["qcovs"],
                row["qstart"],
                row["qend"],
                row["sstart"],
                row["send"],
                row["slen"]
            ])

        else:

            writer.writerow([
                contig,
                "no_hit",
                "NA",
                "NA",
                "NA",
                "NA",
                "0",
                "NA",
                "0",
                "0",
                "0",
                "NA",
                "NA",
                "NA",
                "NA",
                "NA"
            ])