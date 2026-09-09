#!/usr/bin/env python3

import argparse
import gzip
import pysam

parser = argparse.ArgumentParser(description="Extract reads mapping to viral contigs.")
parser.add_argument("--bam", required=True, help="Input BAM file.")
parser.add_argument("--viral-fasta", required=True, help="FASTA containing viral contig IDs.")
parser.add_argument("--out-r1", required=True, help="Output FASTQ.gz for read1.")
parser.add_argument("--out-r2", required=True, help="Output FASTQ.gz for read2.")
args = parser.parse_args()

# --------------------
# Load viral contig IDs
# --------------------
viral_contigs = set()

with open(args.viral_fasta) as f:
    for line in f:
        if line.startswith(">"):
            contig = line[1:].strip().split()[0]
            viral_contigs.add(contig)

print(f"Loaded {len(viral_contigs)} viral contig IDs")

# --------------------
# Open BAM + output FASTQ files
# --------------------
bam = pysam.AlignmentFile(args.bam, "rb")
out_r1 = gzip.open(args.out_r1, "wt")
out_r2 = gzip.open(args.out_r2, "wt")

# --------------------
# Extract read qualities to FASTQ format
# --------------------
def write_fastq(handle, aln):
    seq = aln.query_sequence
    quals = "".join(chr(q + 33) for q in aln.query_qualities)
    handle.write(f"@{aln.query_name}\n{seq}\n+\n{quals}\n")

# --------------------
# Iterate BAM alignments
# --------------------
count = 0

for aln in bam.fetch(until_eof=True):
    if aln.is_unmapped:
        continue

    ref = bam.get_reference_name(aln.reference_id)

    if ref not in viral_contigs:
        continue

    count += 1

    if aln.is_read1:
        write_fastq(out_r1, aln)
    elif aln.is_read2:
        write_fastq(out_r2, aln)
    else:
        # single-end or missing flag treat as R1 by default
        write_fastq(out_r1, aln)

bam.close()
out_r1.close()
out_r2.close()

print(f"Done. Extracted {count} viral-mapped reads.")
