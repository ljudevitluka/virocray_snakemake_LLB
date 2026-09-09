#!/usr/bin/env python3

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True, help='Input FASTA file')
parser.add_argument('--output', required=True, help='Output FASTA file with renamed headers')
parser.add_argument('--prefix', required=True, help='Prefix to prepend to each contig name')
args = parser.parse_args()

with open(args.input) as infile, open(args.output, 'w') as outfile:
    for record in SeqIO.parse(infile, 'fasta'):
        original_id = record.id
        record.id = f"{args.prefix}_{record.id}"
        record.description = ""  # remove description to avoid duplication
        SeqIO.write(record, outfile, 'fasta')
