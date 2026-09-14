#!/usr/bin/env python3

import argparse


def read_fasta(path):
    header = None
    sequence = []

    with open(path) as infile:
        for line in infile:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(sequence)
                header = line[1:].split()[0]
                sequence = []
            elif header is None:
                raise ValueError('FASTA sequence found before the first header')
            else:
                sequence.append(line)

    if header is not None:
        yield header, ''.join(sequence)

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True, help='Input FASTA file')
parser.add_argument('--output', required=True, help='Output FASTA file with renamed headers')
parser.add_argument('--prefix', required=True, help='Prefix to prepend to each contig name')
args = parser.parse_args()

with open(args.output, 'w') as outfile:
    for original_id, sequence in read_fasta(args.input):
        outfile.write(f">{args.prefix}_{original_id}\n")
        for position in range(0, len(sequence), 60):
            outfile.write(f"{sequence[position:position + 60]}\n")
