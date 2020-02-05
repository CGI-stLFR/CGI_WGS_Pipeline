#!/usr/bin/env python3
import sys

if len(sys.argv) > 1:
    fasta = sys.argv[1]
else:
    print("Please supply a fasta file.", file=sys.stderr)
    sys.exit(1)

with open(fasta, "r") as fa:
    length = 0
    for line in fa:
        if line.startswith(">"):
            continue
        length += len(line.strip().replace('N', '').replace('n', ''))
    print(length)
