#!/usr/bin/env python3

import argparse
import os
import re
import dnaio


def extract_sequence(assembly_path, start, end):
    start_pattern = re.compile(start, re.IGNORECASE)
    end_pattern = re.compile(end, re.IGNORECASE)

    table = {}

    with dnaio.open(assembly_path) as reader:
        for record in reader:
            sequence = record.sequence
            name = record.name
            start_match = start_pattern.search(sequence)
            end_match = end_pattern.search(sequence)
            if start_match and end_match:
                table[name] = sequence[start_match.end():end_match.start()]
    return table


def is_valid_flycode(seq):
    pattern = r"^GGTAGT(GCA|GTT|GAT|CCA|GAA|ACT|GGT|TCT|TAC|CTG|TGG|CAG|TTC|AAC){6,8}(TGGCGG|TGGCTGCGG|TGGCAGTCTCGG|TGGCAGGAAGGAGGTCGG)$"
    return re.match(pattern, seq) is not None


def main(assembly, reference):
    gene_start = "GAGGAATTAACC"
    fc_end = "GGCCAAGGGGGT"
    sequences = extract_sequence(assembly, gene_start, fc_end)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly", type=str)
    parser.add_argument("--reference", type=str)

    args = parser.parse_args()
    main(args.assembly, args.reference)
