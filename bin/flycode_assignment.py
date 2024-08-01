#!/usr/bin/env python3

import argparse
import os
import re
import dnaio
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq


def extract_sequence(reference, start, end):
    start_pattern = re.compile(start, re.IGNORECASE)
    end_pattern = re.compile(end, re.IGNORECASE)

    with dnaio.open(reference, mode='r') as reader:
        for record in reader:
            sequence = record.sequence
            start_match = start_pattern.search(sequence)
            end_match = end_pattern.search(sequence)

            if start_match and end_match:
                return sequence[start_match.end():end_match.start()]
            break
    return None


def extract_sequences(assembly, start, end):
    start_pattern = re.compile(start, re.IGNORECASE)
    end_pattern = re.compile(end, re.IGNORECASE)

    table = {}

    with dnaio.open(assembly) as reader:
        for record in reader:
            sequence = record.sequence
            name = record.name
            start_match = start_pattern.search(sequence)
            end_match = end_pattern.search(sequence)
            if start_match and end_match:
                table[name] = sequence[start_match.end():end_match.start()]
    return table


def extract_flycode(sequence, start, end):
    start_match = re.search(start, sequence, re.IGNORECASE)
    end_match = re.search(end, sequence, re.IGNORECASE)
    return sequence[start_match.end():end_match.start()] if start_match and end_match else None


def is_valid_flycode_DNA(seq):
    pattern = r"^GGTAGT(GCA|GTT|GAT|CCA|GAA|ACT|GGT|TCT|TAC|CTG|TGG|CAG|TTC|AAC){6,8}(TGGCGG|TGGCTGCGG|TGGCAGTCTCGG|TGGCAGGAAGGAGGTCGG)$"
    return re.match(pattern, seq) is not None


def is_valid_flycode_AA(seq):
    pattern = r"^GS([ASTNQDEVLFYWGP]{6,8})(WR|WLR|WQSR|WLTVR|WQEGGR)$"
    return re.match(pattern, seq) is not None


def get_mutations(id, reference, sequence):
    mutations = []

    min_length = min(len(reference), len(sequence))
    for position in range(min_length):
        if reference[position] != sequence[position]:
            mutation = f'{reference[position]}{position + 1}{sequence[position]}'
            mutations.append(mutation)

    if (reference[:min_length] == sequence[:min_length] and not mutations) or reference == sequence:
        mutations.append('wt')
    return mutations


def process_tm287_288_fc(assembly, reference):
    orf1_start, orf1_end = "gaggaattaacc", "gatcagaagaag"
    orf2_start, orf2_end = "gggtgatgaacg", "TAAgaacaaaaa"
    tm288_start, tm288_end = "gggtgatgaacg", "GCAAGCTCCAGT"

    tm287_ref = extract_sequence(reference, orf1_start, orf1_end)
    tm288_ref = extract_sequence(reference, tm288_start, tm288_end)
    tm287_ref_aa = str(Seq(tm287_ref).translate())
    tm288_ref_aa = str(Seq(tm288_ref).translate())

    fc_start, fc_end = "LEVLFQGPSR", "GQGGHHHHHH*"

    orf1_dict = extract_sequences(assembly, orf1_start, orf1_end)
    orf2_dict = extract_sequences(assembly, orf2_start, orf2_end)

    flycode_assignment = defaultdict(list)

    for id, orf1 in orf1_dict.items():
        orf2 = orf2_dict.get(id)
        if orf2 and len(orf1) % 3 == 0 and len(orf2) % 3 == 0:
            orf1_aa = str(Seq(orf1).translate())
            orf2_aa = str(Seq(orf2).translate())
            flycode = extract_flycode(orf2_aa, fc_start, fc_end)
            if is_valid_flycode_AA(flycode):
                mutations_287 = get_mutations(id, tm287_ref_aa, orf1_aa,)
                mutations_288 = get_mutations(id, tm288_ref_aa, orf2_aa)
                mutations = f"NL|TM287={' '.join(mutations_287)};TM288={' '.join(mutations_288)}"
                flycode_assignment[flycode].append(mutations)
            else:
                print(id, "no flyocde")

    return (flycode_assignment)


def write_flycode_variant_fasta(file_name, variants):
    variant_dict = defaultdict(list)
    for variant, flycode in variants.items():
        if len(set(flycode)) == 1: # are all variants of a flycode identical
            variant_dict[flycode[0]].append(variant)
        else:
            print(f"flycode {variant}: variants {' '.join(flycode)}")

    with open(file_name, "w") as file:
        for variant, flycode in variant_dict.items():
            file.write(f">{variant}\n")
            file.write(f"{''.join(flycode)}\n")


def process_protein_of_interest(case_name, assembly, reference):
    match case_name:
        case "nanobody_FC":
            print("not implemented")
        case "KDELR_FC_GFP":
            print("not implemented")
        case "TM287/288_FC":
            return process_tm287_288_fc(assembly, reference)


def main(poi, assembly, file_name, reference):
    flycode_assignment = process_protein_of_interest(poi, assembly, reference)
    write_flycode_variant_fasta(file_name, flycode_assignment)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--poi", type=str)
    parser.add_argument("--assembly", type=str)
    parser.add_argument("--file_name", type=str)
    parser.add_argument("--reference", type=str)

    args = parser.parse_args()
    main(args.poi, args.assembly, args.file_name, args.reference)
