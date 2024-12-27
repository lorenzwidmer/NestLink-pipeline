#!/usr/bin/env python3

import argparse
import re
import dnaio
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

def extract_sequences(file_path, flycode, orf1, orf2=None):
    """
    Parses a fasta file, and extracts the flycode, orf1 and orf2 (optional) sequences.

    Parameters:
        file_path (str): Fasta file to parse.
        flycode (tuple(str, str)): Start and end patterns of the flycode.
        orf1 (tuple(str, str)): Start and end patterns of orf1.
        orf2 (tuple(str, str)): Start and end patterns of orf2.

    Returns:
        dict. A dictionary where keys are cluster IDs and values are tuples containing the extracted flycode, ORF1, and optionally ORF2 sequences.
    """
    def is_in_frame(start_match, end_match):
        return (end_match.start() - start_match.end()) % 3 == 0

    def extract_subsequence(sequence, start_match, end_match):
        return sequence[start_match.end():end_match.start()]
    
    sequences = {}

    fc_start = re.compile(flycode[0], re.IGNORECASE)
    fc_end = re.compile(flycode[1], re.IGNORECASE)

    orf1_start = re.compile(orf1[0], re.IGNORECASE)
    orf1_end = re.compile(orf1[1], re.IGNORECASE)

    if orf2:
        orf2_start = re.compile(orf2[0], re.IGNORECASE)
        orf2_end = re.compile(orf2[1], re.IGNORECASE)

    with dnaio.open(file_path, mode="r") as reader:
        for record in reader:
            cluster_id = record.name
            sequence = record.sequence

            # Search for flycodes
            fc_start_match = fc_start.search(sequence)
            fc_end_match = fc_end.search(sequence)

            # Search for orf1
            orf1_start_match = orf1_start.search(sequence)
            orf1_end_match = orf1_end.search(sequence)

            # Search for orf2
            if orf2:
                orf2_start_match = orf2_start.search(sequence)
                orf2_end_match = orf2_end.search(sequence)

            # Extract sequences if matches are found and in-frame
            if fc_start_match and fc_end_match and is_in_frame(fc_start_match, fc_end_match):
                flycode_sequence = extract_subsequence(sequence, fc_start_match, fc_end_match)

                if orf1_start_match and orf1_end_match and is_in_frame(orf1_start_match, orf1_end_match):
                    orf1_sequence = extract_subsequence(sequence, orf1_start_match, orf1_end_match)

                    if orf2 and orf2_start_match and orf2_end_match and is_in_frame(orf2_start_match, orf2_end_match):
                        orf2_sequence = extract_subsequence(sequence, orf2_start_match, orf2_end_match)
                        sequences[cluster_id] = (flycode_sequence, orf1_sequence, orf2_sequence)
                    else:
                        sequences[cluster_id] = (flycode_sequence, orf1_sequence)
    return sequences

def translate(sequence):
    return str(Seq(sequence).translate())

def compare_aa_sequences(sequence, reference):
    changes = []

    if len(sequence) == len(reference):

        for position in range(len(reference)):
            if sequence[position] != reference[position]:
                change = f"{reference[position]}{position + 1}{sequence[position]}"
                changes.append(change)

        if reference == sequence:
            changes.append('wt')

    return changes

def get_aa_changes(sequences, reference, orf2=False):
    """
    Determines the amino acid changes between a reference sequence and consensus sequences.clear

    Parameters:
        sequences (dict):
        refereces (dict):
        orf2 (bool):

    Returns:
        dict
    """
    aa_changes = {}

    # Getting reference sequences on a protein level
    if len(reference) != 1:
        raise ValueError(
            "The reference file does not contain exactly one record."
    )
    _, reference_sequence = next(iter(reference.items()))
    reference_orf1_aa = translate(reference_sequence[1])
    if orf2:
        reference_orf2_aa = translate(reference_sequence[2])

    # Looping through all sequences
    for cluster_id, values in sequences.items():
        flycode_aa = translate(values[0])
        orf1_aa = translate(values[1])
        orf1_aa_changes = compare_aa_sequences(orf1_aa, reference_orf1_aa)
        if orf2:
            orf2_aa = translate(values[2])
            orf2_aa_changes = compare_aa_sequences(orf2_aa, reference_orf2_aa)
            aa_changes[cluster_id]=(flycode_aa, orf1_aa_changes, orf2_aa_changes)
        else:
            aa_changes[cluster_id]=(flycode_aa, orf1_aa_changes)
    return aa_changes

def write_flycode_db(aa_changes, output_path, prefix, orf1, orf2=False):
    """
    Writes a flycode database in fasta format to disk.

    Parameters:
        aa_changes (dict): A dictionairy containing flycodes and amino acid changes.
        ouput_path (str): Filename/ path of the flycode database to be written.
        orf1 (srt): Name of orf1.
        orf2 (str): Name of orf2.
    """
    variants = defaultdict(list)

    for cluster_id, values in aa_changes.items():
        flycode = values[0]
        if orf2:
            orf_changes = tuple(values[1]), tuple(values[2])   # both orfs
        else:
            orf_changes = tuple(values[1]) # only orf1
        variants[orf_changes].append(flycode)

    with open(output_path, "w") as file:
        if orf2:
            for variant, flycodes in variants.items():
                variant_orf1 = ",".join(variant[0])
                variant_orf2 = ",".join(variant[1])
                file.write(f">{prefix}|{orf1}={variant_orf1};{orf2}={variant_orf2}\n")
                file.write(f"{''.join(flycodes)}\n")
        else:
            for variant, flycodes in variants.items():
                variant_orf1 = ",".join(variant)
                file.write(f">{prefix}|{orf1}={variant_orf1}\n")
                file.write(f"{''.join(flycodes)}\n")

def main(assembly_path, reference_path, output_path, flycode, orf1, orf2):
    sequences = extract_sequences(assembly_path, flycode, orf1, orf2)
    reference = extract_sequences(reference_path, flycode, orf1, orf2)

    aa_changes = get_aa_changes(sequences, reference, True)

    write_flycode_db(aa_changes, output_path, "NL22", "TM287", "TM288")

    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly_path", type=str)
    parser.add_argument("--reference_path", type=str)
    parser.add_argument("--output", type=str)
    parser.add_argument("--flycode", nargs="+", type=str)
    parser.add_argument("--orf1", nargs="+", type=str)
    parser.add_argument("--orf2", nargs="+", type=str)


    args = parser.parse_args()
    main(args.assembly_path, args.reference_path, args.output, args.flycode, args.orf1, args.orf2)