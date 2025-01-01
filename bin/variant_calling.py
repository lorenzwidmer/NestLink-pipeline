#!/usr/bin/env python3
import argparse
import re
import dnaio
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq


def extract_sequences(file_path, flycode, orf1, orf2=None):
    """
    Parse a FASTA file to extract flycode, ORF1, and optional ORF2 sequences.

    Args:
        file_path (str): Path to the FASTA file.
        flycode (tuple): A (start_pattern, end_pattern) tuple for the flycode.
        orf1 (tuple): A (start_pattern, end_pattern) tuple for ORF1.
        orf2 (tuple, optional): A (start_pattern, end_pattern) tuple for ORF2 (default: None).

    Returns:
        dict: Keys are cluster IDs; values are tuples of extracted sequences 
              (flycode_seq, orf1_seq, [orf2_seq if provided]) in-frame only.
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
    """
    Translate a nucleotide sequence into an amino acid sequence.

    Args:
        sequence (str): The nucleotide sequence.

    Returns:
        str: The corresponding amino acid sequence.
    """
    return str(Seq(sequence).translate())


def compare_aa_sequences(sequence, reference):
    """
     Compare two amino acid sequences and identify positions where they differ.

    Args:
        sequence (str): The amino acid sequence to compare.
        reference (str): The reference amino acid sequence.

    Returns:
        list: A list of amino acid changes, in the format '<refAA><position><newAA>' or 'wt' if they are identical.
    """
    changes = []

    if len(sequence) == len(reference):

        for position in range(len(reference)):
            if sequence[position] != reference[position]:
                change = f"{reference[position]}{position + 1}{sequence[position]}"
                changes.append(change)

        if reference == sequence:
            changes.append("wt")

    return changes


def get_aa_changes(sequences, reference, orf2=None):
    """
    Determine amino acid changes for each cluster ID relative to a the reference.

    Args:
        sequences (dict): Dictionary mapping cluster IDs to extracted (flycode, ORF1, [ORF2]) sequences.
        reference (dict): Dictionary containing exactly one reference entry of the same format.
        orf2 (str): Name of the second ORF if present, otherwise None.

    Returns:
        dict: Maps each cluster ID to a tuple of (flycode_aa, ORF1_changes, [ORF2_changes]).
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
            aa_changes[cluster_id] = (flycode_aa, orf1_aa_changes, orf2_aa_changes)
        else:
            aa_changes[cluster_id] = (flycode_aa, orf1_aa_changes)
    return aa_changes


def write_flycode_db(aa_changes, output_path, prefix, orf1, orf2=None):
    """
    Write a flycode database in FASTA format, grouping identical variants.

    Args:
        aa_changes (dict): A dictionary of cluster IDs to (flycode_aa, ORF1_changes, [ORF2_changes]).
        output_path (str): Path to the output FASTA file.
        prefix (str): Prefix (e.g., experiment name) to include in FASTA headers.
        orf1 (str): Name of ORF1.
        orf2 (str, optional): Name of ORF2 if present, otherwise None.
    """
    variants = defaultdict(list)

    for cluster_id, values in aa_changes.items():
        flycode = values[0]
        if orf2:
            orf_changes = tuple(values[1]), tuple(values[2])   # both orfs
        else:
            orf_changes = tuple(values[1])  # only orf1
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


def main(args):
    """
    Coordinate the extraction, translation, variant calling, and database writing of flycodes.

    Args:
        args (Namespace): Parsed command-line arguments.
    """
    assembly_path = args.assembly_path
    reference_path = args.reference_path
    output = args.output
    flycode_pattern = args.flycode_pattern
    orf1_name, orf1_pattern = args.orf1_name, args.orf1_pattern
    orf2_name, orf2_pattern = args.orf2_name, args.orf2_pattern
    experiment_name = args.experiment_name

    sequences = extract_sequences(assembly_path, flycode_pattern, orf1_pattern, orf2_pattern)
    reference = extract_sequences(reference_path, flycode_pattern, orf1_pattern, orf2_pattern)

    aa_changes = get_aa_changes(sequences, reference, orf2_name)

    write_flycode_db(aa_changes, output, experiment_name, orf1_name, orf2_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly_path", type=str, required=True)
    parser.add_argument("--reference_path", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)
    parser.add_argument("--experiment_name", type=str, required=True)
    parser.add_argument("--flycode_pattern", nargs=2, type=str, required=True)
    parser.add_argument("--orf1_name", type=str, required=True)
    parser.add_argument("--orf1_pattern", nargs=2, type=str, required=True)
    parser.add_argument("--orf2_name", type=str, default=None)
    parser.add_argument("--orf2_pattern", nargs=2, type=str, default=None)

    args = parser.parse_args()
    main(args)
