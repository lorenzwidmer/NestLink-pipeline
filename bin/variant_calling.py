#!/usr/bin/env python3
import argparse
import re
import dnaio
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq


def extract_sequences(file_path, flycode, orf):
    """
    Parse a FASTA file to extract flycode and ORF sequences.

    Args:
        file_path (str): Path to the FASTA file.
        flycode (tuple): A (start_pattern, end_pattern) tuple for the flycode.
        orf (tuple): A (start_pattern, end_pattern) tuple for the ORF.

    Returns:
        dict: Keys are cluster IDs; values are tuples of extracted sequences 
              (flycode_seq, orf_seq) in-frame only.
    """
    def is_in_frame(start_match, end_match):
        return (end_match.start() - start_match.end()) % 3 == 0

    def extract_subsequence(sequence, start_match, end_match):
        return sequence[start_match.end():end_match.start()]

    sequences = {}

    fc_start = re.compile(flycode[0], re.IGNORECASE)
    fc_end = re.compile(flycode[1], re.IGNORECASE)

    orf_start = re.compile(orf[0], re.IGNORECASE)
    orf_end = re.compile(orf[1], re.IGNORECASE)

    with dnaio.open(file_path, mode="r") as reader:
        for record in reader:
            cluster_id = record.name
            sequence = record.sequence

            # Search for flycodes
            fc_start_match = fc_start.search(sequence)
            fc_end_match = fc_end.search(sequence)

            # Search for orf
            orf_start_match = orf_start.search(sequence)
            orf_end_match = orf_end.search(sequence)

            # Extract sequences if matches are found and in-frame
            if fc_start_match and fc_end_match and is_in_frame(fc_start_match, fc_end_match):
                flycode_sequence = extract_subsequence(sequence, fc_start_match, fc_end_match)

                if orf_start_match and orf_end_match and is_in_frame(orf_start_match, orf_end_match):
                    orf_sequence = extract_subsequence(sequence, orf_start_match, orf_end_match)
                    sequences[cluster_id] = (flycode_sequence, orf_sequence)
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


def get_aa_changes(sequences, reference):
    """
    Determine amino acid changes for each consensus sequence relative to a the reference.

    Args:
        sequences (dict): Dictionary mapping cluster IDs to extracted (flycode, ORF) sequences.
        reference (dict): Dictionary containing exactly one reference entry of the same format.

    Returns:
        dict: Maps each cluster ID to a tuple of (flycode_aa, ORF_changes).
    """
    aa_changes = {}

    # Getting reference sequences on a protein level
    if len(reference) != 1:
        raise ValueError(
            "The reference file does not contain exactly one record."
        )
    _, reference_sequence = next(iter(reference.items()))
    reference_orf_aa = translate(reference_sequence[1])

    # Looping through all sequences
    for cluster_id, values in sequences.items():
        flycode_aa = translate(values[0])
        orf_aa = translate(values[1])
        orf_aa_changes = compare_aa_sequences(orf_aa, reference_orf_aa)
        aa_changes[cluster_id] = (flycode_aa, orf_aa_changes)
    return aa_changes


def write_flycode_db(aa_changes, output_path, prefix, orf):
    """
    Write a flycode database in FASTA format, grouping identical variants.

    Args:
        aa_changes (dict): A dictionary of cluster IDs to (flycode_aa, ORF_changes).
        output_path (str): Path to the output FASTA file.
        prefix (str): Prefix (e.g., experiment name) to include in FASTA headers.
        orf (str): Name of ORF.
    """
    variants = defaultdict(list)

    for cluster_id, values in aa_changes.items():
        flycode = values[0]
        orf_changes = tuple(values[1])
        variants[orf_changes].append(flycode)

    with open(output_path, "w") as file:
        for variant, flycodes in variants.items():
            variant_orf = ",".join(variant)
            file.write(f">{prefix}|{orf}={variant_orf}\n")
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
    orf_name, orf_pattern = args.orf_name, args.orf_pattern
    experiment_name = args.experiment_name

    sequences = extract_sequences(assembly_path, flycode_pattern, orf_pattern)
    reference = extract_sequences(reference_path, flycode_pattern, orf_pattern)

    aa_changes = get_aa_changes(sequences, reference)

    write_flycode_db(aa_changes, output, experiment_name, orf_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly_path", type=str, required=True)
    parser.add_argument("--reference_path", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)
    parser.add_argument("--experiment_name", type=str, required=True)
    parser.add_argument("--flycode_pattern", nargs=2, type=str, required=True)
    parser.add_argument("--orf_name", type=str, required=True)
    parser.add_argument("--orf_pattern", nargs=2, type=str, required=True)

    args = parser.parse_args()
    main(args)
