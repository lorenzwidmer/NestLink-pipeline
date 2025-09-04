#!/usr/bin/env python3
import argparse
import subprocess
from collections import defaultdict
import dnaio
import duckdb
import polars as pl


def barcodes_to_dataframe(file_path):
    """
    Read a FASTA of extracted barcodes into a polars DataFrame with read IDs and their corresponding barcodes.

    Args:
        file_path (str): Path to the FASTA file with extracted barcodes.

    Returns:
        pl.DataFrame: A DataFrame with 'read_id' and 'barcode' columns.
    """
    data = []

    with dnaio.open(file_path) as reader:
        for record in reader:
            row = {"read_id": record.name, "barcode": record.sequence}
            data.append(row)

    if not data:
        raise ValueError(
            f"No barcodes found in file '{file_path}'. "
            "Please check the barcode and or sequence extraction processes."
        )

    return pl.DataFrame(data)


def map_barcodes_to_clusters(file_path):
    """
    Parse a SAM file of barcodes aligned to high-quality barcodes (clusters) into a polars DataFrame mapping read IDs to cluster IDs.

    Args:
        file_path (str): Path to the SAM file with aligned barcodes.

    Returns:
        pl.DataFrame: A DataFrame with 'read_id', 'cluster_id', 'cigar' and 'edit_distance' columns.
    """
    data = []

    with open(file_path, "r") as f:
        for line in f:
            # Skip header lines
            if line.startswith("@"):
                continue

            # Split the SAM fields
            fields = line.strip().split("\t")

            # Extract necessary fields
            read_id = fields[0]     # Barcode name (Query template NAME - QNAME)
            cluster_id = fields[2]  # Reference cluster name (Reference sequence NAME - RNAME)
            flag = int(fields[1])   # SAM flag
            cigar = fields[5]       # CIGAR string

            # Check if the alignment is valid (flag 4 means unaligned)
            if cluster_id != "*" and not (flag & 4):
                edit_distance = fields[12].split(":")[2]
                row = {"read_id": read_id, "cluster_id": cluster_id, "cigar": cigar, "edit_distance": edit_distance}
                data.append(row)

    if not data:
        raise ValueError(
            f"No aligned barcodes found in file '{file_path}'."
        )

    return pl.DataFrame(data, schema_overrides={"edit_distance": pl.UInt8})


def write_references(clusters_df, reference_seq):
    """
    Write reference FASTA file containing all clusters.

    Args:
        clusters_df (pl.DataFrame): DataFrame with cluster IDs and their barcodes.
        reference_seq (str): Path to the reference sequence FASTA file.
    """
    # Open and read the reference sequence.
    with dnaio.open(reference_seq, mode="r") as reader:
        reference_record = None
        for record in reader:
            reference_record = record
            break  # stop after the first record

        if reference_record is None:
            raise ValueError(
                f"No record found in the reference sequence file '{reference_seq}'. "
                "Please check the file."
            )

        reference_sequence = reference_record.sequence

    records = []  # for storing all cluster records.
    for cluster_id, _ in clusters_df.rows():
        record = dnaio.SequenceRecord(cluster_id, reference_sequence)
        records.append(record)

    # Writing a reference file containing all clusters (used by dorado polish).
    file_name = "references.fasta"
    with dnaio.open(file_name, mode="w") as writer:
        for record in records:
            writer.write(record)


def main(sample_id, barcodes, reference_seq, barcode_regex):
    """
    Main function to validate and cluster barcodes and map reads to clusters.

    Args:
        sample_id (str): The sample id.
        barcodes (str): Path to the FASTA file with extracted barcodes.
        reference_seq (str): Path to the reference sequence FASTA file.
        barcode_regex (srt): Regex that matches the used barcode.
    """
    # Reading in the barcodes.
    barcodes_df = barcodes_to_dataframe(barcodes)

    # Add a new column with True/False indicating a valid barcode.
    barcodes_df = barcodes_df.with_columns(
        pl.col("barcode").str.contains(barcode_regex).alias("is_valid_barcode")
    )
    barcodes_df.write_csv(f"{sample_id}_reads.csv")

    # Getting high-quality barcodes, they have to be valid and be seen more than 10 times.
    clusters_df = duckdb.sql(
        """
        SELECT 
            gen_random_uuid() AS cluster_id, 
            barcode
        FROM barcodes_df
        WHERE is_valid_barcode
        GROUP BY barcode
        HAVING COUNT(barcode) >= 10;
        """
    ).pl()
    clusters_df.write_csv(f"{sample_id}_clusters.csv")

    # Writing the high-quality barcodes into clusters.fasta to be used as a reference for alignment.
    with dnaio.open("clusters.fasta", mode="w") as writer:
        for cluster_id, sequence in clusters_df.rows():
            writer.write(dnaio.SequenceRecord(cluster_id, sequence))

    # Indexing the reference.
    subprocess.run(["bwa", "index", "clusters.fasta"], check=True)

    # Aligning all barcodes to the reference.
    with open("barcodes_to_clusters.sai", "w") as sai_file:
        subprocess.run(["bwa", "aln", "-N", "-n 2", "clusters.fasta", barcodes], stdout=sai_file, check=True)

    # Generating alignments in the SAM format,
    with open("barcodes_to_clusters.sam", "w") as sam_file:
        subprocess.run(["bwa", "samse", "clusters.fasta", "barcodes_to_clusters.sai", barcodes], stdout=sam_file, check=True)

    # Mapping barcodes to their clusters.
    mapped_barcodes_df = map_barcodes_to_clusters("barcodes_to_clusters.sam")
    mapped_barcodes_df.write_csv(f"{sample_id}_mapped_reads.csv")

    # Filtering mapped_barcodes_df to only contain up to 100 reads with low edit distance barcodes per cluster.
    mapped_barcodes_df = duckdb.sql(
        """
        WITH ranked_reads AS (
            SELECT 
                read_id,
                cluster_id,
                cigar,
                edit_distance,
                row_number() OVER (PARTITION BY cluster_id ORDER BY edit_distance ASC) AS row_num
            FROM mapped_barcodes_df
        )
        SELECT 
            read_id,
            cluster_id,
            cigar,
            edit_distance
        FROM ranked_reads
        WHERE row_num <= 100;
        """
    ).pl()
    mapped_barcodes_df.write_csv(f"{sample_id}_mapped_reads_filtered.csv")

    # Writing references to disk.
    write_references(clusters_df, reference_seq)


if __name__ == "__main__":
    # Reading arguments and calling the main function.
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_id", type=str)
    parser.add_argument("--barcodes", type=str)
    parser.add_argument("--reference_seq", type=str)
    parser.add_argument("--barcode_regex", type=str)
    args = parser.parse_args()
    main(args.sample_id, args.barcodes, args.reference_seq, args.barcode_regex)
