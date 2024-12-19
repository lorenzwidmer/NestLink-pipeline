#!/usr/bin/env python3
import argparse
import subprocess
from collections import defaultdict
import dnaio
import duckdb
import polars as pl


def flycodes_to_dataframe(file_path):
    """
    Reads a fasta file containing all flycodes from a sequencing experiment extracted by cutadapt.
    Returns a polars dataframe with the read_ids and their corresponding flycodes.
    """
    data = []

    with dnaio.open(file_path) as reader:
        for record in reader:
            row = {"read_id": record.name, "flycode": record.sequence}
            data.append(row)

    return pl.DataFrame(data)


def map_flycodes_to_clusters(file_path):
    """
    Reads a sam file of flycodes aligned to a reference containing high-quality flycodes "clusters".
    Returns a polars dataframe with the read_ids and their corresponding cluster_ids.
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
            read_id = fields[0]     # Barcode name
            cluster_id = fields[2]  # Reference cluster name (aligned to)
            flag = int(fields[1])   # SAM flag

            # Check if the alignment is valid (flag 4 means unaligned)
            if cluster_id != "*" and not (flag & 4):
                row = {"read_id": read_id, "cluster_id": cluster_id}
                data.append(row)

    return pl.DataFrame(data)


def bin_reads_by_flycodes(file_path, flycode_map):
    """
    Returns a dictionary with cluster_ids as key and a list of corresponding reads (records) as value.
    """
    binned_reads = defaultdict(list)
    with dnaio.open(file_path) as reader:
        for record in reader:
            readid = record.name.split('	')[0]
            if readid in flycode_map:
                cluster_id = flycode_map[readid]
                binned_reads[cluster_id].append(record)
    return binned_reads


def write_binned_reads(binned_reads):
    """
    Writes new fastq.gz files containing reads binned by flycodes for alignment. The direction of the reads is all fwd.
    """
    subprocess.run(["mkdir", "clusters"])

    for cluster_id, reads in binned_reads.items():
        file_path = f"clusters/{cluster_id}.fastq.gz"
        with dnaio.open(file_path, mode="w") as writer:
            for record in reads:
                writer.write(record)


def write_references(clusters_df, reference_seq):
    """
    Writes a reference sequence fasta file for every cluster.
    Writes a fasta file containing the reference sequences of all clusters.
    """
    subprocess.run(["mkdir", "references"])

    # Open and read the reference sequence.
    with dnaio.open(reference_seq, mode="r") as reader:
        for record in reader:
            reference = record.sequence
            break
    
    # Splitting the reference intwo two parts using the flycode as delimiter.
    reference = reference.split("GGTAGTNNNNNNNNNNNNNNNNNNNNNTGGcgg")

    records = [] # for storing all cluster records.

    # Writing individual reference files for each cluster.
    for cluster_id, flycode in clusters_df.rows():
        file_name = f"references/{cluster_id}.fasta"
        reference_sequence = f"{reference[0]}{flycode}{reference[1]}"
        record = dnaio.SequenceRecord(cluster_id, reference_sequence)
        records.append(record)
        with dnaio.open(file_name, mode="w") as writer:
            writer.write(record)

    # Writing a reference file containing all clusters (used by medaka sequence).
    file_name = "references.fasta"
    with dnaio.open(file_name, mode="w") as writer:
        for record in records:
            writer.write(record)


def main(flycodes, sequences, reference_seq):
    """
    Main function for binning reads by their flycodes and preparing them for the alignment process.
    """
    # Reading in the flycodes
    flycodes_df = flycodes_to_dataframe(flycodes)

    # Add a new column with True/False indicating a valid flycode
    valid_flycode = r"^GGTAGT(GCA|GTT|GAT|CCA|GAA|ACT|GGT|TCT|TAC|CTG|TGG|CAG|TTC|AAC){6,8}(TGGCGG|TGGCTGCGG|TGGCAGTCTCGG|TGGCAGGAAGGAGGTCGG)$"
    flycodes_df = flycodes_df.with_columns(
        pl.col("flycode").str.contains(valid_flycode).alias("is_valid_flycode")
    )

    # Getting high-quality flycodes, they have to be valid and be seen more than 10 times.
    clusters_df = duckdb.sql(
        """
        SELECT 
            gen_random_uuid() AS cluster_id, 
            flycode
        FROM 
            flycodes_df
        WHERE 
            is_valid_flycode
        GROUP BY 
            flycode
        HAVING 
            COUNT(flycode) >= 10;
        """
    ).pl()

    # Writing the high-quality flycodes into clusters.fasta to be used as a reference for alignment.
    with dnaio.open("clusters.fasta", mode="w") as writer:
        for cluster_id, sequence in clusters_df.rows():
            writer.write(dnaio.SequenceRecord(cluster_id, sequence))

    # Indexing the reference.
    subprocess.run(["bwa", "index", "clusters.fasta"])

    # Alining all flycodes to the reference.
    with open("flycodes_to_clusters.sai", "w") as sai_file:
        subprocess.run(["bwa", "aln", "clusters.fasta", flycodes], stdout=sai_file)

    # Generating alignments in the SAM format
    with open("flycodes_to_clusters.sam", "w") as sam_file:
        subprocess.run(["bwa", "samse", "clusters.fasta", "flycodes_to_clusters.sai", flycodes], stdout=sam_file)

    # Mapping flycodes to their clusters.
    mapped_flycodes_df = map_flycodes_to_clusters("flycodes_to_clusters.sam")

    # Filtering mapped_flycodes_df to only contain up to 100 reads per cluster.
    mapped_flycodes_df = duckdb.sql(
        """
        WITH ranked_reads AS (
            SELECT 
                read_id,
                cluster_id,
                row_number() OVER (PARTITION BY cluster_id ORDER BY read_id) AS row_num
            FROM mapped_flycodes_df
        )
        SELECT 
            read_id,
            cluster_id
        FROM ranked_reads
        WHERE row_num <= 100;
        """
    ).pl()

    # Converting the flycode map into a dict.
    flycode_map = {row["read_id"]: row["cluster_id"] for row in mapped_flycodes_df.iter_rows(named=True)}

    # Binning reads by flycode cluster, storing them in a dict.
    binned_reads = bin_reads_by_flycodes(sequences, flycode_map)

    # Writing binned reads and references to disk.
    write_binned_reads(binned_reads)
    write_references(clusters_df, reference_seq)

    # Writing dataframes to disk as csv.
    flycodes_df.write_csv("flycodes.csv")
    clusters_df.write_csv("clusters.csv")
    mapped_flycodes_df.write_csv("mapped_flycodes.csv")


if __name__ == "__main__":
    # Reading arguments and calling the main function.
    parser = argparse.ArgumentParser()
    parser.add_argument("--flycodes", type=str)
    parser.add_argument("--sequences", type=str)
    parser.add_argument("--reference_seq", type=str)
    args = parser.parse_args()
    main(args.flycodes, args.sequences, args.reference_seq)
