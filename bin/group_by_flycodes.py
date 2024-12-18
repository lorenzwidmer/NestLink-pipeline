#!/usr/bin/env python3
import argparse
import subprocess
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


def main(flycodes):
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

    # Writing dataframe to disk as csv.
    flycodes_df.write_csv("flycodes.csv")
    clusters_df.write_csv("clusters.csv")
    mapped_flycodes_df.write_csv("mapped_flycodes.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--flycodes", type=str)
    args = parser.parse_args()
    main(args.flycodes)
