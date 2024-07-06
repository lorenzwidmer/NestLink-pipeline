#!/usr/bin/env python3
import argparse
import os
import dnaio
from collections import defaultdict


def parse_centroids(file_path, min_cluster_sizes):
    """
    Reads the flycode_centroids.fasta file and returns a list of clusterids with clustersizes equal or larger than min_cluster_sizes.
    """
    centroids = []
    with dnaio.open(file_path) as reader:
        for record in reader:
            _, clusterid, size = parse_header(record.name)
            if size >= min_cluster_sizes:
                centroids.append(clusterid)
    return centroids


def parse_header(fasta_header):
    """
    Parses the header of a flycode_centroids file and returns a tuple containg the readid, direction of where the flycode was found, the clusterid and the size of the cluster.
    """
    parts = fasta_header.split(';')
    readid = parts[0]
    clusterid = int(parts[1].split('=')[1])
    size = int(parts[2].split('=')[1])
    return (readid, clusterid, size)


def create_read_dictionary(folder_path, cluster_ids):
    """
    Returns a dictionary with read ids as key and their corresponding cluster_id as a value for a list of cluster_ids.
    """
    reads = {}
    for cluster_id in cluster_ids:
        file_name = 'cluster.fasta' + str(cluster_id)
        file_path = os.path.join(folder_path, file_name)
        with dnaio.open(file_path) as reader:
            for record in reader:
                readid, _ = record.name.split(';')
                reads[readid] = cluster_id
    return reads


def bin_reads_by_flycodes(reads_path, reads_dict):
    """
    Returns a dictionary with clusterids as key and a list of corresponding reads as value.
    """
    reads_binned = defaultdict(list)
    with dnaio.open(reads_path) as reader:
        for record in reader:
            readid = record.name.split('	')[0]
            if readid in reads_dict:
                cluster_id = reads_dict[readid]
                reads_binned[cluster_id].append(record)
    return reads_binned


def write_binned_reads(folder_path, binned_reads):
    '''
    Writes new fastq files containing reads binned by flycodes. The direction of the reads is all fwd.
    '''
    for cluster_id, reads in binned_reads.items():
        file_name = 'reads_fc' + str(cluster_id) + '.fastq.gz'
        file_path = os.path.join(folder_path, file_name)
        with dnaio.open(file_path, mode="w") as writer:
            for record in reads:
                writer.write(record)


def main(centroids, clusters, sequences):
    cluster_ids = parse_centroids(centroids, 30)
    reads_dict = create_read_dictionary(clusters, cluster_ids)
    binned_reads = bin_reads_by_flycodes(sequences, reads_dict)
    write_binned_reads("binned_sequences", binned_reads)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--centroids", type=str)
    parser.add_argument("--clusters", type=str)
    parser.add_argument("--sequences", type=str)

    args = parser.parse_args()
    main(args.centroids, args.clusters, args.sequences)
