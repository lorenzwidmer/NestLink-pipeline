#!/usr/bin/env bash
set -e
set -u

# Set the number of threads
threads=$1

mkdir alignments

# Loop through all the clustered reads fastq.gz files
for reads_file in clusters/*.fastq.gz; do
    uuid=$(basename "$reads_file" .fastq.gz)
 
    reference_file="references/${uuid}.fasta"
    output_bam="alignments/cluster_${uuid}"

    echo "Alligning reads for flycode: $uuid."
    mini_align.sh -i $reads_file -r $reference_file -m -p $output_bam -t $threads

done