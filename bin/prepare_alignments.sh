#!/usr/bin/env bash
set -e
set -u

# Set the number of threads to use
threads=$1

mkdir alignments

# Loop through all the reads fastq files in the input path
for reads_file in clusters/*.fastq.gz; do
    uuid=$(basename "$reads_file" .fastq.gz)
 
    # Generate the corresponding seed file and output bam file names
    reference_file="references/${uuid}.fasta"
    output_bam="alignments/calls_to_draft_fc${uuid}"

    echo "Processing reads file: $reads_file."
    mini_align.sh -i $reads_file -r $reference_file -m -p $output_bam -t $threads

done