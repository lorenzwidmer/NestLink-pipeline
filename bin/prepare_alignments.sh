#!/usr/bin/env bash
set -e
set -u
set -o pipefail

# Set the number of threads
threads=$1

mkdir alignments

# Loop through all the clustered reads fastq.gz files
for reads_file in clusters/*.fastq.gz; do
    uuid=$(basename "$reads_file" .fastq.gz)
 
    reference_file="references/${uuid}.fasta"
    output_bam="alignments/${uuid}.bam"

    # Indexing reference
    samtools faidx $reference_file
    minimap2 -x map-ont -d ${reference_file}.mmi $reference_file

    # Aligning reads to references; filtering unaligned, secondary and supplementary reads; sorting bam
    minimap2 -x map-ont --secondary=no -L --MD -t $threads -a ${reference_file}.mmi $reads_file |
        samtools view -@ 8 -T $reference_file -F 2308 -b - |
        samtools sort -@ 8 -o $output_bam -
done