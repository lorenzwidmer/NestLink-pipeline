#!/usr/bin/env bash
set -e
set -u
set -o pipefail

threads=1

while getopts "i:r:o:t:" opt; do
  case $opt in
    i) input_folder=$OPTARG;;
    r) reference_folder=$OPTARG;;
    o) output_folder=$OPTARG;;
    t) threads=$OPTARG;;
  esac
done

mkdir "$output_folder"

# Loop through all the clustered reads fastq.gz files
for reads_file in ${input_folder}/*.fastq.gz; do
    uuid=$(basename "$reads_file" .fastq.gz)
 
    reference_file="${reference_folder}/${uuid}.fasta"
    output_bam="${output_folder}/${uuid}.bam"

    # Indexing reference
    samtools faidx "$reference_file"
    minimap2 -x map-ont -d "${reference_file}.mmi" "$reference_file"

    # Aligning reads to references; filtering unaligned, secondary and supplementary reads; sorting bam
    minimap2 -x map-ont --secondary=no -L --MD -t "$threads" -a "${reference_file}.mmi" "$reads_file" |
        samtools view -@ 8 -T "$reference_file" -F 2308 -b - |
        samtools sort -@ 8 -o "$output_bam" -
done