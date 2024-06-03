#!/usr/bin/env bash
set -e
set -u

reads=${1}

nb06="GACTACTTTCTGCCTTTGCGAGAA...TTCTCGCAAAGGCAGAAAGTAGTC"
output_file="${reads}_cut.fastq.gz"
cutadapt -g "${nb06}" -m "1700" -M "1900" --discard-untrimmed  --cores=0 -o "${output_file}" "${reads}.fastq.gz"