#!/usr/bin/env python3
import os
import subprocess
import glob


def merge_bam_files(max_files=256):
    temp_dir = "temp"
    os.makedirs(temp_dir, exist_ok=True)

    # Get all BAM files
    bam_files = glob.glob("alignments/*.bam")
    temp_files = []

    # Merge in batches of max_files
    for i in range(0, len(bam_files), max_files):
        batch_files = bam_files[i:i + max_files]
        temp_file_path = os.path.join(temp_dir, f"temp_merged_{i // max_files}.bam")
        merge_bams(batch_files, temp_file_path)
        temp_files.append(temp_file_path)

    # Merge all temporary files
    merge_bams(temp_files, "merged.bam")


def merge_bams(files_list, output_file):
    cmd = ["samtools", "merge", output_file] + files_list
    subprocess.run(cmd, check=True)


def main():
    merge_bam_files()

    cmd_sort = ["samtools", "sort", "merged.bam", "-o", "merged.sorted.bam"]
    cmd_index = ["samtools", "index", "merged.sorted.bam"]
    subprocess.run(cmd_sort)
    subprocess.run(cmd_index)


if __name__ == "__main__":
    main()
