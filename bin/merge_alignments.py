#!/usr/bin/env python3
import argparse
import os
import subprocess
import glob


def merge_bam_files(output_path, max_files=256):
    # Create a directory for temporary merged files if it does not exist
    temp_dir = os.path.join(output_path, "temp_merged")
    os.makedirs(temp_dir, exist_ok=True)

    # Get all BAM files
    bam_files = glob.glob(os.path.join(output_path, "*.bam"))
    temp_files = []

    # Merge in batches of max_files
    for i in range(0, len(bam_files), max_files):
        batch_files = bam_files[i:i + max_files]
        temp_file_path = os.path.join(temp_dir, f"temp_merged_{i // max_files}.bam")
        merge_files(batch_files, temp_file_path)
        temp_files.append(temp_file_path)

    # Merge all temporary files
    merge_files(temp_files, "merged.bam")


def merge_files(files_list, output_file):
    # Prepare the samtools merge command
    cmd = ["samtools", "merge", output_file] + files_list
    # Run the command
    subprocess.run(cmd, check=True)


def sort_index_bam(bam_path):
    cmd_sort = ["samtools", "sort", "merged.bam", "-o", "merged.sorted.bam"]
    cmd_index = ["samtools", "index", "merged.sorted.bam"]
    subprocess.run(cmd_sort)
    subprocess.run(cmd_index)


def main(bam_files):
    merge_bam_files(bam_files)
    sort_index_bam(bam_files)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam_files", type=str)

    args = parser.parse_args()
    main(args.bam_files)
