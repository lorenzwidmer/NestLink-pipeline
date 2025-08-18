#!/usr/bin/env python3
import os
import subprocess
import glob


def merge_bam_files(max_files=256):
    """
    Merge all BAM files from the 'alignments' folder in batches, then produce a final merged file.

    Args:
        max_files (int, optional): Maximum number of files to merge in one batch. Default is 256.
    """
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
    """
    Use samtools to merge a list of BAM files into a single output BAM file.

    Args:
        files_list (list[str]): Paths to the BAM files to merge.
        output_file (str): Path where the merged BAM file will be created.
    """
    cmd = ["samtools", "merge", output_file] + files_list
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    merge_bam_files()
