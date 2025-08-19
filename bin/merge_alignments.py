#!/usr/bin/env python3
import argparse
from pathlib import Path
import subprocess
import tempfile


def merge_bam_files(input_folder, output_bam, threads):
    """
    Merge all BAM files from the folder 'input_folder' in batches, then produce a final merged file 'output_bam'.

    Args:
        input_folder (str): Folder containing BAM files to merge.
        output_bam (str): Path to the final merged BAM file.
        threads (int): Number of threads to use. 
    """
    max_files = 200

    input_folder = Path(input_folder)
    if not input_folder.exists() or not input_folder.is_dir():
        raise ValueError(f"Input directory '{input_folder}' does not exist or is not a directory.")

    output_bam = Path(output_bam)

    # Get all BAM files
    bam_files = sorted(input_folder.glob("*.bam"))
    if not bam_files:
        raise ValueError(f"No BAM files found in '{input_folder}.'")

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = Path(temp_dir)
        temp_files = []

        print(f"Merging {len(bam_files)} BAMs with samtools merge.")

        # Merge in batches of max_files
        for i in range(0, len(bam_files), max_files):
            batch_files = bam_files[i:i + max_files]
            temp_file = temp_dir / f"temp_merged_{i // max_files}.bam"
            print(f"({i // max_files + 1}/{(len(bam_files) - 1) // max_files + 1}) Merging {len(batch_files)} BAMs into {temp_file.name}.")
            merge_bams(batch_files, temp_file, threads)
            temp_files.append(temp_file)

        # Merge all temporary files
        print(f"Merging {len(temp_files)} temporary BAMs into {output_bam}.")
        merge_bams(temp_files, output_bam, threads)


def merge_bams(files_list, output_file, threads):
    """
    Use samtools to merge a list of BAM files into a single output BAM file.

    Args:
        files_list (list[Path]): Paths to the BAM files to merge.
        output_file (Path): Path where the merged BAM file will be created.
        threads (int): Number of threads to use. 
    """
    cmd = ["samtools", "merge", "-@", str(threads), output_file] + files_list
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=str)
    parser.add_argument("-o", "--output", required=True, type=str)
    parser.add_argument("-t", "--threads", required=False, type=int, default=1)
    args = parser.parse_args()
    merge_bam_files(args.input, args.output, args.threads)
