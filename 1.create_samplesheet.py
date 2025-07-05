#!/usr/bin/env python

import argparse
import os

def get_fastq_files(directory):
    fastq_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".fastq.gz") or file.endswith(".fq.gz"):
                fastq_files.append(os.path.join(root, file))
    return fastq_files

def generate_sample_sheet(directory, sample_sheet_file):
    abs_directory = os.path.abspath(directory)
    fastq_files = get_fastq_files(abs_directory)
    with open(sample_sheet_file, "w") as f:
        f.write("sample_name\tread1\tread2\n")
        for fastq in fastq_files:
            sample_name = fastq.split("/")[-1].split("_")[0]
            if "_R1_" in fastq:
                read1 = fastq
                read2 = fastq.replace("_R1_", "_R2_")
                if read2 in fastq_files:
                    f.write(f"{sample_name}\t{read1}\t{read2}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate sample sheet from fastq files')
    parser.add_argument('directory', type=str, help='Directory containing fastq files')
    parser.add_argument('sample_sheet_file', type=str, help='Output sample sheet file')
    args = parser.parse_args()

    directory = args.directory
    sample_sheet_file = args.sample_sheet_file

    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory.")
        exit(1)

    if os.path.exists(sample_sheet_file):
        overwrite = input(f"The file {sample_sheet_file} already exists. Do you want to overwrite it? (y/n) ")
        if overwrite.lower() != "y":
            print("Sample sheet generation canceled.")
            exit(0)

    generate_sample_sheet(directory, sample_sheet_file)
