#!/usr/bin/env python

import argparse
import os
import subprocess
import sys
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_command(command, step_name):
    """Execute a shell command and handle errors."""
    logging.info(f"Running {step_name}: {command}")
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        logging.error(f"{step_name} failed with return code {process.returncode}: {stderr.decode('utf-8')}")
        sys.exit(1)
    logging.info(f"{step_name} output: {stdout.decode('utf-8')}")
    if stderr:
        logging.info(f"{step_name} stderr: {stderr.decode('utf-8')}")
    return stdout.decode('utf-8'), stderr.decode('utf-8')

def check_file_exists(file_path, description):
    """Check if a file exists, exit if not."""
    if not os.path.exists(file_path):
        logging.error(f"{description} not found: {file_path}")
        sys.exit(1)

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Automate bioinformatics pipeline for variant calling and annotation.')
    parser.add_argument('--input_dir', type=str, required=True, help='Directory containing FASTQ files.')
    parser.add_argument('--reference_fasta', type=str, required=True, help='Path to reference FASTA file.')
    parser.add_argument('--genbank_file', type=str, required=True, help='Path to GenBank file for SnpEff database.')
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory for all results.')
    parser.add_argument('--database_name', type=str, default='denv1', help='Name of the SnpEff database (default: denv1).')
    args = parser.parse_args()

    # Validate input paths
    check_file_exists(args.input_dir, "Input directory")
    check_file_exists(args.reference_fasta, "Reference FASTA file")
    check_file_exists(args.genbank_file, "GenBank file")

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    sam_files_dir = os.path.join(args.output_dir, "sam_files")
    os.makedirs(sam_files_dir, exist_ok=True)

    # Define paths for intermediate files
    sample_sheet = os.path.join(args.output_dir, "samplesheet.tsv")

    # Step 1: Generate sample sheet
    cmd1 = f"python 1.create_samplesheet.py {args.input_dir} {sample_sheet}"
    run_command(cmd1, "Sample sheet generation")

    # Verify sample sheet exists
    check_file_exists(sample_sheet, "Sample sheet")

    # Step 2: Map reads
    cmd2 = f"python 2.map_reads.py --samplesheet {sample_sheet} --reference {args.reference_fasta}"
    run_command(cmd2, "Read mapping")

    # Step 3: Convert SAM to BAM
    cmd3 = f"python 3.samtobamdenv.py --input_dir {sam_files_dir} --reference_fasta {args.reference_fasta} --output_dir {args.output_dir}"
    run_command(cmd3, "SAM to BAM conversion")

    # Step 4: Create SnpEff database
    cmd4 = f"python 5.create_snpeff_database.py --genbank_file {args.genbank_file} --reference_fasta {args.reference_fasta} --output_dir {args.output_dir} --database_name {args.database_name}"
    run_command(cmd4, "SnpEff database creation")

    # Step 5: Generate consensus sequences and initial VCFs
    cmd5 = f"python 4.sam2consensus_test2_ivar.py --input_dir {args.output_dir} --reference_fasta {args.reference_fasta} --output_dir {args.output_dir}"
    run_command(cmd5, "Consensus sequence generation")

    # Step 6: Variant calling and SnpEff annotation
    cmd6 = f"python 6.variant_calling_consensus.py --input_dir {args.output_dir} --reference_fasta {args.reference_fasta} --output_dir {args.output_dir} --database_name {args.database_name}"
    run_command(cmd6, "Variant calling and annotation")

    # Step 7: Summarize coverage and FASTA results
    cmd7 = f"python 7.summarize_result.py --input_dir {args.output_dir} --output_dir {args.output_dir}"
    run_command(cmd7, "Coverage and FASTA summarization")

    # Step 8: Summarize SnpEff annotations
    cmd8 = f"python 8.summarize_snpEff.py"
    run_command(cmd8, "SnpEff annotation summarization")

    logging.info("Pipeline completed successfully!")

if __name__ == "__main__":
    main()