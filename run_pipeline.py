#!/usr/bin/env python

import argparse
import os
import sys
import logging
import shutil

from virus_pipeline import (
    create_samplesheet,
    map_reads,
    samtobamdenv,
    create_snpeff_database,
    sam2consensus_test2_ivar,
    variant_calling_consensus,
    summarize_result,
    summarize_snpEff,
)

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def check_file_exists(file_path, description):
    if not os.path.exists(file_path):
        logging.error(f"{description} not found: {file_path}")
        sys.exit(1)

def check_write_permission(directory):
    try:
        test_file = os.path.join(directory, '.test_write')
        with open(test_file, 'w') as f:
            f.write('test')
        os.remove(test_file)
    except (PermissionError, OSError) as e:
        logging.error(f"No write permission for {directory}: {e}")
        sys.exit(1)

def check_tools():
    tools = ['bwa-mem2', 'samtools', 'fastp', 'fastqc', 'gatk', 'ivar', 'bcftools']
    # Check for both snpeff and snpEff
    snpeff_found = False
    for tool in ['snpeff', 'snpEff']:
        if shutil.which(tool):
            snpeff_found = True
            break
    if not snpeff_found:
        logging.error("Required tool not found: snpeff/snpEff")
        sys.exit(1)
    for tool in tools:
        if not shutil.which(tool):
            logging.error(f"Required tool not found: {tool}")
            sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Automate variant calling pipeline.')
    parser.add_argument('--input_dir', type=str, required=True)
    parser.add_argument('--reference_fasta', type=str, required=True)
    parser.add_argument('--genbank_file', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)
    parser.add_argument('--database_name', type=str, default='denv1')
    args = parser.parse_args()

    # Check inputs and tools
    check_file_exists(args.input_dir, "Input directory")
    check_file_exists(args.reference_fasta, "Reference FASTA")
    check_file_exists(args.genbank_file, "GenBank file")
    check_tools()

    # Create and check output directories
    os.makedirs(args.output_dir, exist_ok=True)
    check_write_permission(args.output_dir)
    sam_files_dir = os.path.join(args.output_dir, "sam_files")
    os.makedirs(sam_files_dir, exist_ok=True)
    check_write_permission(sam_files_dir)

    sample_sheet = os.path.join(args.output_dir, "samplesheet.tsv")

    # Run pipeline steps
    try:
        logging.info("Starting create_samplesheet")
        create_samplesheet([args.input_dir, sample_sheet])
        check_file_exists(sample_sheet, "Sample sheet")
        logging.info("Completed create_samplesheet")
    except Exception as e:
        logging.error(f"Failed in create_samplesheet: {e}")
        sys.exit(1)

    try:
        logging.info("Starting map_reads")
        map_reads(['--samplesheet', sample_sheet, '--reference', args.reference_fasta])
        logging.info("Completed map_reads")
    except Exception as e:
        logging.error(f"Failed in map_reads: {e}")
        sys.exit(1)

    try:
        logging.info("Starting samtobamdenv")
        samtobamdenv(['--input_dir', sam_files_dir, '--reference_fasta', args.reference_fasta, '--output_dir', args.output_dir])
        logging.info("Completed samtobamdenv")
    except Exception as e:
        logging.error(f"Failed in samtobamdenv: {e}")
        sys.exit(1)

    try:
        logging.info("Starting create_snpeff_database")
        create_snpeff_database(['--genbank_file', args.genbank_file, '--reference_fasta', args.reference_fasta, '--output_dir', args.output_dir, '--database_name', args.database_name])
        logging.info("Completed create_snpeff_database")
    except Exception as e:
        logging.error(f"Failed in create_snpeff_database: {e}")
        sys.exit(1)

    try:
        logging.info("Starting sam2consensus_test2_ivar")
        sam2consensus_test2_ivar(['--input_dir', args.output_dir, '--reference_fasta', args.reference_fasta, '--output_dir', args.output_dir])
        logging.info("Completed sam2consensus_test2_ivar")
    except Exception as e:
        logging.error(f"Failed in sam2consensus_test2_ivar: {e}")
        sys.exit(1)

    try:
        logging.info("Starting variant_calling_consensus")
        variant_calling_consensus(['--input_dir', args.output_dir, '--reference_fasta', args.reference_fasta, '--output_dir', args.output_dir, '--database_name', args.database_name])
        logging.info("Completed variant_calling_consensus")
    except Exception as e:
        logging.error(f"Failed in variant_calling_consensus: {e}")
        sys.exit(1)

    try:
        logging.info("Starting summarize_result")
        summarize_result(['--input_dir', args.output_dir, '--output_dir', args.output_dir])
        logging.info("Completed summarize_result")
    except Exception as e:
        logging.error(f"Failed in summarize_result: {e}")
        sys.exit(1)

    try:
        logging.info("Starting summarize_snpEff")
        summarize_snpEff(['--input_dir', args.output_dir, '--output_dir', args.output_dir])
        logging.info("Completed summarize_snpEff")
    except Exception as e:
        logging.error(f"Failed in summarize_snpEff: {e}")
        sys.exit(1)

    logging.info("Pipeline completed successfully!")

if __name__ == "__main__":
    main()
