#!/usr/bin/env python

import argparse
import os
import sys
import logging

# Add the pipeline module to the Python path if running locally
sys.path.append(os.path.join(os.path.dirname(__file__), 'pipeline'))

from pipeline import (
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

def main():
    parser = argparse.ArgumentParser(description='Automate variant calling pipeline.')
    parser.add_argument('--input_dir', type=str, required=True)
    parser.add_argument('--reference_fasta', type=str, required=True)
    parser.add_argument('--genbank_file', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)
    parser.add_argument('--database_name', type=str, default='denv1')
    args = parser.parse_args()

    check_file_exists(args.input_dir, "Input directory")
    check_file_exists(args.reference_fasta, "Reference FASTA")
    check_file_exists(args.genbank_file, "GenBank file")

    os.makedirs(args.output_dir, exist_ok=True)
    sam_files_dir = os.path.join(args.output_dir, "sam_files")
    os.makedirs(sam_files_dir, exist_ok=True)

    sample_sheet = os.path.join(args.output_dir, "samplesheet.tsv")

    # Run steps by calling modules directly
    create_samplesheet.main([args.input_dir, sample_sheet])
    check_file_exists(sample_sheet, "Sample sheet")

    map_reads.main([
        '--samplesheet', sample_sheet,
        '--reference', args.reference_fasta
    ])

    samtobamdenv.main([
        '--input_dir', sam_files_dir,
        '--reference_fasta', args.reference_fasta,
        '--output_dir', args.output_dir
    ])

    create_snpeff_database.main([
        '--genbank_file', args.genbank_file,
        '--reference_fasta', args.reference_fasta,
        '--output_dir', args.output_dir,
        '--database_name', args.database_name
    ])

    sam2consensus_test2_ivar.main([
        '--input_dir', args.output_dir,
        '--reference_fasta', args.reference_fasta,
        '--output_dir', args.output_dir
    ])

    variant_calling_consensus.main([
        '--input_dir', args.output_dir,
        '--reference_fasta', args.reference_fasta,
        '--output_dir', args.output_dir,
        '--database_name', args.database_name
    ])

    summarize_result.main([
        '--input_dir', args.output_dir,
        '--output_dir', args.output_dir
    ])

    summarize_snpEff.main([])

    logging.info("Pipeline completed successfully!")

if __name__ == "__main__":
    main()

