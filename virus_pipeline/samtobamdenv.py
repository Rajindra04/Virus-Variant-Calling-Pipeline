import sys
import argparse
import glob
import os
import subprocess
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_command(command):
    logging.info(f"Running command: {command}")
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise Exception(f"Command execution failed with return code {process.returncode}, stderr: {stderr.decode('utf-8')}")
    return stdout.decode("utf-8"), stderr.decode("utf-8")

def validate_bam(bam_file):
    """Validate BAM file using samtools quickcheck and check alignment count."""
    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"BAM file {bam_file} does not exist")
    validate_command = f"samtools quickcheck {bam_file}"
    stdout, stderr = run_command(validate_command)
    if stderr:
        raise Exception(f"BAM validation failed for {bam_file}: {stderr}")
    # Check alignment count
    count_command = f"samtools view -c {bam_file}"
    stdout, stderr = run_command(count_command)
    alignment_count = int(stdout.strip())
    if alignment_count == 0:
        raise ValueError(f"BAM file {bam_file} contains no alignments")
    logging.info(f"BAM file {bam_file} validated successfully with {alignment_count} alignments")

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, help='Path to input directory containing SAM files.')
    parser.add_argument('--reference_fasta', type=str, help='Path to reference FASTA file.')
    parser.add_argument('--output_dir', type=str, help='Path to output directory for BAM files.')
    args = parser.parse_args(argv)

    input_dir = args.input_dir
    reference_fasta = args.reference_fasta
    output_dir = args.output_dir

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    sam_files = glob.glob(os.path.join(input_dir, "*_aln.sam"))
    if not sam_files:
        logging.error(f"No SAM files found in {input_dir}")
        sys.exit(1)

    for sam_file in sam_files:
        logging.info(f"Processing: {sam_file}")
        sample_name = os.path.basename(sam_file).replace('_aln.sam', '')

        try:
            # Convert SAM to BAM
            bam_file = os.path.join(output_dir, f"{sample_name}.bam")
            sam_to_bam_command = f"samtools view -S -b {sam_file} > {bam_file}"
            stdout, stderr = run_command(sam_to_bam_command)
            validate_bam(bam_file)

            # Sort BAM file
            sorted_bam_file = os.path.join(output_dir, f"{sample_name}.sorted.bam")
            sort_command = f"samtools sort -T {os.path.join(output_dir, f'temp_{sample_name}')} -o {sorted_bam_file} {bam_file}"
            stdout, stderr = run_command(sort_command)
            validate_bam(sorted_bam_file)

            # Index sorted BAM file
            index_command = f"samtools index {sorted_bam_file}"
            stdout, stderr = run_command(index_command)
            logging.info(f"Sorted BAM file generated and indexed: {sorted_bam_file}")

            # Clean up intermediate BAM file
            os.remove(bam_file)
        except Exception as e:
            logging.error(f"Error processing {sam_file}: {str(e)}")
            continue

if __name__ == '__main__':
    main()