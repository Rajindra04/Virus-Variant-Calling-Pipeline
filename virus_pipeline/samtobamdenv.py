import sys
import argparse
import glob
import os
import subprocess

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

def run_command(command):
    print(f"Running command: {command}")
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise Exception(f"Command execution failed with return code {process.returncode}, stderr: {stderr.decode('utf-8')}")
    return stdout.decode("utf-8"), stderr.decode("utf-8")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, help='Path to input directory containing SAM files.')
    parser.add_argument('--reference_fasta', type=str, help='Path to reference FASTA file.')
    parser.add_argument('--output_dir', type=str, help='Path to output directory for BAM files.')
    args = parser.parse_args()

    input_dir = args.input_dir
    reference_fasta = args.reference_fasta
    output_dir = args.output_dir

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    sam_files = glob.glob(os.path.join(input_dir, "*_aln.sam"))

    for sam_file in sam_files:
        print(f"Processing: {sam_file}")
        sample_name = os.path.basename(sam_file).replace('_aln.sam', '')

        # Convert SAM to BAM
        bam_file = os.path.join(output_dir, f"{sample_name}.bam")
        sam_to_bam_command = f"samtools view -S -b {sam_file} > {bam_file}"
        run_command(sam_to_bam_command)

        # Sort BAM file
        sorted_bam_file = os.path.join(output_dir, f"{sample_name}.sorted.bam")
        sort_command = f"samtools sort -T {os.path.join(output_dir, f'temp_{sample_name}')} -o {sorted_bam_file} {bam_file}"
        run_command(sort_command)

        # Index BAM file
        index_command = f"samtools index {sorted_bam_file}"
        run_command(index_command)

        print(f"BAM file generated: {sorted_bam_file}")