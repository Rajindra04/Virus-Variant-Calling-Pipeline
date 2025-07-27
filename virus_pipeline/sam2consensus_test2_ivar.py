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

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, help='Path to input directory containing BAM files.')
    parser.add_argument('--reference_fasta', type=str, help='Path to reference FASTA file.')
    parser.add_argument('--output_dir', type=str, help='Path to output directory for VCF files.')
    args = parser.parse_args(argv)

    input_dir = args.input_dir
    reference_fasta = args.reference_fasta
    output_dir = args.output_dir

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    bam_files = glob.glob(os.path.join(input_dir, "*.sorted.bam"))
    if not bam_files:
        logging.error(f"No sorted BAM files found in {input_dir}")
        sys.exit(1)

    for bam_file in bam_files:
        logging.info(f"Processing: {bam_file}")
        sample_name = os.path.basename(bam_file).replace('.sorted.bam', '')

        try:
            # Ensure reference FASTA has an index
            if not os.path.exists(f"{reference_fasta}.fai"):
                logging.info(f"Indexing reference FASTA: {reference_fasta}")
                run_command(f"samtools faidx {reference_fasta}")

            # Variant calling with samtools + bcftools
            vcf_gz_file = os.path.join(output_dir, f"{sample_name}.vcf.gz")
            vcf_file = os.path.join(output_dir, f"{sample_name}.vcf")
            bcftools_command = (
                f"samtools mpileup -f {reference_fasta} {bam_file} | "
                f"bcftools call -mv --ploidy 1 -Oz -o {vcf_gz_file}"
            )

            stdout, stderr = run_command(bcftools_command)
            logging.info(f"Variant calling output: {stdout}")
            logging.info(f"Variant calling error: {stderr}")

            # Index the compressed VCF
            run_command(f"bcftools index {vcf_gz_file}")

            # Convert to plain .vcf
            run_command(f"bcftools view {vcf_gz_file} > {vcf_file}")

            # Clean up .vcf.gz
            if os.path.exists(vcf_file):
                os.remove(vcf_gz_file)
                logging.info(f"VCF file created: {vcf_file}")
            else:
                logging.error(f"VCF file not found: {vcf_file}")
                sys.exit(1)

        except Exception as e:
            logging.error(f"Error occurred during processing {sample_name}: {str(e)}")
            continue

if __name__ == '__main__':
    main()
