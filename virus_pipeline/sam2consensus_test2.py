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

    bam_files = glob.glob(os.path.join(input_dir, "*_rg.bam"))
    if not bam_files:
        logging.error(f"No read-group BAM files found in {input_dir}")
        sys.exit(1)

    for bam_file in bam_files:
        logging.info(f"Processing: {bam_file}")
        sample_name = os.path.basename(bam_file).replace('_rg.bam', '')
        
        try:
            vcf_output = os.path.join(output_dir, f"{sample_name}.vcf.gz")
            pileup_command = (
                f"samtools mpileup -f {reference_fasta} {bam_file} | "
                f"bcftools call -mv --ploidy 1 -Oz -o {vcf_output}"
            )
            stdout, stderr = run_command(pileup_command)
            logging.info(f"Variant-calling output: {stdout}")
            logging.info(f"Variant-calling error: {stderr}")
        except Exception as e:
            logging.error(f"Error occurred during processing {sample_name}: {str(e)}")
            continue

if __name__ == '__main__':
    main()