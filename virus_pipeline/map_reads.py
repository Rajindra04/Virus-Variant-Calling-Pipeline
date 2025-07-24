import sys
import argparse
import pandas as pd
import subprocess
import os
import logging

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_command(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise Exception(f"Command execution failed with return code {process.returncode}, stderr: {stderr.decode('utf-8')}")
    return stdout.decode("utf-8"), stderr.decode("utf-8")

def samplesheet_verify(samplesheet):
    try:
        samples = pd.read_csv(samplesheet, sep='\t')
        if not all(col in samples.columns for col in ['sample_name', 'read1', 'read2']):
            raise ValueError("Sample sheet must have columns 'sample_name', 'read1', 'read2'.")
    except Exception as e:
        raise ValueError(f"Error reading sample sheet: {e}")
    return samples

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--samplesheet', type=str, help='Path to sample sheet with columns "sample_name", "read1", "read2".')
    parser.add_argument('--reference', type=str, help='Path to reference FASTA file')
    args = parser.parse_args()

    samples = samplesheet_verify(args.samplesheet)
    base_dir = os.path.abspath(os.path.dirname(args.samplesheet))
    sam_files_dir = os.path.join(base_dir, "sam_files")
    os.makedirs(sam_files_dir, exist_ok=True)

    # Reference indexing with bwa-mem2
    index_command = f"bwa-mem2 index {args.reference}"
    print(f"Indexing reference {args.reference}...")
    stdout, stderr = run_command(index_command)
    print(stdout)
    print(stderr)
          
    for i, row in samples.iterrows():
        sample_name = row['sample_name']
        read1 = row['read1']
        read2 = row['read2']

        sample_folder = os.path.join(base_dir, f"{sample_name}_output")
        os.makedirs(sample_folder, exist_ok=True)

        # fastp
        trimmomatic_command = f"fastp -i {read1} -I {read2} -o {os.path.join(sample_folder, f'{sample_name}_trimmed_R1.fastq.gz')} -O {os.path.join(sample_folder, f'{sample_name}_trimmed_R2.fastq.gz')}"
        print(f"Running fastp for sample {sample_name}...")
        stdout, stderr = run_command(trimmomatic_command)
        print(stdout)
        print(stderr)  

        # FastQC
        fastqc_command = f"fastqc {os.path.join(sample_folder, f'{sample_name}_trimmed_R1.fastq.gz')} {os.path.join(sample_folder, f'{sample_name}_trimmed_R2.fastq.gz')} --outdir={sample_folder}"
        print(f"Running FastQC for sample {sample_name}...")
        stdout, stderr = run_command(fastqc_command)
        print(stdout)
        print(stderr)
        
        # Mapping with bwa-mem2
        sam_file = os.path.join(sample_folder, f"{sample_name}_aln.sam")
        bwa_command = f"bwa-mem2 mem {args.reference} {os.path.join(sample_folder, f'{sample_name}_trimmed_R1.fastq.gz')} {os.path.join(sample_folder, f'{sample_name}_trimmed_R2.fastq.gz')} > {sam_file}"
        print(f"Running bwa-mem2 for sample {sample_name}...")
        stdout, stderr = run_command(bwa_command)
        print(stdout)
        print(stderr)

        # Move SAM files to sam_files_dir
        run_command(f"mv {os.path.join(sample_folder, f'{sample_name}*.sam')} {sam_files_dir}")
        print(f"SAM files moved to {sam_files_dir}")
        
    logging.info("Processing completed!")