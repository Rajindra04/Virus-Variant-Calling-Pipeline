import sys
import argparse
import glob
import os
import subprocess
import matplotlib.pyplot as plt
import re
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_command(command):
    logging.info(f"Running command: {command}")
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise Exception(f"Command execution failed with return code {process.returncode}, stderr: {stderr.decode('utf-8')}")
    return stdout.decode("utf-8"), stderr.decode("utf-8")

def validate_fasta(fasta_file):
    valid_bases = set("ACGTNacgtnRYKMSWBDHVrykmswbvh")
    with open(fasta_file, "r") as f:
        sequence = ""
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip()
    invalid_bases = set(sequence) - valid_bases
    if invalid_bases:
        raise ValueError(f"Non-IUPAC bases found in {fasta_file}: {invalid_bases}")
    logging.info(f"FASTA file {fasta_file} validated: contains only IUPAC bases.")

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

def plot_coverage(coverage_file, output_dir):
    positions = []
    depths = []
    with open(coverage_file, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            position = int(columns[1])
            depth = int(columns[2])
            positions.append(position)
            depths.append(depth)

    plt.plot(positions, depths)
    plt.xlabel('Position')
    plt.ylabel('Coverage')
    plt.title('Read Coverage')
    plt.grid(True)

    coverage_plot_file = os.path.join(output_dir, os.path.splitext(os.path.basename(coverage_file))[0] + '.png')
    plt.savefig(coverage_plot_file)
    plt.close()

    logging.info(f"Coverage plot generated: {coverage_plot_file}")

def prepare_reference(reference_fasta, output_dir):
    validate_fasta(reference_fasta)
    if not os.path.exists(f"{reference_fasta}.fai"):
        logging.info("Indexing reference FASTA...")
        run_command(f"samtools faidx {reference_fasta}")

    dict_file = os.path.splitext(reference_fasta)[0] + ".dict"
    if not os.path.exists(dict_file):
        logging.info("Creating sequence dictionary for reference...")
        run_command(f"gatk CreateSequenceDictionary -R {reference_fasta} -O {dict_file}")

def add_read_groups(bam_file, sample_name, output_dir):
    rg_bam = os.path.join(output_dir, f"{sample_name}_rg.bam")
    rg_command = (
        f"samtools addreplacerg -r 'ID:{sample_name}\tSM:{sample_name}\tLB:lib1\tPL:ILLUMINA' "
        f"-o {rg_bam} {bam_file}"
    )
    stdout, stderr = run_command(rg_command)
    print("Add read groups output:", stdout)
    print("Add read groups error:", stderr)

    # Validate and sort the BAM file
    run_command(f"samtools sort -T {os.path.join(output_dir, f'temp_{sample_name}')} -o {rg_bam}.sorted {rg_bam}")
    os.rename(f"{rg_bam}.sorted", rg_bam)
    run_command(f"samtools index {rg_bam}")
    return rg_bam

def run_variant_calling(bam_file, reference_fasta, sample_name, output_dir):
    validate_bam(bam_file)
    raw_vcf = os.path.join(output_dir, f"{sample_name}.vcf")
    gatk_command = (
        f"gatk HaplotypeCaller "
        f"-R {reference_fasta} "
        f"-I {bam_file} "
        f"-O {raw_vcf} "
        f"--standard-min-confidence-threshold-for-calling 30 "
        f"--min-base-quality-score 20"
    )
    stdout, stderr = run_command(gatk_command)
    logging.info(f"GATK HaplotypeCaller output: {stdout}")
    logging.info(f"GATK HaplotypeCaller error: {stderr}")
    return raw_vcf

def run_snpeff_annotation(raw_vcf, sample_name, output_dir, reference_name="denv1"):
    annotated_vcf = os.path.join(output_dir, f"{sample_name}_annotated.vcf")
    summary_html = os.path.join(output_dir, f"{sample_name}_snpEff_summary.html")
    summary_csv = os.path.join(output_dir, f"{sample_name}_snpEff_summary.csv")
    snpeff_command = (
        f"snpEff -Xmx4g "
        f"-c {os.path.join(output_dir, 'snpEff.config')} "
        f"-v {reference_name} "
        f"-s {summary_html} "
        f"-csvStats {summary_csv} "
        f"{raw_vcf} > {annotated_vcf}"
    )
    stdout, stderr = run_command(snpeff_command)
    logging.info(f"SnpEff annotation output: {stdout}")
    logging.info(f"SnpEff annotation error: {stderr}")
    return annotated_vcf, summary_html, summary_csv, os.path.join(output_dir, f"{sample_name}_snpEff_genes.txt")

def run_snpsift_extract(annotated_vcf, sample_name, output_dir):
    snpsift_output = os.path.join(output_dir, f"{sample_name}_snpSift.txt")
    snpsift_command = (
        f"SnpSift extractFields {annotated_vcf} "
        f"CHROM POS REF ALT "
        f"\"ANN[*].EFFECT\" \"ANN[*].IMPACT\" \"ANN[*].GENE\" \"ANN[*].GENEID\" "
        f"\"ANN[*].FEATURE\" \"ANN[*].HGVS_C\" \"ANN[*].HGVS_P\" \"ANN[*].AA_POS\" "
        f"\"EFF[*].CODON\" \"EFF[*].AA\" \"EFF[*].GENE\" > {snpsift_output}"
    )
    stdout, stderr = run_command(snpsift_command)
    logging.info(f"SnpSift extract output: {stdout}")
    logging.info(f"SnpSift extract error: {stderr}")
    return snpsift_output

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, help='Path to input directory containing BAM files.')
    parser.add_argument('--reference_fasta', type=str, help='Path to reference FASTA file.')
    parser.add_argument('--output_dir', type=str, help='Path to output directory for consensus FASTA and VCF files.')
    parser.add_argument('--database_name', type=str, default="denv1", help='Name of the SnpEff database (default: denv1).')
    args = parser.parse_args(argv)

    input_dir = args.input_dir
    reference_fasta = args.reference_fasta
    output_dir = args.output_dir
    database_name = args.database_name
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    bam_files = glob.glob(os.path.join(input_dir, "*.sorted.bam"))
    if not bam_files:
        logging.error(f"No sorted BAM files found in {input_dir}")
        sys.exit(1)

    prepare_reference(reference_fasta, output_dir)
    bam_files = glob.glob(os.path.join(input_dir, "*sorted.bam"))
    for bam_file in bam_files:
        print("Processing:", bam_file)
        sample_name = os.path.splitext(os.path.basename(bam_file))[0].replace('.sorted', '')
        try:
            rg_bam = add_read_groups(bam_file, sample_name, output_dir)

        
            # Add read groups before GATK
            rg_bam = add_read_groups(bam_file, sample_name, output_dir)
            coverage_file = os.path.join(output_dir, f"{sample_name}_coverage.txt")
            coverage_qual_file = os.path.join(output_dir, f"{sample_name}_qual.txt")
            coverage_command = f"samtools depth {rg_bam} > {coverage_file} && samtools depth -Q 20 {rg_bam} > {coverage_qual_file}"
            stdout, stderr = run_command(coverage_command)
            logging.info(f"Coverage command output: {stdout}")
            logging.info(f"Coverage command error: {stderr}")
            plot_coverage(coverage_file, output_dir)
            pileup_command = f"samtools mpileup -d 1000 -A -Q 0 {rg_bam} | ivar consensus -p {sample_name} -q 20 -t 0"
            stdout, stderr = run_command(pileup_command)
            logging.info(f"Consensus command output: {stdout}")
            logging.info(f"Consensus command error: {stderr}")
            consensus_fasta = f"{sample_name}.fa"
            os.rename(consensus_fasta, os.path.join(output_dir, consensus_fasta))
            raw_vcf = run_variant_calling(rg_bam, reference_fasta, sample_name, output_dir)
            annotated_vcf, summary_html, summary_csv, summary_txt = run_snpeff_annotation(raw_vcf, sample_name, output_dir, database_name)
            snpsift_output = run_snpsift_extract(annotated_vcf, sample_name, output_dir)
            logging.info(f"Consensus sequence, coverage file, coverage plot, annotated VCF, SnpEff reports, and SnpSift output: {os.path.join(output_dir, consensus_fasta)}, {coverage_file}, {raw_vcf}, {annotated_vcf}, {summary_html}, {summary_csv}, {summary_txt}, {snpsift_output}")
        except Exception as e:
            print(f"Error occurred during processing {sample_name}: {str(e)}")
            continue

if __name__ == '__main__':
    main()
