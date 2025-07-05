import argparse
import glob
import os
import subprocess
import matplotlib.pyplot as plt
import re

# Usage: python 6.variant_calling_consensus.py --input_dir test/ --reference_fasta references/denv1.fasta --output_dir test_out/ --database_name denv1

def run_command(command):
    print("Running command:", command)
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    process.wait()

    stdout = process.stdout.read().decode("utf-8")
    stderr = process.stderr.read().decode("utf-8")

    if process.returncode != 0:
        error_message = stderr.strip()
        raise Exception(f"Command execution failed with return code {process.returncode}: {error_message}")

    return stdout, stderr

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
    print(f"FASTA file {fasta_file} validated: contains only IUPAC bases.")

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

    print("Coverage plot generated:", coverage_plot_file)

def prepare_reference(reference_fasta, output_dir):
    validate_fasta(reference_fasta)
    if not os.path.exists(f"{reference_fasta}.fai"):
        print("Indexing reference FASTA...")
        run_command(f"samtools faidx {reference_fasta}")

    dict_file = os.path.splitext(reference_fasta)[0] + ".dict"
    if not os.path.exists(dict_file):
        print("Creating sequence dictionary for reference...")
        run_command(f"gatk CreateSequenceDictionary -R {reference_fasta} -O {dict_file}")

def add_read_groups(bam_file, sample_name, output_dir):
    rg_bam = os.path.join(output_dir, f"{sample_name}_rg.bam")
    rg_command = (
        f"samtools addreplacerg -r 'ID:{sample_name}\tSM:{sample_name}\tLB:lib1\tPL:ILLUMINA' "
        f"{bam_file} -o {rg_bam}"
    )
    stdout, stderr = run_command(rg_command)
    print("Add read groups output:", stdout)
    print("Add read groups error:", stderr)

    run_command(f"samtools index {rg_bam}")
    return rg_bam

def run_variant_calling(bam_file, reference_fasta, sample_name, output_dir):
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
    print("GATK HaplotypeCaller output:", stdout)
    print("GATK HaplotypeCaller error:", stderr)
    return raw_vcf

def run_snpeff_annotation(raw_vcf, sample_name, output_dir, reference_name="denv1"):
    annotated_vcf = os.path.join(output_dir, f"{sample_name}_annotated.vcf")
    summary_html = os.path.join(output_dir, f"{sample_name}_snpEff_summary.html")    
    snpeff_command = (
        f"snpEff -Xmx4g "
        f"-c {os.path.join(output_dir, 'snpEff.config')} "
        f"-v {reference_name} "
        f"-s {summary_html} "        
        f"{raw_vcf} > {annotated_vcf}"
    )
    stdout, stderr = run_command(snpeff_command)
    print("SnpEff annotation output:", stdout)
    print("SnpEff annotation error:", stderr)
    return annotated_vcf, summary_html, os.path.join(output_dir, f"{sample_name}_snpEff_genes.txt")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, help='Path to input directory containing BAM files.')
    parser.add_argument('--reference_fasta', type=str, help='Path to reference FASTA file.')
    parser.add_argument('--output_dir', type=str, help='Path to output directory for consensus FASTA and VCF files.')
    parser.add_argument('--database_name', type=str, default="denv1", help='Name of the SnpEff database (default: denv1).')
    args = parser.parse_args()

    input_dir = args.input_dir
    reference_fasta = args.reference_fasta
    output_dir = args.output_dir
    database_name = args.database_name
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    prepare_reference(reference_fasta, output_dir)
    bam_files = glob.glob(os.path.join(input_dir, "*sorted.bam"))

    for bam_file in bam_files:
        print("Processing:", bam_file)
        sample_name = os.path.splitext(os.path.basename(bam_file))[0]

        try:
            rg_bam = add_read_groups(bam_file, sample_name, output_dir)
            coverage_file = f"{sample_name}_coverage.txt"
            coverage_qual_file = f"{sample_name}_qual.txt"
            coverage_command = f"samtools depth {rg_bam} > {os.path.join(output_dir, coverage_file)} && samtools depth -Q 20 {rg_bam} > {os.path.join(output_dir, coverage_qual_file)}"
            stdout, stderr = run_command(coverage_command)
            print("Command output:", stdout)
            print("Command error:", stderr)
            plot_coverage(os.path.join(output_dir, coverage_file), output_dir)
            pileup_command = f"samtools mpileup -d 1000 -A -Q 0 {rg_bam} | ivar consensus -p {sample_name} -q 20 -t 0"
            stdout, stderr = run_command(pileup_command)
            print("Command output:", stdout)
            print("Command error:", stderr)
            consensus_fasta = f"{sample_name}.fa"
            os.rename(consensus_fasta, os.path.join(output_dir, consensus_fasta))
            raw_vcf = run_variant_calling(rg_bam, reference_fasta, sample_name, output_dir)
            annotated_vcf, summary_html, summary_txt = run_snpeff_annotation(raw_vcf, sample_name, output_dir, database_name)
            print("Consensus sequence, coverage file, coverage plot, annotated VCF, and SnpEff reports generated:", 
                  os.path.join(output_dir, consensus_fasta), 
                  os.path.join(output_dir, coverage_file),
                  annotated_vcf,  # Unpacked tuple
                  summary_html,
                  summary_txt)
        except Exception as e:
            print("Error occurred during processing:", str(e))
            continue
