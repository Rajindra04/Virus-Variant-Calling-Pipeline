import sys

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
import argparse
import glob
import os
import subprocess
import matplotlib.pyplot as plt

# Usage: python sam_to_consensus.py --input_dir input_dir --reference_fasta reference.fasta --output_dir output_dir 

def run_command(command):
    print("Running command:", command)
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    process.wait()

    stdout = process.stdout.read().decode("utf-8")
    stderr = process.stderr.read().decode("utf-8")

    if process.returncode != 0:
        error_message = stderr.strip()
        raise Exception("Command execution failed with return code {}: {}".format(process.returncode, error_message))

    return stdout, stderr

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

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, help='Path to input directory containing BAM files.')
    parser.add_argument('--reference_fasta', type=str, help='Path to reference FASTA file.')
    parser.add_argument('--output_dir', type=str, help='Path to output directory for consensus FASTA files and VCF.')
    args = parser.parse_args()

    input_dir = args.input_dir
    reference_fasta = args.reference_fasta
    output_dir = args.output_dir

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    bam_files = glob.glob(os.path.join(input_dir, "*sorted.bam"))

    for bam_file in bam_files:
        print("Processing:", bam_file)

        sample_name = os.path.splitext(os.path.basename(bam_file))[0]

        try:
            # Generate coverage files
            coverage_file = f"{sample_name}_coverage.txt"
            coverage_qual_file = f"{sample_name}_qual.txt"
            coverage_command = f"samtools depth {bam_file} > {os.path.join(output_dir, coverage_file)} && samtools depth -Q 20 {bam_file} > {os.path.join(output_dir, coverage_qual_file)}"
            stdout, stderr = run_command(coverage_command)
            print("Command output:", stdout)
            print("Command error:", stderr)

            # Plot coverage
            plot_coverage(os.path.join(output_dir, coverage_file), output_dir)

            # Generate pileup for variant calling
            pileup_file = os.path.join(output_dir, f"{sample_name}.pileup")
            mpileup_command = f"samtools mpileup -d 1000 -A -Q 0 -f {reference_fasta} {bam_file} > {pileup_file}"
            stdout, stderr = run_command(mpileup_command)
            print("Pileup created:", pileup_file)

            # Variant calling with ivar
            # Define output VCF path (.vcf.gz initially)
            # Ensure the reference FASTA has an index (.fai)
            if not os.path.exists(reference_fasta + ".fai"):
                print("Indexing reference FASTA...")
                index_command = f"samtools faidx {reference_fasta}"
                stdout, stderr = run_command(index_command)

            # Define output file paths
            vcf_gz_file = os.path.join(output_dir, sample_name + ".vcf.gz")
            vcf_file = os.path.join(output_dir, sample_name + ".vcf")

            # Variant calling using samtools + bcftools
            bcftools_command = (
                f"samtools mpileup -f {reference_fasta} {bam_file} | "
                f"bcftools call -mv --ploidy 1 -Oz -o {vcf_gz_file}"
            )

            try:
                # Run the variant calling pipeline
                stdout, stderr = run_command(bcftools_command)

                # Index the compressed VCF
                index_command = f"bcftools index {vcf_gz_file}"
                stdout, stderr = run_command(index_command)

                # Convert to plain .vcf for compatibility
                decompress_command = f"bcftools view {vcf_gz_file} > {vcf_file}"
                stdout, stderr = run_command(decompress_command)

                # Confirm output
                if os.path.exists(vcf_file):
                    print("✅ VCF file created:", vcf_file)
                else:
                    print("⚠️ VCF file not found (but no errors occurred):", vcf_file)

            except Exception as e:
                print("❌ Error occurred during bcftools variant calling:", str(e))

            # Generate consensus sequence
            ivar_consensus_command = f"samtools mpileup -d 1000 -A -Q 0 -f {reference_fasta} {bam_file} | ivar consensus -p {sample_name} -q 20 -t 0.75"
            stdout, stderr = run_command(ivar_consensus_command)
            print("Consensus calling complete.")

            # Move consensus file to output
            consensus_fasta = f"{sample_name}.fa"
            if os.path.exists(consensus_fasta):
                os.rename(consensus_fasta, os.path.join(output_dir, consensus_fasta))

            print("Generated:", os.path.join(output_dir, consensus_fasta), os.path.join(output_dir, vcf_file))

        except Exception as e:
            print("Error occurred during processing:", str(e))
            continue

