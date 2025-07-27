import sys
import os
import argparse
import pandas as pd
import logging
pd.set_option('future.no_silent_downcasting', True)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def summarize_coverage(input_dir, output_file):
    sample_names = []
    average_coverage = []

    for file_name in os.listdir(input_dir):
        if file_name.endswith('_coverage.txt'):
            coverage_file = os.path.join(input_dir, file_name)
            sample_name = file_name.replace('_coverage.txt', '')
            sample_names.append(sample_name)

            coverage_values = []
            try:
                with open(coverage_file, 'r') as file:
                    for line in file:
                        columns = line.strip().split('\t')
                        depth = int(columns[2])
                        coverage_values.append(depth)
                if coverage_values:
                    average_coverage.append(sum(coverage_values) / len(coverage_values))
                else:
                    average_coverage.append(0)
            except Exception as e:
                logging.error(f"Error processing coverage file {coverage_file}: {e}")
                average_coverage.append(0)

    coverage_data = {
        'Sample': sample_names,
        'Average Coverage': average_coverage,
    }
    coverage_df = pd.DataFrame(coverage_data)
    coverage_df.to_excel(output_file, index=False)
    logging.info(f"Coverage summary saved to {output_file}")

def summarize_fasta(input_dir, output_file, database_name):
    sample_names = []
    n_counts = []
    total_bases = []

    for file_name in os.listdir(input_dir):
        if file_name.endswith(f'_{database_name}_aln.sorted.fa'):
            fasta_file = os.path.join(input_dir, file_name)
            sample_name = file_name.replace(f'_{database_name}_aln.sorted.fa', '')
            sample_names.append(sample_name)

            n_count = 0
            total_base_count = 0
            try:
                with open(fasta_file, 'r') as file:
                    for line in file:
                        if not line.startswith('>'):
                            sequence = line.strip()
                            n_count += sequence.count('N')
                            total_base_count += len(sequence)
                n_counts.append(n_count)
                total_bases.append(total_base_count)
            except Exception as e:
                logging.error(f"Error processing FASTA file {fasta_file}: {e}")
                n_counts.append(0)
                total_bases.append(0)

    fasta_data = {
        'Sample': sample_names,
        'Number of Ns': n_counts,
        'Total Bases': total_bases
    }
    fasta_df = pd.DataFrame(fasta_data)
    fasta_df.to_excel(output_file, index=False)
    logging.info(f"FASTA summary saved to {output_file}")

def merge_excel_files(coverage_file, fasta_file, output_file):
    try:
        coverage_df = pd.read_excel(coverage_file)
        fasta_df = pd.read_excel(fasta_file)
        merged_df = pd.merge(coverage_df, fasta_df, on='Sample', how='outer')
        merged_df.fillna(0, inplace=True)
        merged_df.to_excel(output_file, index=False)
        logging.info(f"Merged summary saved to {output_file}")
    except Exception as e:
        logging.error(f"Error merging Excel files: {e}")
        sys.exit(1)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, help='Path to input directory containing coverage and FASTA files.')
    parser.add_argument('--output_dir', type=str, help='Path to output directory for summarized information.')
    parser.add_argument('--database_name', type=str, default='denv1', help='Name of the SnpEff database (default: denv1).')
    args = parser.parse_args(argv)

    input_dir = args.input_dir
    output_dir = args.output_dir
    database_name = args.database_name

    if not os.path.exists(input_dir):
        logging.error(f"Input directory {input_dir} does not exist")
        sys.exit(1)

    coverage_output_file = os.path.join(output_dir, 'coverage_summary.xlsx')
    summarize_coverage(input_dir, coverage_output_file)

    fasta_output_file = os.path.join(output_dir, 'fasta_summary.xlsx')
    summarize_fasta(input_dir, fasta_output_file, database_name)

    merged_output_file = os.path.join(output_dir, 'merged_summary.xlsx')
    merge_excel_files(coverage_output_file, fasta_output_file, merged_output_file)

if __name__ == '__main__':
    main()
