import os
import argparse
import pandas as pd

def summarize_coverage(input_dir, output_file):
    # Initialize data lists
    sample_names = []
    average_coverage = []

    # Iterate through the files in the input directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith('_coverage.txt'):
            # Process coverage file
            coverage_file = os.path.join(input_dir, file_name)
            sample_name = file_name.split('_')[0]
            sample_names.append(sample_name)

            # Calculate average coverage
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
            except FileNotFoundError:
                print(f"Coverage file not found for sample: {sample_name}")
                average_coverage.append(0)

    # Create a DataFrame with the summarized coverage information
    coverage_data = {
        'Sample': sample_names,
        'Average Coverage': average_coverage,
    }
    coverage_df = pd.DataFrame(coverage_data)

    # Save the coverage DataFrame to an Excel file
    coverage_df.to_excel(output_file, index=False)

def summarize_fasta(input_dir, output_file):
    # Initialize data lists
    sample_names = []
    n_counts = []
    total_bases = []

    # Iterate through the files in the input directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith('.fa'):
            # Process FASTA file
            fasta_file = os.path.join(input_dir, file_name)
            sample_name = file_name.replace('_denv1_aln.sorted.fa', '')
            sample_names.append(sample_name)

            # Calculate number of 'N' bases and total bases
            n_count = 0
            total_base_count = 0
            with open(fasta_file, 'r') as file:
                for line in file:
                    if not line.startswith('>'):
                        sequence = line.strip()
                        n_count += sequence.count('N')
                        total_base_count += len(sequence)
            n_counts.append(n_count)
            total_bases.append(total_base_count)

    # Create a DataFrame with the summarized FASTA information
    fasta_data = {
        'Sample': sample_names,
        'Number of Ns': n_counts,
        'Total Bases': total_bases
    }
    fasta_df = pd.DataFrame(fasta_data)

    # Save the FASTA DataFrame to an Excel file
    fasta_df.to_excel(output_file, index=False)

def merge_excel_files(coverage_file, fasta_file, output_file):
    # Read the coverage and FASTA Excel files
    coverage_df = pd.read_excel(coverage_file)
    fasta_df = pd.read_excel(fasta_file)

    # Merge the two DataFrames based on the 'Sample' column
    merged_df = pd.merge(coverage_df, fasta_df, on='Sample', how='outer')

    # Fill NaN values in the merged DataFrame
    merged_df.fillna(0, inplace=True)

    # Save the merged DataFrame to an Excel file
    merged_df.to_excel(output_file, index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, help='Path to input directory containing coverage and FASTA files.')
    parser.add_argument('--output_dir', type=str, help='Path to output directory for summarized information.')
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir

    # Summarize coverage
    coverage_output_file = os.path.join(output_dir, 'coverage_summary.xlsx')
    summarize_coverage(input_dir, coverage_output_file)

    # Summarize FASTA
    fasta_output_file = os.path.join(output_dir, 'fasta_summary.xlsx')
    summarize_fasta(input_dir, fasta_output_file)

    # Merge the two summary Excel files
    merged_output_file = os.path.join(output_dir, 'merged_summary.xlsx')
    merge_excel_files(coverage_output_file, fasta_output_file, merged_output_file)
