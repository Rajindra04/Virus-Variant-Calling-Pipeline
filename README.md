# Dengue Variant Calling Pipeline

This repository contains a bioinformatics pipeline for processing paired-end FASTQ files to perform read mapping, variant calling, consensus sequence generation, and annotation for dengue virus sequences. The pipeline is automated through a master script, `run_pipeline.py`, which orchestrates the execution of individual scripts in the correct order.

## Table of Contents
- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
- [Directory Structure](#directory-structure)
- [Input Files](#input-files)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)

## Overview
The pipeline processes paired-end FASTQ files to:
1. Generate a sample sheet from FASTQ files.
2. Map reads to a reference genome using `bwa-mem2`.
3. Convert SAM to BAM files and sort/index them.
4. Create a SnpEff database for annotation.
5. Generate consensus sequences and initial VCFs using `ivar` and `bcftools`.
6. Perform variant calling and annotation with `GATK` and `SnpEff`.
7. Summarize coverage and FASTA statistics.
8. Summarize SnpEff annotations for variant impact analysis.

The master script, `run_pipeline.py`, automates these steps, ensuring proper input/output handling and error checking.

## Prerequisites
- **Python 3.6+** with the following packages:
  - `pandas`
  - `matplotlib`
  - `argparse`
  - `subprocess`
  - `os`
  - `glob`
- **Bioinformatics Tools** (ensure they are in your system PATH):
  - `bwa-mem2`
  - `samtools`
  - `fastp`
  - `fastqc`
  - `gatk`
  - `snpEff`
  - `ivar`
  - `bcftools`
- **Input Files**:
  - Paired-end FASTQ files (e.g., `sample_R1.fastq.gz`, `sample_R2.fastq.gz`).
  - Reference FASTA file (e.g., `denv1.fasta`).
  - GenBank file for SnpEff database (e.g., `denv1.gb`).

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/dengue-pipeline.git
   cd dengue-pipeline
   ```
2. Install Python dependencies:
   ```bash
   pip install pandas matplotlib
   ```
3. Ensure all bioinformatics tools are installed and accessible. For example:
   - Install `samtools`: `sudo apt-get install samtools` (Ubuntu) or equivalent.
   - Install `snpEff`: Follow instructions at [SnpEff documentation](http://pcingola.github.io/SnpEff/).
   - Verify tools are in PATH: `bwa-mem2 version`, `samtools --version`, etc.

## Usage
Run the pipeline with a single command using `run_pipeline.py`:

```bash
python run_pipeline.py --input_dir /path/to/fastq_files --reference_fasta /path/to/denv1.fasta --genbank_file /path/to/denv1.gb --output_dir /path/to/output --database_name denv1
```

### Command-Line Arguments
- `--input_dir`: Directory containing paired-end FASTQ files.
- `--reference_fasta`: Path to the reference FASTA file.
- `--genbank_file`: Path to the GenBank file for SnpEff annotation.
- `--output_dir`: Directory to store all output files.
- `--database_name`: Name of the SnpEff database (default: `denv1`).

### Example
```bash
python run_pipeline.py --input_dir ./fastq_data --reference_fasta ./references/denv1.fasta --genbank_file ./references/denv1.gb --output_dir ./output --database_name denv1
```

## Pipeline Steps
The pipeline executes the following scripts in order:
1. **`1.create_samplesheet.py`**: Generates a sample sheet (`samplesheet.tsv`) from FASTQ files.
2. **`2.map_reads.py`**: Maps reads to the reference using `bwa-mem2`, trims with `fastp`, and runs `FastQC`.
3. **`3.samtobamdenv.py`**: Converts SAM to BAM, sorts, and indexes BAM files.
4. **`5.create_snpeff_database.py`**: Creates a SnpEff database from the GenBank and FASTA files.
5. **`4.sam2consensus_test2_ivar.py`**: Generates consensus sequences and initial VCFs using `ivar` and `bcftools`.
6. **`6.variant_calling_consensus.py`**: Performs variant calling with `GATK` and annotates with `SnpEff`.
7. **`7.summarize_result.py`**: Summarizes coverage and FASTA statistics in Excel files.
8. **`8.summarize_snpEff.py`**: Summarizes SnpEff annotations into a CSV and JSON for visualization.

## Directory Structure
```plaintext
dengue-pipeline/
├── 1.create_samplesheet.py
├── 2.map_reads.py
├── 3.samtobamdenv.py
├── 4.sam2consensus_test2_ivar.py
├── 5.create_snpeff_database.py
├── 6.variant_calling_consensus.py
├── 7.summarize_result.py
├── 8.summarize_snpEff.py
├── run_pipeline.py
├── fastq_data/               # Input directory for FASTQ files
├── references/               # Directory for reference FASTA and GenBank files
└── output/                   # Output directory for results
```

## Input Files
- **FASTQ Files**: Paired-end files named with `_R1_` and `_R2_` (e.g., `sample1_R1.fastq.gz`, `sample1_R2.fastq.gz`).
- **Reference FASTA**: A single FASTA file (e.g., `denv1.fasta`) containing the reference genome.
- **GenBank File**: A GenBank file (e.g., `denv1.gb`) for SnpEff annotation.

## Output Files
The pipeline generates the following outputs in the `output_dir`:
- `samplesheet.tsv`: Sample sheet with sample names and FASTQ paths.
- `sam_files/`: Directory containing SAM files from read mapping.
- `*_output/`: Sample-specific directories with trimmed FASTQ and FastQC reports.
- `*.sorted.bam`: Sorted and indexed BAM files.
- `*.fa`: Consensus FASTA sequences.
- `*.vcf`: Variant call files (both raw and annotated).
- `*_coverage.txt`: Coverage data files.
- `*_coverage.png`: Coverage plots.
- `*_snpEff_summary.html`, `*_snpEff_genes.txt`: SnpEff annotation reports.
- `coverage_summary.xlsx`, `fasta_summary.xlsx`, `merged_summary.xlsx`: Summary statistics.
- `summary_table.csv`, `chart_data.json`: SnpEff annotation summaries.

## Troubleshooting
- **File Not Found Errors**: Ensure all input files and directories exist. Check paths provided to `--input_dir`, `--reference_fasta`, and `--genbank_file`.
- **Tool Not Found**: Verify that all bioinformatics tools are installed and in your system PATH (`which samtools`, `which bwa-mem2`, etc.).
- **SnpEff Directory Issue**: The `8.summarize_snpEff.py` script uses a hardcoded input directory. Update the `input_dir` variable in the script to match your `output_dir` or modify the script to accept it as an argument.
- **Permission Issues**: Ensure you have write permissions for the output directory.
- **Column Mismatch in SnpEff Summaries**: If `8.summarize_snpEff.py` reports column mismatches, verify that SnpEff output files (`*_snpEff_summary.genes.txt`) are correctly generated and contain expected columns.

For detailed error logs, check the console output or log files generated by `run_pipeline.py`.

## Contributing
Contributions are welcome! Please submit a pull request or open an issue for bug reports, feature requests, or improvements.

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
