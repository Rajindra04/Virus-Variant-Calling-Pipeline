Virus Variant Calling Pipeline

This repository contains a bioinformatics pipeline for processing paired-end FASTQ files to perform read mapping, variant calling, consensus sequence generation, and annotation for dengue virus sequences. The pipeline is automated through a master script, run_pipeline.py, or callable via a Conda-installed CLI command, dengue_pipeline.
Table of Contents

    Overview

    Prerequisites

    Installation

    Usage

    Pipeline Steps

    Directory Structure

    Input Files

    Output Files

    Troubleshooting

    Contributing

    License

Overview

The pipeline processes paired-end FASTQ files to:

    Generate a sample sheet from FASTQ files.

    Map reads to a reference genome using bwa-mem2.

    Convert SAM to BAM, sort, and index the alignments.

    Generate consensus sequences using ivar and bcftools.

    Create a SnpEff database for annotation.

    Perform variant calling using GATK, annotate variants with SnpEff, and extract annotation fields using SnpSift.

    Summarize coverage and FASTA stats.

    Summarize SnpEff annotations for variant impact analysis.

The run_pipeline.py script automates these steps with robust error handling and path management.
Prerequisites
Python

    Python ≥ 3.8

    Required packages:

        pandas

        matplotlib

        pyyaml

        argparse

Bioinformatics Tools

(Install via Bioconda)

    bwa-mem2

    samtools

    fastp

    fastqc

    gatk4

    snpeff

    snpsift

    ivar

    bcftools

Input Files

    Paired-end FASTQ files (e.g., sample1_R1.fastq.gz, sample1_R2.fastq.gz)

    Reference FASTA (e.g., denv1.fasta)

    GenBank file for building SnpEff database (e.g., denv1.gb)

Installation
Option 1: Manual

git clone https://github.com/Rajindra04/Virus-Variant-Calling-Pipeline.git
cd Virus-Variant-Calling-Pipeline
pip install -r requirements.txt

Ensure all external tools (listed above) are installed and in your PATH.
Option 2: Using Conda (recommended)

Once released via Anaconda or Bioconda:

conda install -c bioconda dengue-pipeline

Then, the pipeline can be run with:

dengue_pipeline --input_dir ...  # (see below)

Usage
From source:

python run_pipeline.py \
  --input_dir ./fastq_data \
  --reference_fasta ./references/denv1.fasta \
  --genbank_file ./references/denv1.gb \
  --output_dir ./output \
  --database_name denv1

From Conda-installed CLI (after packaging):

dengue_pipeline \
  --input_dir ./fastq_data \
  --reference_fasta ./references/denv1.fasta \
  --genbank_file ./references/denv1.gb \
  --output_dir ./output \
  --database_name denv1

Pipeline Steps
Script	Description
1.create_samplesheet.py	Auto-generates a sample sheet from FASTQ filenames
2.map_reads.py	Maps reads using bwa-mem2, trims with fastp, and runs FastQC
3.samtobamdenv.py	Converts SAM → BAM, sorts, indexes
4.sam2consensus_test2_ivar.py	Generates consensus FASTAs using ivar and bcftools
5.create_snpeff_database.py	Builds SnpEff database from GenBank and FASTA
6.variant_calling_consensus.py	Calls variants with GATK, annotates with SnpEff, extracts with SnpSift
7.summarize_result.py	Summarizes coverage and FASTA statistics
8.summarize_snpEff.py	Summarizes SnpEff annotations for variant effect interpretation

Directory Structure

Virus-Variant-Calling-Pipeline/
├── virus_pipeline/
│   ├── __init__.py
│   ├── pipeline/
│   │   ├── 1.create_samplesheet.py
│   │   ├── 2.map_reads.py
│   │   ├── 3.samtobamdenv.py
│   │   ├── 4.sam2consensus_test2_ivar.py
│   │   ├── 5.create_snpeff_database.py
│   │   ├── 6.variant_calling_consensus.py
│   │   ├── 7.summarize_result.py
│   │   └── 8.summarize_snpEff.py
│   └── run_pipeline.py
├── output/                   # Output folder
├── fastq_data/              # Input FASTQ files
├── references/              # Reference FASTA and GenBank
├── LICENSE
├── README.md
├── environment.yml
└── setup.py

Input Files
File Type	Description
FASTQ	Paired-end reads (e.g., sample1_R1.fastq.gz)
FASTA	Reference genome (e.g., denv1.fasta)
GenBank	Genomic annotations for SnpEff (e.g., denv1.gb)

Output Files
File	Description
samplesheet.tsv	Sample names and input file paths
*_output/	FastQC reports and trimmed FASTQ
*.sorted.bam, .bai	Sorted BAM alignments
*.fa, *.vcf	Consensus FASTA and variant calls
*_coverage.txt, *.png	Coverage summaries and plots
*_snpEff_summary.csv	Annotated variant summaries
summary_table.csv, chart_data.json	JSON/CSV summaries for plotting
Troubleshooting

    ❗ File not found
    Double-check --input_dir, --reference_fasta, and --genbank_file.

    ❗ Missing tools
    Ensure tools like samtools, snpEff, SnpSift are in your PATH:

which samtools
which SnpSift

❗ SnpEff DB not found
Make sure create_snpeff_database.py is run before annotation steps.

❗ Permission errors
Ensure the output directory is writable:

    chmod -R u+rw output/

    ❗ Unexpected VCF columns
    Inspect annotated VCFs to ensure proper ANN fields are present.

Contributing

Contributions are welcome!
Feel free to submit a pull request or file an issue for bugs, improvements, or feature ideas.
License

This project is licensed under the MIT License.
