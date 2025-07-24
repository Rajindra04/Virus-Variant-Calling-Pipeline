# Virus Variant Calling Pipeline

This repository provides a complete bioinformatics pipeline for processing paired-end FASTQ files to perform read mapping, variant calling, consensus sequence generation, and functional annotation of dengue virus genomes. The pipeline is orchestrated via the `run_pipeline` command, which automates the workflow using modular Python scripts.

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

The pipeline performs the following key steps:

- Generates a sample sheet from FASTQ files.
- Maps reads to a reference genome using `bwa-mem2`.
- Converts SAM to sorted/indexed BAM files.
- Creates a SnpEff database using GenBank annotations.
- Generates consensus sequences and raw VCFs using `ivar` and `bcftools`.
- Performs variant calling with GATK, annotates variants using SnpEff, and extracts fields using SnpSift.
- Summarizes alignment coverage and consensus FASTA statistics.
- Summarizes functional annotations (impacts and gene effects).

The pipeline is executed by the `run_pipeline` entry point, which ensures proper sequencing, I/O handling, and error checking.

## Prerequisites

- Python 3.6+
  - Packages: `pandas`, `matplotlib`, `argparse`
- Bioinformatics Tools (must be in PATH):
  - `bwa-mem2`, `samtools`, `fastp`, `fastqc`, `gatk`, `snpEff`, `SnpSift`, `ivar`, `bcftools`

## Installation

### Option 1: From GitHub

```bash
git clone https://github.com/Rajindra04/Virus-Variant-Calling-Pipeline.git
cd Virus-Variant-Calling-Pipeline
conda create -n dengue_pipeline python=3.10
conda activate dengue_pipeline
pip install -r requirements.txt
```

### Option 2: From Conda (once published)

```bash
conda install -c your-channel virus-variant-calling-pipeline
```

### Bioinformatics Tools via Conda

```bash
conda install -c bioconda bwa-mem2 samtools fastp fastqc gatk4 snpeff snpsift ivar bcftools
```

Verify installations:

```bash
which bwa-mem2
samtools --version
SnpSift
```

## Usage

```bash
run_pipeline \
  --input_dir ./fastq_data \
  --reference_fasta ./references/denv1.fasta \
  --genbank_file ./references/denv1.gb \
  --output_dir ./output \
  --database_name denv1
```

### Command-Line Arguments

| Argument            | Description                                              |
|--------------------|----------------------------------------------------------|
| `--input_dir`      | Directory containing paired-end FASTQ files              |
| `--reference_fasta`| Path to the reference FASTA file                         |
| `--genbank_file`   | Path to the GenBank file for SnpEff annotation           |
| `--output_dir`     | Directory where output will be saved                     |
| `--database_name`  | Name of the custom SnpEff database (default: `denv1`)    |

## Pipeline Steps

| Step | Script                          | Description                                                       |
|------|---------------------------------|-------------------------------------------------------------------|
| 1    | `1.create_samplesheet.py`       | Create `samplesheet.tsv` from FASTQ files                         |
| 2    | `2.map_reads.py`                | Read mapping with `bwa-mem2`, trimming, and quality control       |
| 3    | `3.samtobamdenv.py`             | Convert SAM to sorted/indexed BAM                                 |
| 4    | `5.create_snpeff_database.py`   | Build custom SnpEff database                                     |
| 5    | `4.sam2consensus_test2_ivar.py` | Generate consensus FASTA and raw VCFs                            |
| 6    | `6.variant_calling_consensus.py`| Call variants, annotate with SnpEff + SnpSift                    |
| 7    | `7.summarize_result.py`         | Summarize coverage, depth, and consensus FASTAs                  |
| 8    | `8.summarize_snpEff.py`         | Summarize annotation impacts and gene-level effects              |

## Directory Structure

```
Virus-Variant-Calling-Pipeline/
├── pipeline/                  # Python package directory
│   ├── __init__.py
│   ├── run_pipeline.py
│   ├── 1.create_samplesheet.py
│   ├── 2.map_reads.py
│   ├── 3.samtobamdenv.py
│   ├── 4.sam2consensus_test2_ivar.py
│   ├── 5.create_snpeff_database.py
│   ├── 6.variant_calling_consensus.py
│   ├── 7.summarize_result.py
│   └── 8.summarize_snpEff.py
├── fastq_data/
├── references/
└── output/
```

## Input Files

- **FASTQ Files**: Paired-end reads (`*_R1.fastq.gz`, `*_R2.fastq.gz`)
- **Reference Genome**: FASTA file (`denv1.fasta`)
- **GenBank File**: Annotation file for SnpEff (`denv1.gb`)

## Output Files

- `samplesheet.tsv`: List of samples and file paths
- `*.sorted.bam`: Sorted and indexed BAM files
- `*.fa`: Consensus FASTA sequences
- `*.vcf`: Raw and annotated variant call files
- `*_snpEff_summary.*`: HTML, CSV, and TXT reports from SnpEff
- `*_snpSift.txt`: Extracted fields from annotated VCFs
- `coverage_summary.xlsx`: Read depth summaries
- `fasta_summary.xlsx`: Summary of consensus sequences
- `summary_table.csv`, `.json`: Combined annotation summaries

## Troubleshooting

### File Errors

- Make sure all input files exist and correct paths are provided.
- Check that the output directory has write permissions:
  ```bash
  chmod -R u+rw output
  ```

### Tool Issues

- Check tool availability:
  ```bash
  which samtools
  which SnpSift
  ```

### SnpEff Database Problems

- Missing or incomplete `snpEffectPredictor.bin` or `genes.gbk`:
  Rebuild using:
  ```bash
  python pipeline/5.create_snpeff_database.py \
    --genbank_file references/denv1.gb \
    --reference_fasta references/denv1.fasta \
    --output_dir output \
    --database_name denv1
  ```

- Replace broken GenBank:
  ```bash
  wget -O references/denv1.gb \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NC_001477.1&rettype=gb&retmode=text"
  ```

## Contributing

Contributions, issues, and feature requests are welcome. Please submit a pull request or open an issue on GitHub.

## License

This project is licensed under the MIT License.
