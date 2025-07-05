import argparse
import os
import subprocess

# Usage: python create_snpeff_database.py --genbank_file denv1.gb --reference_fasta references/denv1.fasta --output_dir test_out --database_name denv1

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

def validate_genbank_file(genbank_file):
    try:
        with open(genbank_file, "r") as f:
            content = f.read()
            if not content.startswith("LOCUS"):
                raise ValueError("GenBank file does not start with 'LOCUS'. Invalid format.")
            if "FEATURES" not in content:
                raise ValueError("GenBank file lacks 'FEATURES' section. Annotations may be missing.")
        print(f"GenBank file {genbank_file} appears valid.")
    except Exception as e:
        raise Exception(f"GenBank file validation failed: {str(e)}")

def validate_fasta_genbank_match(genbank_file, reference_fasta):
    try:
        with open(reference_fasta, "r") as f:
            fasta_id = f.readline().strip().lstrip(">").split()[0]
        
        with open(genbank_file, "r") as f:
            for line in f:
                if line.startswith("LOCUS"):
                    gbk_id = line.split()[1]
                    break
            else:
                raise ValueError("No LOCUS line found in GenBank file.")
        
        if fasta_id != gbk_id:
            raise ValueError(f"Sequence ID mismatch: FASTA ID '{fasta_id}' does not match GenBank ID '{gbk_id}'.")
        print(f"Sequence IDs match: {fasta_id}")
    except Exception as e:
        raise Exception(f"FASTA-GenBank validation failed: {str(e)}")

def create_snpeff_config(output_dir, database_name, reference_fasta):
    data_dir = os.path.abspath(os.path.join(output_dir, "data"))
    custom_config = os.path.join(output_dir, "snpEff.config")
    config_content = f"""# SnpEff configuration file
data.dir = {data_dir}

# Database entry
{database_name}.genome = {database_name}
{database_name}.reference = {os.path.abspath(reference_fasta)}
"""
    with open(custom_config, "w") as config_file:
        config_file.write(config_content)
    
    print(f"SnpEff configuration file created: {custom_config}")
    return custom_config

def prepare_files(genbank_file, reference_fasta, output_dir, database_name):
    data_dir = os.path.join(output_dir, "data", database_name)
    os.makedirs(data_dir, exist_ok=True)
    fasta_dest = os.path.join(data_dir, "sequences.fa")
    genbank_dest = os.path.join(data_dir, "genes.gbk")
    run_command(f"cp {reference_fasta} {fasta_dest}")
    run_command(f"cp {genbank_file} {genbank_dest}")
    # Verify files exist
    if not os.path.exists(fasta_dest):
        raise FileNotFoundError(f"FASTA file not found at {fasta_dest}")
    if not os.path.exists(genbank_dest):
        raise FileNotFoundError(f"GenBank file not found at {genbank_dest}")
    print(f"Files prepared: {fasta_dest}, {genbank_dest}")

def build_snpeff_database(database_name, config_file, output_dir):
    data_dir = os.path.abspath(os.path.join(output_dir, "data"))
    build_command = (
        f"snpEff build -genbank -v {database_name} "
        f"-c {config_file} "
        f"-dataDir {data_dir}"
    )
    stdout, stderr = run_command(build_command)
    print("SnpEff build output:", stdout)
    print("SnpEff build error:", stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--genbank_file", type=str, help="Path to GenBank (.gb or .gbk) file.")
    parser.add_argument("--reference_fasta", type=str, help="Path to reference FASTA file.")
    parser.add_argument("--output_dir", type=str, help="Path to output directory for SnpEff database.")
    parser.add_argument("--database_name", type=str, default="denv1", help="Name of the SnpEff database (default: denv1).")
    parser.add_argument("--skip_id_validation", action="store_true", help="Skip validation of sequence ID matching.")
    args = parser.parse_args()

    genbank_file = args.genbank_file
    reference_fasta = args.reference_fasta
    output_dir = args.output_dir
    database_name = args.database_name

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    try:
        validate_genbank_file(genbank_file)
        if not args.skip_id_validation:
            validate_fasta_genbank_match(genbank_file, reference_fasta)
        config_file = create_snpeff_config(output_dir, database_name, reference_fasta)
        prepare_files(genbank_file, reference_fasta, output_dir, database_name)
        build_snpeff_database(database_name, config_file, output_dir)
        print(f"SnpEff database '{database_name}' created successfully in {output_dir}")
    except Exception as e:
        print("Error occurred during database creation:", str(e))
