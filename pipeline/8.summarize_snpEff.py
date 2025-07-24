import sys

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
import pandas as pd
import glob
import os
import json

# Directory containing the input files
input_dir = "/home/rajindra/Documents/dengue/denv3_out"

# File pattern
file_pattern = "*_aln.sorted_snpEff_summary.genes.txt"

# Expected columns to keep (E to M)
columns_to_keep = [
    "TranscriptId",
    "variants_impact_HIGH",
    "variants_impact_LOW",
    "variants_impact_MODERATE",
    "variants_impact_MODIFIER",
    "variants_effect_downstream_gene_variant",
    "variants_effect_frameshift_variant",
    "variants_effect_missense_variant",
    "variants_effect_splice_region_variant",
    "variants_effect_synonymous_variant"
]

# Full column list for files with headers (includes all possible columns)
full_columns = [
    "GeneName", "GeneId", "TranscriptId", "BioType",
    "variants_impact_HIGH", "variants_impact_LOW", "variants_impact_MODERATE", "variants_impact_MODIFIER",
    "variants_effect_conservative_inframe_deletion",  # Added for D-26_denv1
    "variants_effect_disruptive_inframe_insertion",   # Added for D1-123_denv1, etc.
    "variants_effect_downstream_gene_variant", "variants_effect_frameshift_variant",
    "variants_effect_missense_variant", "variants_effect_splice_region_variant",
    "variants_effect_stop_gained", "variants_effect_synonymous_variant", "variants_effect_upstream_gene_variant"
]

# Initialize an empty list to store dataframes
all_data = []

# Check if directory exists
if not os.path.exists(input_dir):
    print(f"Error: Directory '{input_dir}' does not exist.")
    exit(1)

# Get list of files
files = glob.glob(os.path.join(input_dir, file_pattern))
if not files:
    print(f"Error: No files found matching pattern '{file_pattern}' in '{input_dir}'.")
    print(f"Files in directory: {os.listdir(input_dir)}")
    exit(1)

# Iterate over all files matching the pattern
for file_path in files:
    # Extract sample name from filename
    sample_name = os.path.basename(file_path).replace("_aln.sorted_snpEff_summary.genes.txt", "")
    print(f"\nProcessing file: {file_path} (Sample: {sample_name})")
    
    try:
        # Read first few lines for debugging
        with open(file_path, 'r') as f:
            print(f"First 3 lines of {file_path}:")
            lines = []
            for i in range(3):
                try:
                    line = next(f).strip()
                    lines.append(line)
                    print(f"Line {i+1}: {line}")
                    fields = line.split("\t")
                    print(f"Fields in line {i+1}: {len(fields)} (Expected: 13, 14, or 15)")
                    if i > 0 and len(fields) not in [13, 14, 15]:
                        print(f"Warning: Incorrect field count in line {i+1}")
                except StopIteration:
                    break
        
        # Try reading with tab delimiter, second row as header
        try:
            df = pd.read_csv(file_path, sep="\t", comment="#", header=1)
        except Exception as e:
            print(f"Tab delimiter failed: {str(e)}. Trying comma delimiter...")
            try:
                df = pd.read_csv(file_path, sep=",", comment="#", header=1)
            except Exception as e:
                print(f"Comma delimiter failed: {str(e)}. Skipping file.")
                continue
        
        # Print available columns for debugging
        print(f"Columns in {file_path}: {list(df.columns)}")
        
        # Validate column count
        if len(df.columns) not in [13, 14, 15]:
            print(f"Error: Column count mismatch in {file_path}. Expected 13, 14, or 15, found {len(df.columns)}")
            continue
        
        # Check if header matches expected columns
        expected_cols = set(full_columns)
        actual_cols = set(df.columns)
        if not any(col in actual_cols for col in columns_to_keep):
            print(f"Error: No expected columns found in {file_path}. Attempting manual column assignment.")
            if len(df.columns) in [13, 14, 15]:
                print(f"Assigning columns: {full_columns[:len(df.columns)]}")
                df.columns = full_columns[:len(df.columns)]
            else:
                continue
        
        # Handle missing columns in columns_to_keep
        for col in columns_to_keep:
            if col not in df.columns:
                print(f"Warning: Column {col} missing in {file_path}. Filling with zeros.")
                df[col] = 0
        
        # Verify required columns
        missing_cols = [col for col in columns_to_keep if col not in df.columns]
        if missing_cols:
            print(f"Error: Missing columns in {file_path}: {missing_cols}")
            continue
        
        # Select only the requested columns
        df = df[columns_to_keep].copy()
        
        # Add sample name as a column
        df["Sample"] = sample_name
        
        # Append to the list
        all_data.append(df)
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        continue

# Check if any data was collected
if not all_data:
    print("Error: No valid dataframes to concatenate. Check file contents or column names.")
    exit(1)

# Concatenate all dataframes
combined_df = pd.concat(all_data, ignore_index=True)

# Group by TranscriptId and aggregate
summary_df = combined_df.groupby("TranscriptId").agg({
    "variants_impact_HIGH": "sum",
    "variants_impact_LOW": "sum",
    "variants_impact_MODERATE": "sum",
    "variants_impact_MODIFIER": "sum",
    "variants_effect_downstream_gene_variant": "sum",
    "variants_effect_frameshift_variant": "sum",
    "variants_effect_missense_variant": "sum",
    "variants_effect_splice_region_variant": "sum",
    "variants_effect_synonymous_variant": "sum",
    "Sample": lambda x: ";".join(x)
}).reset_index()

# Rename columns for clarity
summary_df.columns = [
    "TranscriptId",
    "Total_variants_impact_HIGH",
    "Total_variants_impact_LOW",
    "Total_variants_impact_MODERATE",
    "Total_variants_impact_MODIFIER",
    "Total_variants_effect_downstream_gene_variant",
    "Total_variants_effect_frameshift_variant",
    "Total_variants_effect_missense_variant",
    "Total_variants_effect_splice_region_variant",
    "Total_variants_effect_synonymous_variant",
    "Samples"
]

# Save the summary table to a CSV file
output_file = os.path.join(input_dir, "summary_table.csv")
summary_df.to_csv(output_file, index=False)
print(f"Summary table saved to {output_file}")

# Prepare data for Chart.js visualization (focusing on impact columns)
chart_data = {
    "labels": summary_df["TranscriptId"].tolist(),
    "datasets": [
        {
            "label": "High Impact",
            "data": summary_df["Total_variants_impact_HIGH"].tolist(),
            "backgroundColor": "rgba(255, 99, 132, 0.7)"
        },
        {
            "label": "Low Impact",
            "data": summary_df["Total_variants_impact_LOW"].tolist(),
            "backgroundColor": "rgba(54, 162, 235, 0.7)"
        },
        {
            "label": "Moderate Impact",
            "data": summary_df["Total_variants_impact_MODERATE"].tolist(),
            "backgroundColor": "rgba(255, 206, 86, 0.7)"
        },
        {
            "label": "Modifier Impact",
            "data": summary_df["Total_variants_impact_MODIFIER"].tolist(),
            "backgroundColor": "rgba(75, 192, 192, 0.7)"
        }
    ]
}

# Save chart data to JSON
chart_json_file = os.path.join(input_dir, "chart_data.json")
with open(chart_json_file, "w") as f:
    json.dump(chart_data, f, indent=2)
print(f"Chart data saved to {chart_json_file}")
