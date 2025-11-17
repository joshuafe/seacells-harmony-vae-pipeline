"""
Generates a CSV template for sample annotation.

This script scans the data directories to find all original samples,
checks their processing status from the pipeline output, and creates
a CSV file with placeholder columns for clinical metadata.
"""
import pandas as pd
from pathlib import Path
import numpy as np

def generate_annotation_template():
    """
    Generates and saves a sample annotation template CSV.
    """
    print("Generating sample annotation template...")

    # --- Step 1: Find all original samples ---
    sample_dirs = {
        "Normal": Path("data/normalsgated"),
        "Abnormal": Path("data/deidentify_abngated"),
        "PTCy": Path("data/ptcy_fcs/PTCy Lymphocytes")
    }
    
    all_sample_ids = []
    for sample_type, directory in sample_dirs.items():
        if directory.exists():
            for i, fcs_path in enumerate(sorted(directory.glob("*.fcs"))):
                all_sample_ids.append(f"{sample_type}_{i+1}")
    
    if not all_sample_ids:
        print("Error: No FCS files found. Cannot generate template.")
        return

    print(f"Found {len(all_sample_ids)} total original samples.")

    # --- Step 2: Check processing status ---
    output_dir = Path("seacells_pipeline_output_17_markers")
    processed_samples = set()
    if output_dir.exists():
        for f in output_dir.glob("*_metacells.h5ad"):
            # Extract sample_id from filename like 'Normal_1_metacells.h5ad'
            sample_id = '_'.join(f.name.split('_')[:-1])
            processed_samples.add(sample_id)
    
    print(f"Found {len(processed_samples)} successfully processed samples.")

    # --- Step 3: Create the DataFrame ---
    df = pd.DataFrame(index=all_sample_ids)
    df.index.name = 'sample_id'
    
    # Add user-requested columns
    df['biobank_id'] = ''
    df['tissue_type'] = ''
    df['day_relative_to_transplant'] = ''
    df['conditioning'] = ''
    df['gvh_prophylaxis'] = ''
    
    # Add processing status column
    df['processing_status'] = 'Skipped - Missing CD279+CD24'
    
    # Update status for processed samples
    # Use .loc for safe assignment
    processed_mask = df.index.isin(processed_samples)
    df.loc[processed_mask, 'processing_status'] = 'Processed'

    # --- Step 4: Pre-fill 'NA' for normal samples ---
    normal_mask = df.index.str.contains('Normal')
    transplant_cols = ['day_relative_to_transplant', 'conditioning', 'gvh_prophylaxis']
    df.loc[normal_mask, transplant_cols] = 'NA'

    # --- Step 5: Save the template ---
    output_path = "sample_annotation_template.csv"
    df.to_csv(output_path)
    
    print(f"\nSuccessfully created annotation template: {output_path}")
    print("Please fill in the required metadata in this file for downstream analysis.")

if __name__ == "__main__":
    generate_annotation_template()
