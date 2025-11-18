import pandas as pd
import anndata as ad
import os
from pathlib import Path

def create_sample_id_mapping():
    """
    Creates a mapping from filename to the generated sample_id from the pipeline.
    """
    sample_dirs = {
        "Normal": Path("data/normalsgated"),
        "Abnormal": Path("data/deidentify_abngated"),
        "PTCy": Path("data/ptcy_fcs/PTCy Lymphocytes")
    }

    mapping = {}
    for sample_type, directory in sample_dirs.items():
        if directory.exists():
            fcs_files = sorted([f.name for f in directory.glob("*.fcs")])
            for i, fcs_filename in enumerate(fcs_files):
                generated_id = f"{sample_type}_{i+1}"
                mapping[fcs_filename] = generated_id
    return mapping

def get_sample_files(filename_to_id_map):
    """
    Reads the sample_annotations.csv file and returns lists of sample_ids for
    the Normal, d28 PTCy, and d28 Rux cohorts.
    """
    annotations_df = pd.read_csv("sample_annotations.csv")
    annotations_df.columns = annotations_df.columns.str.strip()

    # Normal samples - directly use the 'sample_id' column
    normal_sample_ids = annotations_df[annotations_df["Notes"] == "Normal"]["sample_id"].dropna().tolist()

    # d28 PTCy samples
    ptcy_filenames = annotations_df[
        (annotations_df["tissue_type"] == "BM") &
        (annotations_df["gvh_prophylaxis"] == "PTCy") &
        (annotations_df["day_relative_to_transplant"] >= 28) &
        (annotations_df["day_relative_to_transplant"] <= 50)
    ]["filename"].dropna().tolist()
    ptcy_sample_ids = [filename_to_id_map[f] for f in ptcy_filenames if f in filename_to_id_map]

    # d28 Rux samples
    rux_filenames = annotations_df[
        (annotations_df["tissue_type"] == "BM") &
        (annotations_df["gvh_prophylaxis"] == "Rux") &
        (annotations_df["day_relative_to_transplant"] >= 28) &
        (annotations_df["day_relative_to_transplant"] <= 50)
    ]["filename"].dropna().tolist()
    rux_sample_ids = [filename_to_id_map[f] for f in rux_filenames if f in filename_to_id_map]

    return normal_sample_ids, ptcy_sample_ids, rux_sample_ids

def main():
    """
    Main function to prepare the data for the d28 comparative analysis.
    """
    # Create the mapping from filename to sample_id
    filename_to_id_map = create_sample_id_mapping()

    # Get the sample filenames
    normal_samples, ptcy_samples, rux_samples = get_sample_files(filename_to_id_map)
    
    print("Selected samples:")
    print(f"Normal: {normal_samples}")
    print(f"PTCy: {ptcy_samples}")
    print(f"Rux: {rux_samples}")

    # Define the input and output file paths
    input_h5ad_path = "seacells_pipeline_output_17_markers/integrated_metacells_17_markers.h5ad"
    output_h5ad_path = "d28_comparative_analysis.h5ad"

    # Check if the input file exists
    if not os.path.exists(input_h5ad_path):
        print(f"Error: Input file not found at {input_h5ad_path}")
        print("Please run the data regeneration pipeline first (scripts/pipeline_step1_seacells.py).")
        return

    # Load the integrated metacells data
    print(f"Loading data from {input_h5ad_path}...")
    adata = ad.read_h5ad(input_h5ad_path)

    # Combine the sample lists
    selected_samples = normal_samples + ptcy_samples + rux_samples

    # Filter the anndata object
    print(f"Filtering for {len(selected_samples)} selected samples...")
    adata_filtered = adata[adata.obs["sample_id"].isin(selected_samples)].copy()

    # Save the filtered data
    print(f"Saving filtered data to {output_h5ad_path}...")
    adata_filtered.write_h5ad(output_h5ad_path)

    print("Data preparation complete.")

if __name__ == "__main__":
    main()