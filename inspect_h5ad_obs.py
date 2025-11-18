
import anndata as ad
import sys

def inspect_h5ad(file_path):
    """
    Loads an h5ad file and prints its obs columns and the first 5 rows.
    """
    try:
        adata = ad.read_h5ad(file_path)
        print("Columns in adata.obs:")
        print(adata.obs.columns)
        print("\nFirst 5 rows of adata.obs:")
        print(adata.obs.head())
    except Exception as e:
        print(f"Error reading {file_path}: {e}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        file_path = sys.argv[1]
        inspect_h5ad(file_path)
    else:
        print("Usage: python inspect_h5ad_obs.py <path_to_h5ad_file>")
        # As a fallback, try to inspect the known file
        known_file = "seacells_pipeline_output_17_markers/integrated_metacells_17_markers.h5ad"
        print(f"\nAttempting to inspect default file: {known_file}")
        inspect_h5ad(known_file)
