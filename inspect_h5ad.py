import scanpy as sc
import sys
from pathlib import Path

# Take the first h5ad file path from command line arguments
if len(sys.argv) < 2:
    print("Usage: python inspect_h5ad.py <path_to_h5ad_file>")
    sys.exit(1)

h5ad_path = Path(sys.argv[1])

if not h5ad_path.exists():
    print(f"Error: File not found at {h5ad_path}")
    sys.exit(1)

try:
    print(f"Inspecting {h5ad_path}...")
    adata = sc.read_h5ad(h5ad_path)

    print("\n--- adata.obs columns (assignments) ---")
    print(adata.obs.columns.tolist())

    print("\n--- adata.uns keys (config and metrics) ---")
    print(list(adata.uns.keys()))

    # Check for specific keys I expect to use
    print("\n--- Checking for specific data ---")
    if 'archetype_id_hard' in adata.obs.columns and 'archetype_entropy' in adata.obs.columns:
        print("✓ Found 'archetype_id_hard' and 'archetype_entropy' in adata.obs")
    else:
        print("✗ Did not find expected assignment/entropy columns in adata.obs")

    if 'training_metrics' in adata.uns:
        print("✓ Found 'training_metrics' in adata.uns")
        # print(adata.uns['training_metrics']) # This might be too verbose
    else:
        print("✗ Did not find 'training_metrics' in adata.uns")

    if 'config' in adata.uns:
        print("✓ Found 'config' in adata.uns")
        # print(adata.uns['config'])
    else:
        print("✗ Did not find 'config' in adata.uns")

except Exception as e:
    print(f"An error occurred: {e}")
