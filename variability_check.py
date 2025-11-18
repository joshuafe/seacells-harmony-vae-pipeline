import scanpy as sc
import pandas as pd

try:
    adata = sc.read_h5ad('seacells_output/three_groups_phenograph/integrated_metacells_harmony_phenograph.h5ad')
    
    # Check if any samples were run multiple times
    sample_summary = adata.obs.groupby('sample_id').size()
    print("Metacells per sample_id:")
    print(sample_summary)

    # Check for confounding
    confound_check = adata.obs.groupby('sample_id')['sample_type'].unique()
    print("\nUnique sample_types per sample_id:")
    print(confound_check)

    is_confounded = all(len(types) == 1 for types in confound_check)
    if is_confounded:
        print("\nResult: Confounding is confirmed. Each sample_id maps to exactly one sample_type.")
    else:
        print("\nResult: No perfect confounding found. Some sample_ids map to multiple sample_types.")

except FileNotFoundError:
    print("Error: The file 'seacells_output/three_groups_phenograph/integrated_metacells_harmony_phenograph.h5ad' was not found.")
except Exception as e:
    print(f"An error occurred: {e}")
