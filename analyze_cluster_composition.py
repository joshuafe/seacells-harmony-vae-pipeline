
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse

def analyze_composition(adata, cluster_key, output_dir):
    """
    Analyzes the composition of clusters by sample type (Normal, PTCy, Rux).

    Args:
        adata: AnnData object with clustering results.
        cluster_key: The key in adata.obs that contains the cluster labels.
        output_dir: Directory to save the plots.
    """
    print(f"--- Analyzing composition for {cluster_key} ---")

    # Create a DataFrame for analysis
    df = adata.obs[[cluster_key, 'sample_id']].copy()

    # Map sample_id to cohort
    def get_cohort(sample_id):
        if 'Normal' in sample_id:
            return 'Normal'
        elif 'PTCy' in sample_id:
            return 'PTCy'
        elif 'Abnormal' in sample_id: # Rux samples are in the 'Abnormal' group
            return 'Rux'
        return 'Unknown'

    df['cohort'] = df['sample_id'].apply(get_cohort)

    # Calculate the number of metacells from each cohort in each cluster
    composition = df.groupby([cluster_key, 'cohort']).size().unstack(fill_value=0)

    # Normalize to get the proportion of each cohort in each cluster
    composition_proportions = composition.div(composition.sum(axis=1), axis=0)

    # Save the composition table
    composition_proportions.to_csv(output_dir / f"{cluster_key}_composition.csv")
    print(f"  Saved composition table to {output_dir / f'{cluster_key}_composition.csv'}")

    # Plot the composition
    fig, ax = plt.subplots(figsize=(12, max(8, len(composition_proportions) * 0.4)))
    composition_proportions.plot(kind='barh', stacked=True, ax=ax, colormap='viridis')

    ax.set_xlabel("Proportion of Metacells")
    ax.set_ylabel("Cluster")
    ax.set_title(f"Cluster Composition by Cohort ({cluster_key})")
    plt.tight_layout()
    plt.savefig(output_dir / f"{cluster_key}_composition.png", dpi=300)
    print(f"  Saved composition plot to {output_dir / f'{cluster_key}_composition.png'}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Analyze cluster composition by cohort.")
    parser.add_argument(
        "--input",
        type=Path,
        default="d28_clustering_analysis/metacells_with_clusters.h5ad",
        help="Path to the h5ad file with clustering results."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default="d28_clustering_analysis",
        help="Output directory for plots and results."
    )
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load the data
    adata = ad.read_h5ad(args.input)

    # Analyze composition for Phenograph clusters
    if 'phenograph' in adata.obs:
        analyze_composition(adata, 'phenograph', args.output_dir)

    # Analyze composition for Leiden clusters
    if 'leiden_optimal' in adata.obs:
        analyze_composition(adata, 'leiden_optimal', args.output_dir)

    print("Composition analysis complete.")

if __name__ == "__main__":
    main()
