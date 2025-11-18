
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import numpy as np

def analyze_interstitial_cells(adata, output_dir, entropy_threshold_percentile=75):
    """
    Analyzes cells that are in the "interstitial" spaces between archetypes.

    Args:
        adata: AnnData object with archetype probabilities and entropy.
        output_dir: Directory to save the plots.
        entropy_threshold_percentile: The percentile of entropy to use as a threshold
                                      for defining interstitial cells.
    """
    print("--- Analyzing interstitial cells ---")

    # Calculate the entropy threshold
    entropy_threshold = np.percentile(adata.obs['archetype_entropy'], entropy_threshold_percentile)
    print(f"  Using entropy threshold: {entropy_threshold:.4f} ({entropy_threshold_percentile}th percentile)")

    # Identify interstitial cells
    adata.obs['is_interstitial'] = (adata.obs['archetype_entropy'] > entropy_threshold).astype(str)

    # --- UMAP Plot ---
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        x=adata.obsm['X_umap'][:, 0],
        y=adata.obsm['X_umap'][:, 1],
        hue=adata.obs['is_interstitial'],
        palette={'True': 'red', 'False': 'lightgray'},
        s=20,
        alpha=0.7
    )
    plt.title(f"Interstitial Cells (Entropy > {entropy_threshold_percentile}th percentile)")
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.legend(title="Interstitial")
    plt.tight_layout()
    plt.savefig(output_dir / "interstitial_cells_umap.png", dpi=300)
    print(f"  Saved interstitial cells UMAP plot to {output_dir / 'interstitial_cells_umap.png'}")
    plt.close()

    # --- Cohort Composition Analysis ---
    interstitial_cells = adata[adata.obs['is_interstitial'] == 'True'].copy()

    def get_cohort(sample_id):
        if 'Normal' in sample_id:
            return 'Normal'
        elif 'PTCy' in sample_id:
            return 'PTCy'
        elif 'Abnormal' in sample_id:
            return 'Rux'
        return 'Unknown'

    interstitial_cells.obs['cohort'] = interstitial_cells.obs['sample_id'].apply(get_cohort)

    cohort_composition = interstitial_cells.obs['cohort'].value_counts(normalize=True)
    print("\n  Cohort composition of interstitial cells:")
    print(cohort_composition)

    # Plot the composition
    plt.figure(figsize=(6, 6))
    cohort_composition.plot(kind='bar', color=['blue', 'orange', 'green'])
    plt.title("Cohort Composition of Interstitial Cells")
    plt.ylabel("Proportion")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_dir / "interstitial_cells_composition.png", dpi=300)
    print(f"  Saved interstitial cells composition plot to {output_dir / 'interstitial_cells_composition.png'}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Analyze interstitial cells between archetypes.")
    parser.add_argument(
        "--input",
        type=Path,
        default="d28_vae_analysis/metacells_with_archetypes.h5ad",
        help="Path to the h5ad file with archetype results."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default="d28_vae_analysis",
        help="Output directory for plots and results."
    )
    parser.add_argument(
        "--entropy-threshold",
        type=int,
        default=75,
        help="The percentile of entropy to use as a threshold for defining interstitial cells."
    )
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load the data
    adata = ad.read_h5ad(args.input)

    analyze_interstitial_cells(adata, args.output_dir, args.entropy_threshold)

    print("\nInterstitial cell analysis complete.")

if __name__ == "__main__":
    main()
