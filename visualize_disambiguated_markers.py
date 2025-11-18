
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import numpy as np

def visualize_disambiguated_markers(adata, profiles_file, output_dir):
    """
    Visualizes the expression of disambiguated markers on a UMAP plot.

    Args:
        adata: AnnData object with archetype probabilities and UMAP coordinates.
        profiles_file: Path to the archetype profiles (disambiguated) CSV file.
        output_dir: Directory to save the plots.
    """
    print("--- Visualizing disambiguated markers ---")

    profiles_df = pd.read_csv(profiles_file, index_col=0)
    archetype_cols = [c for c in adata.obs.columns if 'archetype_' in c and '_prob' in c]
    archetype_probs = adata.obs[archetype_cols].values

    # Calculate the disambiguated expression for each cell
    disambiguated_expression = archetype_probs @ profiles_df.values

    # Add the disambiguated expression to the AnnData object
    for i, marker in enumerate(profiles_df.columns):
        adata.obs[f"disambiguated_{marker}"] = disambiguated_expression[:, i]

    # Create a directory for the plots
    plots_dir = output_dir / "disambiguated_marker_umaps"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # Generate UMAP plots for each disambiguated marker
    for marker in profiles_df.columns:
        plt.figure(figsize=(8, 6))
        sns.scatterplot(
            x=adata.obsm['X_umap'][:, 0],
            y=adata.obsm['X_umap'][:, 1],
            hue=adata.obs[f"disambiguated_{marker}"],
            palette='viridis',
            s=10,
            linewidth=0,
            alpha=0.7
        )
        plt.title(f"Disambiguated Expression of {marker}")
        plt.xlabel("UMAP 1")
        plt.ylabel("UMAP 2")
        plt.legend(title="Expression")
        plt.tight_layout()
        plt.savefig(plots_dir / f"{marker}_umap.png", dpi=300)
        plt.close()

    print(f"  Saved UMAP plots to {plots_dir}")


def main():
    parser = argparse.ArgumentParser(description="Visualize disambiguated marker expression.")
    parser.add_argument(
        "--input",
        type=Path,
        default="d28_vae_analysis/metacells_with_archetypes.h5ad",
        help="Path to the h5ad file with archetype results."
    )
    parser.add_argument(
        "--profiles-file",
        type=Path,
        default="d28_vae_analysis/archetype_profiles_disambiguated.csv",
        help="Path to the archetype profiles (disambiguated) CSV file."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default="d28_vae_analysis",
        help="Output directory for plots and results."
    )
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load the data
    adata = ad.read_h5ad(args.input)

    visualize_disambiguated_markers(adata, args.profiles_file, args.output_dir)

    print("\nDisambiguated marker visualization complete.")

if __name__ == "__main__":
    main()
