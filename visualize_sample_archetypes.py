
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse

def visualize_sample_archetypes(adata, output_dir):
    """
    Creates a heatmap of archetype usage for each individual sample.

    Args:
        adata: AnnData object with archetype probabilities.
        output_dir: Directory to save the plots.
    """
    print("--- Visualizing sample-level archetype usage ---")

    archetype_cols = [c for c in adata.obs.columns if 'archetype_' in c and '_prob' in c]
    df = adata.obs[archetype_cols + ['sample_id']].copy()

    # Calculate the mean archetype probability for each sample
    sample_archetype_usage = df.groupby('sample_id')[archetype_cols].mean()

    # Add cohort information for color-coding the sample labels
    def get_cohort(sample_id):
        if 'Normal' in sample_id:
            return 'Normal'
        elif 'PTCy' in sample_id:
            return 'PTCy'
        elif 'Abnormal' in sample_id:
            return 'Rux'
        return 'Unknown'

    sample_archetype_usage['cohort'] = sample_archetype_usage.index.map(get_cohort)
    sample_archetype_usage = sample_archetype_usage.sort_values('cohort')

    # Create cohort color mapping
    cohort_colors = sample_archetype_usage['cohort'].map({'Normal': 'g', 'PTCy': 'b', 'Rux': 'r'})

    # Create the heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(
        sample_archetype_usage[archetype_cols],
        cmap='viridis',
        annot=True,
        fmt=".2f",
        linewidths=.5
    )
    plt.title("Mean Archetype Usage by Sample")
    plt.xlabel("Archetype")
    plt.ylabel("Sample")
    
    # Color the y-axis labels by cohort
    for tick_label in plt.gca().get_yticklabels():
        tick_label.set_color(cohort_colors[tick_label.get_text()])

    plt.tight_layout()
    plt.savefig(output_dir / "sample_archetype_usage_heatmap.png", dpi=300)
    print(f"  Saved sample-level archetype usage heatmap to {output_dir / 'sample_archetype_usage_heatmap.png'}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Visualize sample-level archetype usage.")
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
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load the data
    adata = ad.read_h5ad(args.input)

    visualize_sample_archetypes(adata, args.output_dir)

    print("\nSample-level archetype visualization complete.")

if __name__ == "__main__":
    main()
