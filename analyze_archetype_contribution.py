
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse

def analyze_archetype_contribution(adata, output_dir):
    """
    Analyzes the contribution of archetypes to each cohort.

    Args:
        adata: AnnData object with archetype probabilities.
        output_dir: Directory to save the plots.
    """
    print("--- Analyzing archetype contribution ---")

    # Create a DataFrame for analysis
    archetype_cols = [c for c in adata.obs.columns if 'archetype_' in c and '_prob' in c]
    df = adata.obs[archetype_cols + ['sample_id']].copy()

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

    # Calculate the mean archetype probability for each cohort
    cohort_archetype_usage = df.groupby('cohort')[archetype_cols].mean()

    # Save the usage table
    cohort_archetype_usage.to_csv(output_dir / "cohort_archetype_usage.csv")
    print(f"  Saved archetype usage table to {output_dir / 'cohort_archetype_usage.csv'}")

    # Plot the heatmap
    plt.figure(figsize=(12, 6))
    sns.heatmap(cohort_archetype_usage, cmap='viridis', annot=True, fmt=".2f")
    plt.title("Mean Archetype Usage by Cohort")
    plt.xlabel("Archetype")
    plt.ylabel("Cohort")
    plt.tight_layout()
    plt.savefig(output_dir / "cohort_archetype_usage.png", dpi=300)
    print(f"  Saved archetype usage heatmap to {output_dir / 'cohort_archetype_usage.png'}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Analyze archetype contribution by cohort.")
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

    analyze_archetype_contribution(adata, args.output_dir)

    print("Archetype contribution analysis complete.")

if __name__ == "__main__":
    main()
