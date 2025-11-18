
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import re
import numpy as np

def plot_archetype_usage_by_cohort(adata, stats_file, output_dir):
    """
    Generates publication-ready bar plots for each archetype's usage by cohort,
    with error bars and statistical significance markers.

    Args:
        adata: AnnData object with archetype probabilities.
        stats_file: Path to the archetype_usage_statistics.txt file.
        output_dir: Directory to save the generated plots.
    """
    print("--- Generating archetype usage plots by cohort ---")

    archetype_cols = [c for c in adata.obs.columns if 'archetype_' in c and '_prob' in c]
    df = adata.obs[archetype_cols + ['sample_id']].copy()

    # Map sample_id to cohort
    def get_cohort(sample_id):
        if 'Normal' in sample_id:
            return 'Normal'
        elif 'PTCy' in sample_id:
            return 'PTCy'
        elif 'Abnormal' in sample_id:
            return 'Rux'
        return 'Unknown'

    df['cohort'] = df['sample_id'].apply(get_cohort)

    # Read statistical results
    with open(stats_file, 'r') as f:
        stats_content = f.read()

    # Prepare output directory for plots
    plots_dir = output_dir / "archetype_usage_plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    sns.set_style("whitegrid")

    for archetype_col in archetype_cols:
        archetype_id = archetype_col.split('_')[1] # e.g., 'archetype_0_prob' -> '0'
        
        # Extract p-values for this archetype
        archetype_pattern = re.compile(rf"--- Archetype {re.escape(archetype_col)} ---\n.*?P-value: ([\d.]+)", re.DOTALL)
        match = archetype_pattern.search(stats_content)
        kruskal_p_value = float(match.group(1)) if match else 1.0

        dunn_p_values = {}
        dunn_pattern = re.compile(rf"--- Archetype {re.escape(archetype_col)} ---\n.*?Post-hoc Dunn's test \(p-values\):\n(.*?)\n\n", re.DOTALL)
        dunn_match = dunn_pattern.search(stats_content)
        if dunn_match:
            dunn_table_str = dunn_match.group(1)
            dunn_df = pd.read_csv(pd.io.common.StringIO(dunn_table_str), sep=r'\s+', index_col=0)
            for c1 in dunn_df.index:
                for c2 in dunn_df.columns:
                    if c1 != c2:
                        dunn_p_values[tuple(sorted((c1, c2)))] = dunn_df.loc[c1, c2]

        # Create the bar plot
        plt.figure(figsize=(6, 5))
        ax = sns.barplot(x='cohort', y=archetype_col, data=df, errorbar='sd', capsize=0.1, palette='viridis')
        sns.stripplot(x='cohort', y=archetype_col, data=df, color='black', size=3, jitter=0.2, ax=ax)

        plt.title(f"Archetype {archetype_id} Usage by Cohort", fontsize=14, fontweight='bold')
        plt.xlabel("Cohort", fontsize=12, fontweight='bold')
        plt.ylabel(f"Mean Probability of Archetype {archetype_id}", fontsize=12, fontweight='bold')
        plt.ylim(0, df[archetype_col].max() * 1.2) # Adjust y-limit for significance stars

        # Add significance stars
        if kruskal_p_value < 0.05:
            y_max = df[archetype_col].max()
            h = y_max * 0.05
            
            cohorts = ['Normal', 'PTCy', 'Rux']
            for i, c1 in enumerate(cohorts):
                for j, c2 in enumerate(cohorts):
                    if i < j:
                        p_val = dunn_p_values.get(tuple(sorted((c1, c2))), 1.0)
                        if p_val < 0.001:
                            stars = '***'
                        elif p_val < 0.01:
                            stars = '**'
                        elif p_val < 0.05:
                            stars = '*'
                        else:
                            stars = ''
                        
                        if stars:
                            x1, x2 = i, j
                            y = y_max * 1.05
                            ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='k')
                            ax.text((x1+x2)*.5, y+h, stars, ha='center', va='bottom', color='k', fontsize=10)

        plt.tight_layout()
        plot_path = plots_dir / f"archetype_{archetype_id}_usage.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved plot for Archetype {archetype_id} to {plot_path}")

    print("\nArchetype usage plots generation complete.")


def main():
    parser = argparse.ArgumentParser(description="Plot archetype usage by cohort.")
    parser.add_argument(
        "--input",
        type=Path,
        default="d28_vae_analysis/metacells_with_archetypes.h5ad",
        help="Path to the h5ad file with archetype results."
    )
    parser.add_argument(
        "--stats-file",
        type=Path,
        default="d28_vae_analysis/archetype_usage_statistics.txt",
        help="Path to the archetype_usage_statistics.txt file."
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

    plot_archetype_usage_by_cohort(adata, args.stats_file, args.output_dir)

if __name__ == "__main__":
    main()
