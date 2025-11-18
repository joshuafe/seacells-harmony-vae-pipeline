
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import re

def plot_interstitial_proportions(stats_file, output_plot_file):
    """
    Generates a publication-ready bar plot of interstitial cell proportions by cohort.

    Args:
        stats_file: Path to the interstitial_cells_statistics.txt file.
        output_plot_file: Path to save the generated plot.
    """
    print("--- Generating interstitial cell proportions plot ---")

    # Read the statistics file
    with open(stats_file, 'r') as f:
        content = f.read()

    # Extract proportions
    proportions_match = re.search(r"Proportion of Interstitial Cells within each Cohort:\n(.*?)\n\n", content, re.DOTALL)
    if not proportions_match:
        raise ValueError("Could not find 'Proportion of Interstitial Cells within each Cohort' in the stats file.")
    
    proportions_str = proportions_match.group(1)
    proportions_df = pd.read_csv(pd.io.common.StringIO(proportions_str), sep=r'\s+', index_col=0)
    
    # Extract p-value
    p_value_match = re.search(r"P-value: ([\d.]+)", content)
    if not p_value_match:
        raise ValueError("Could not find 'P-value' in the stats file.")
    p_value = float(p_value_match.group(1))

    # Prepare data for plotting
    plot_data = proportions_df['True'].reset_index()
    plot_data.columns = ['Cohort', 'Proportion of Interstitial Cells']

    # Create the bar plot
    sns.set_style("whitegrid")
    plt.figure(figsize=(8, 6))
    ax = sns.barplot(x='Cohort', y='Proportion of Interstitial Cells', data=plot_data, palette='viridis')

    plt.title("Proportion of Interstitial Cells by Cohort", fontsize=16, fontweight='bold')
    plt.xlabel("Cohort", fontsize=14, fontweight='bold')
    plt.ylabel("Proportion of Interstitial Cells", fontsize=14, fontweight='bold')
    plt.ylim(0, plot_data['Proportion of Interstitial Cells'].max() * 1.2) # Adjust y-limit for p-value

    # Add p-value to the plot
    plt.text(0.05, 0.95, f"Chi-squared p-value: {p_value:.4f}", transform=ax.transAxes,
             fontsize=12, verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5))

    plt.tight_layout()
    plt.savefig(output_plot_file, dpi=300, bbox_inches='tight')
    print(f"  Saved plot to {output_plot_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Plot interstitial cell proportions.")
    parser.add_argument(
        "--stats-file",
        type=Path,
        default="d28_vae_analysis/interstitial_cells_statistics.txt",
        help="Path to the interstitial_cells_statistics.txt file."
    )
    parser.add_argument(
        "--output-plot-file",
        type=Path,
        default="d28_vae_analysis/interstitial_cells_proportions.png",
        help="Output file for the plot."
    )
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    args.output_plot_file.parent.mkdir(parents=True, exist_ok=True)

    plot_interstitial_proportions(args.stats_file, args.output_plot_file)

    print("\nInterstitial cell proportions plot generation complete.")

if __name__ == "__main__":
    main()
