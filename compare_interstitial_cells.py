
import pandas as pd
import anndata as ad
from pathlib import Path
import argparse
import numpy as np
from scipy.stats import chi2_contingency

def compare_interstitial_cells(adata, output_file, entropy_threshold_percentile=75):
    """
    Performs a statistical comparison of the proportion of interstitial cells
    between cohorts.

    Args:
        adata: AnnData object with archetype probabilities and entropy.
        output_file: Path to save the statistical report.
        entropy_threshold_percentile: The percentile of entropy to use as a threshold
                                      for defining interstitial cells.
    """
    print("--- Statistically comparing interstitial cells ---")

    # Calculate the entropy threshold
    entropy_threshold = np.percentile(adata.obs['archetype_entropy'], entropy_threshold_percentile)
    adata.obs['is_interstitial'] = (adata.obs['archetype_entropy'] > entropy_threshold)

    # Map sample_id to cohort
    def get_cohort(sample_id):
        if 'Normal' in sample_id:
            return 'Normal'
        elif 'PTCy' in sample_id:
            return 'PTCy'
        elif 'Abnormal' in sample_id:
            return 'Rux'
        return 'Unknown'

    adata.obs['cohort'] = adata.obs['sample_id'].apply(get_cohort)

    # Create a contingency table
    contingency_table = pd.crosstab(adata.obs['cohort'], adata.obs['is_interstitial'])

    # Calculate proportions
    proportions = contingency_table.div(contingency_table.sum(axis=1), axis=0)

    # Perform the chi-squared test
    chi2, p, dof, expected = chi2_contingency(contingency_table)

    # Generate the report
    with open(output_file, "w") as report:
        report.write("Statistical Comparison of Interstitial Cells by Cohort\n")
        report.write("=" * 60 + "\n\n")
        report.write(f"Entropy Threshold: {entropy_threshold:.4f} ({entropy_threshold_percentile}th percentile)\n\n")
        report.write("Contingency Table:\n")
        report.write(contingency_table.to_string())
        report.write("\n\n")
        report.write("Proportion of Interstitial Cells within each Cohort:\n")
        report.write(proportions.to_string())
        report.write("\n\n")
        report.write("Chi-Squared Test Results:\n")
        report.write(f"  Chi-squared statistic: {chi2:.4f}\n")
        report.write(f"  P-value: {p:.4f}\n")
        report.write(f"  Degrees of freedom: {dof}\n")
        report.write("\nExpected Frequencies:\n")
        report.write(str(expected))

    print(f"Statistical report saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Statistically compare interstitial cells.")
    parser.add_argument(
        "--input",
        type=Path,
        default="d28_vae_analysis/metacells_with_archetypes.h5ad",
        help="Path to the h5ad file with archetype results."
    )
    parser.add_argument(
        "--output-file",
        type=Path,
        default="d28_vae_analysis/interstitial_cells_statistics.txt",
        help="Output file for the statistical report."
    )
    parser.add_argument(
        "--entropy-threshold",
        type=int,
        default=75,
        help="The percentile of entropy to use as a threshold for defining interstitial cells."
    )
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    args.output_file.parent.mkdir(parents=True, exist_ok=True)

    # Load the data
    adata = ad.read_h5ad(args.input)

    compare_interstitial_cells(adata, args.output_file, args.entropy_threshold)

    print("\nInterstitial cell comparison complete.")

if __name__ == "__main__":
    main()
