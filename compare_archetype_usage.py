
import pandas as pd
import anndata as ad
from pathlib import Path
import argparse
from scipy.stats import kruskal
import scikit_posthocs as sp

def compare_archetype_usage(adata, output_file):
    """
    Performs statistical comparison of archetype usage between cohorts.

    Args:
        adata: AnnData object with archetype probabilities.
        output_file: Path to save the statistical report.
    """
    print("--- Statistically comparing archetype usage ---")

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

    with open(output_file, "w") as report:
        report.write("Statistical Comparison of Archetype Usage by Cohort\n")
        report.write("=" * 60 + "\n\n")

        for archetype in archetype_cols:
            report.write(f"--- Archetype {archetype} ---")

            # Prepare data for Kruskal-Wallis test
            groups = [df[df['cohort'] == 'Normal'][archetype],
                      df[df['cohort'] == 'PTCy'][archetype],
                      df[df['cohort'] == 'Rux'][archetype]]

            # Kruskal-Wallis test
            h_stat, p_val = kruskal(*groups)
            report.write(f"  Kruskal-Wallis H-statistic: {h_stat:.4f}\n")
            report.write(f"  P-value: {p_val:.4f}\n")

            # Post-hoc Dunn's test if significant
            if p_val < 0.05:
                report.write("\n  Post-hoc Dunn's test (p-values):\n")
                dunn_results = sp.posthoc_dunn(df, val_col=archetype, group_col='cohort')
                report.write(dunn_results.to_string())
                report.write("\n")
            
            report.write("\n")

    print(f"Statistical report saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Statistically compare archetype usage.")
    parser.add_argument(
        "--input",
        type=Path,
        default="d28_vae_analysis/metacells_with_archetypes.h5ad",
        help="Path to the h5ad file with archetype results."
    )
    parser.add_argument(
        "--output-file",
        type=Path,
        default="d28_vae_analysis/archetype_usage_statistics.txt",
        help="Output file for the statistical report."
    )
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    args.output_file.parent.mkdir(parents=True, exist_ok=True)

    # Install scikit-posthocs if not already installed
    try:
        import scikit_posthocs
    except ImportError:
        print("Installing scikit-posthocs...")
        import subprocess
        import sys
        subprocess.check_call([sys.executable, "-m", "pip", "install", "scikit-posthocs"])

    # Load the data
    adata = ad.read_h5ad(args.input)

    compare_archetype_usage(adata, args.output_file)

    print("\nArchetype usage comparison complete.")

if __name__ == "__main__":
    main()
