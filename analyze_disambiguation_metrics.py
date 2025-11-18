
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import argparse

def analyze_disambiguation(profiles_file, output_dir):
    """
    Analyzes the VAE's disambiguation of merged channels.

    Args:
        profiles_file: Path to the archetype profiles (disambiguated) CSV file.
        output_dir: Directory to save the plots and report.
    """
    print("--- Analyzing channel disambiguation ---")

    profiles_df = pd.read_csv(profiles_file, index_col=0)

    merged_to_unmerged = {
        'CD16+TIGIT': ['CD16', 'TIGIT'],
        'CD8+CD14': ['CD8', 'CD14'],
        'CD4+CD33': ['CD4', 'CD33'],
        'CD279+CD24': ['CD279', 'CD24'],
        'CD34+CD223': ['CD34', 'CD223'],
        'CD3+CD19': ['CD3', 'CD19'],
    }

    with open(output_dir / "disambiguation_report.txt", "w") as report:
        report.write("Channel Disambiguation Analysis Report\n")
        report.write("=" * 40 + "\n\n")

        for merged, unmerged in merged_to_unmerged.items():
            marker1, marker2 = unmerged[0], unmerged[1]
            
            if marker1 not in profiles_df.columns or marker2 not in profiles_df.columns:
                print(f"  Skipping {merged}: Markers not found in profiles.")
                continue

            # --- Correlation ---
            correlation = profiles_df[marker1].corr(profiles_df[marker2])
            report.write(f"--- {merged} -> {marker1} vs {marker2} ---")
            report.write(f"  Correlation: {correlation:.4f}\n")

            # --- Scatter Plot ---
            plt.figure(figsize=(6, 6))
            sns.scatterplot(data=profiles_df, x=marker1, y=marker2)
            plt.title(f"Disambiguation of {merged}")
            plt.xlabel(marker1)
            plt.ylabel(marker2)
            plt.grid(True)
            plt.tight_layout()
            plot_path = output_dir / f"disambiguation_{marker1}_vs_{marker2}.png"
            plt.savefig(plot_path, dpi=300)
            plt.close()
            report.write(f"  Scatter plot saved to: {plot_path}\n\n")

    print(f"Disambiguation analysis complete. Report saved to {output_dir / 'disambiguation_report.txt'}")


def main():
    parser = argparse.ArgumentParser(description="Analyze VAE channel disambiguation.")
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

    analyze_disambiguation(args.profiles_file, args.output_dir)


if __name__ == "__main__":
    main()
