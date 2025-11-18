
import pandas as pd
import argparse

def analyze_phenotypes(composition_file, profiles_file, output_file, top_n=5):
    """
    Analyzes the phenotypes of clusters based on their composition and marker expression.

    Args:
        composition_file: Path to the cluster composition CSV file.
        profiles_file: Path to the cluster characterization profiles CSV file.
        output_file: Path to save the summary report.
        top_n: The number of top markers to report for each cluster.
    """
    # Load the data
    composition_df = pd.read_csv(composition_file)
    profiles_df = pd.read_csv(profiles_file)

    # Determine the dominant cohort for each cluster
    def get_dominant_cohort(row):
        if row['Normal'] > 0.7:
            return 'Normal-dominant'
        if row['Rux'] > 0.7:
            return 'Rux-dominant'
        if row['PTCy'] > 0.7:
            return 'PTCy-dominant'
        return 'Mixed'

    composition_df['dominant_cohort'] = composition_df.apply(get_dominant_cohort, axis=1)

    # Get the top markers for each cluster
    profiles_df = profiles_df.set_index('cluster')
    top_markers = {}
    for cluster_id, row in profiles_df.iterrows():
        top_markers[cluster_id] = row.nlargest(top_n).index.tolist()

    # Create the summary report
    with open(output_file, 'w') as f:
        f.write("Cluster Phenotype Analysis Summary\n")
        f.write("=" * 40 + "\n\n")

        for cohort_type in ['Normal-dominant', 'Rux-dominant', 'PTCy-dominant', 'Mixed']:
            f.write(f"--- {cohort_type} Clusters ---\n")
            clusters_in_cohort = composition_df[composition_df['dominant_cohort'] == cohort_type]['phenograph'].tolist()

            if not clusters_in_cohort:
                f.write("None\n\n")
                continue

            for cluster_id in clusters_in_cohort:
                f.write(f"  Cluster {cluster_id}:\n")
                f.write(f"    Top Markers: {', '.join(top_markers.get(cluster_id, []))}\n")
            f.write("\n")

    print(f"Summary report saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Analyze cluster phenotypes.")
    parser.add_argument(
        "--composition-file",
        type=str,
        default="d28_clustering_analysis/phenograph_composition.csv",
        help="Path to the cluster composition CSV file."
    )
    parser.add_argument(
        "--profiles-file",
        type=str,
        default="d28_clustering_analysis/phenograph_characterization_profiles.csv",
        help="Path to the cluster characterization profiles CSV file."
    )
    parser.add_argument(
        "--output-file",
        type=str,
        default="d28_clustering_analysis/phenograph_phenotype_summary.txt",
        help="Path to save the summary report."
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=5,
        help="The number of top markers to report for each cluster."
    )
    args = parser.parse_args()

    analyze_phenotypes(args.composition_file, args.profiles_file, args.output_file, args.top_n)

if __name__ == "__main__":
    main()
