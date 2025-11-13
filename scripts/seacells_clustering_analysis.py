"""
PhenoGraph and Leiden Clustering Analysis of SEACells Metacells
Identifies optimal clusters and characterizes populations.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import scanpy as sc
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
import argparse

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 11

def ensure_fraction_abnormal(adata):
    """
    Guarantee that adata.obs contains a usable 'fraction_abnormal' column.

    Returns:
        available (bool): True if downstream analyses can use the column.
        note (str): Description of the data source or reason it's unavailable.
    """
    if 'fraction_abnormal' in adata.obs.columns:
        series = pd.to_numeric(adata.obs['fraction_abnormal'], errors='coerce')
        if series.notna().sum() == 0:
            return False, "fraction_abnormal column is present but empty."
        adata.obs['fraction_abnormal'] = series.clip(0, 1)
        return True, "Using existing fraction_abnormal values from the input AnnData."

    if 'sample_type' in adata.obs.columns:
        sample_type = adata.obs['sample_type'].astype(str).str.lower()
        allowed = {'normal', 'abnormal'}
        if set(sample_type.unique()).issubset(allowed):
            adata.obs['fraction_abnormal'] = (sample_type == 'abnormal').astype(float)
            return True, "Derived fraction_abnormal from binary sample_type labels."
        return False, ("sample_type exists but contains values beyond normal/abnormal; "
                       "cannot derive fraction_abnormal automatically.")

    return False, "No abnormality metadata available to derive fraction_abnormal."

def summarize_sample_types(adata, cluster_key, category_key='sample_type'):
    """Return a table of category counts per cluster for reporting."""
    if cluster_key not in adata.obs or category_key not in adata.obs:
        return None

    table = (adata.obs.groupby([cluster_key, category_key], observed=True)
                     .size()
                     .unstack(fill_value=0)
                     .astype(int))
    table['total'] = table.sum(axis=1)
    return table

import anndata

def load_metacells(input_path):
    """
    Load metacells AnnData object.

    Args:
        input_path: Path to integrated metacells h5ad file or directory with individual files

    Returns:
        AnnData object with metacells
    """
    input_path = Path(input_path)

    if input_path.is_file():
        # Single integrated file
        print(f"Loading integrated metacells from {input_path}...")
        adata = sc.read_h5ad(input_path)
        print(f"  Loaded {adata.n_obs} total metacells with {adata.n_vars} features")
        return adata
    elif input_path.is_dir():
        # Directory with multiple files - concatenate them
        h5ad_paths = list(input_path.glob('*.h5ad'))
        if not h5ad_paths:
            raise ValueError(f"No h5ad files found in {input_path}")

        print(f"Loading {len(h5ad_paths)} metacell files from directory...")
        adatas = []
        for path in h5ad_paths:
            print(f"  Loading {path.name}...")
            adata = sc.read_h5ad(path)
            adatas.append(adata)

        print("Concatenating AnnData objects...")
        concatenated_adata = anndata.concat(adatas, join='outer', merge='same')
        print(f"  Loaded {concatenated_adata.n_obs} total metacells with {concatenated_adata.n_vars} features")
        return concatenated_adata
    else:
        raise ValueError(f"Input path {input_path} is neither a file nor a directory")


def run_phenograph(adata, k=30, output_dir=None):
    """
    Run PhenoGraph clustering.

    Args:
        adata: AnnData object
        k: Number of nearest neighbors for PhenoGraph
        output_dir: Directory to save results
    """
    print(f"\n--- Running PhenoGraph Clustering ---")
    print(f"  k (nearest neighbors): {k}")

    try:
        import phenograph

        # Run PhenoGraph on the expression matrix
        X = adata.X
        if hasattr(X, 'toarray'):
            X = X.toarray()

        communities, graph, Q = phenograph.cluster(X, k=k)

        print(f"  Found {len(np.unique(communities))} communities")
        print(f"  Modularity Q: {Q:.4f}")

        # Add to adata
        adata.obs['phenograph'] = pd.Categorical(communities)

        # Calculate cluster sizes
        cluster_sizes = pd.Series(communities).value_counts().sort_index()
        print(f"  Cluster sizes:")
        for cluster, size in cluster_sizes.items():
            print(f"    Cluster {cluster}: {size} metacells")

        return adata, Q

    except ImportError:
        print("  WARNING: PhenoGraph not installed. Skipping PhenoGraph analysis.")
        print("  Install with: pip install phenograph")
        return adata, None

def run_leiden_clustering(adata, resolutions=[0.5, 1.0, 1.5, 2.0]):
    """
    Run Leiden clustering at multiple resolutions to find optimal.

    Args:
        adata: AnnData object
        resolutions: List of resolution parameters to test
    """
    print(f"\n--- Running Leiden Clustering ---")

    # Compute neighbors if not already done
    if 'neighbors' not in adata.uns:
        print("  Computing nearest neighbors...")
        sc.pp.neighbors(adata, n_neighbors=15, use_rep='X')

    results = []

    for res in resolutions:
        print(f"\n  Testing resolution: {res}")
        sc.tl.leiden(
            adata,
            resolution=res,
            key_added=f'leiden_{res}',
            n_iterations=2,
            directed=False
        )

        n_clusters = len(adata.obs[f'leiden_{res}'].unique())
        print(f"    Number of clusters: {n_clusters}")

        # Calculate clustering metrics
        X = adata.X
        if hasattr(X, 'toarray'):
            X = X.toarray()

        labels = pd.Categorical(adata.obs[f'leiden_{res}'])
        label_codes = labels.codes
        valid_mask = label_codes != -1

        if n_clusters > 1 and valid_mask.sum() > 1:
            silhouette = silhouette_score(X[valid_mask], label_codes[valid_mask])
            calinski = calinski_harabasz_score(X[valid_mask], label_codes[valid_mask])
            davies_bouldin = davies_bouldin_score(X[valid_mask], label_codes[valid_mask])

            print(f"    Silhouette score: {silhouette:.4f}")
            print(f"    Calinski-Harabasz score: {calinski:.2f}")
            print(f"    Davies-Bouldin score: {davies_bouldin:.4f}")

            results.append({
                'resolution': res,
                'n_clusters': n_clusters,
                'silhouette': silhouette,
                'calinski_harabasz': calinski,
                'davies_bouldin': davies_bouldin
            })
        else:
            print(f"    Skipping metrics (requires >=2 clusters with valid assignments)")

    results_df = pd.DataFrame(results)

    # Find optimal resolution (highest silhouette)
    if len(results_df) > 0:
        optimal_idx = results_df['silhouette'].idxmax()
        optimal_res = results_df.loc[optimal_idx, 'resolution']
        print(f"\n  Optimal resolution (by silhouette): {optimal_res}")
        print(f"    {results_df.loc[optimal_idx, 'n_clusters']} clusters")

        # Set optimal as main clustering
        adata.obs['leiden_optimal'] = adata.obs[f'leiden_{optimal_res}']
    else:
        optimal_res = None

    return adata, results_df, optimal_res

def plot_clustering_metrics(results_df, output_path):
    """Plot clustering quality metrics across resolutions."""
    if results_df.empty:
        print("  No Leiden results to plot (results dataframe is empty).")
        return

    print(f"\n--- Creating Clustering Metrics Plot ---")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Number of clusters
    ax = axes[0, 0]
    ax.plot(results_df['resolution'], results_df['n_clusters'], 'o-', linewidth=2, markersize=8, color='steelblue')
    ax.set_xlabel('Resolution', fontweight='bold')
    ax.set_ylabel('Number of Clusters', fontweight='bold')
    ax.set_title('Number of Clusters vs Resolution', fontweight='bold')
    ax.grid(alpha=0.3)

    # Silhouette score
    ax = axes[0, 1]
    ax.plot(results_df['resolution'], results_df['silhouette'], 'o-', linewidth=2, markersize=8, color='darkgreen')
    optimal_idx = results_df['silhouette'].idxmax()
    ax.axvline(results_df.loc[optimal_idx, 'resolution'], color='red', linestyle='--', linewidth=2, alpha=0.7)
    ax.set_xlabel('Resolution', fontweight='bold')
    ax.set_ylabel('Silhouette Score', fontweight='bold')
    ax.set_title('Silhouette Score (higher is better)', fontweight='bold')
    ax.grid(alpha=0.3)

    # Calinski-Harabasz score
    ax = axes[1, 0]
    ax.plot(results_df['resolution'], results_df['calinski_harabasz'], 'o-', linewidth=2, markersize=8, color='darkorange')
    ax.set_xlabel('Resolution', fontweight='bold')
    ax.set_ylabel('Calinski-Harabasz Score', fontweight='bold')
    ax.set_title('Calinski-Harabasz Score (higher is better)', fontweight='bold')
    ax.grid(alpha=0.3)

    # Davies-Bouldin score
    ax = axes[1, 1]
    ax.plot(results_df['resolution'], results_df['davies_bouldin'], 'o-', linewidth=2, markersize=8, color='darkred')
    ax.set_xlabel('Resolution', fontweight='bold')
    ax.set_ylabel('Davies-Bouldin Score', fontweight='bold')
    ax.set_title('Davies-Bouldin Score (lower is better)', fontweight='bold')
    ax.grid(alpha=0.3)

    plt.suptitle('Leiden Clustering Quality Metrics', fontweight='bold', fontsize=16, y=1.00)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  Saved clustering metrics plot to {output_path}")
    plt.close()

def plot_umap_clusters(adata, cluster_key, output_path, title=""):
    """Plot UMAP colored by cluster assignments."""
    print(f"\n--- Creating UMAP Plot for {cluster_key} ---")

    # Compute neighbors if not already done
    if 'neighbors' not in adata.uns:
        print("  Computing nearest neighbors...")
        sc.pp.neighbors(adata, n_neighbors=15, use_rep='X')

    # Compute UMAP if not already done
    if 'X_umap' not in adata.obsm:
        print("  Computing UMAP...")
        sc.tl.umap(adata, random_state=42)

    fig, ax = plt.subplots(figsize=(10, 8))

    sc.pl.umap(adata, color=cluster_key, ax=ax, show=False,
               title=title, legend_loc='right margin',
               frameon=False, size=100)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  Saved UMAP plot to {output_path}")
    plt.close()

def plot_cluster_composition(adata, cluster_key, output_path):
    """Plot composition of each cluster by normal/abnormal fraction."""
    print(f"\n--- Creating Cluster Composition Plot ---")

    # Create DataFrame
    df = pd.DataFrame({
        'cluster': adata.obs[cluster_key],
        'fraction_abnormal': adata.obs['fraction_abnormal']
    })

    # Group by cluster and get mean, std
    stats = df.groupby('cluster', observed=True)['fraction_abnormal'].agg(['mean', 'std', 'count']).reset_index()
    stats = stats.sort_values('mean', ascending=False)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Bar plot of mean fraction abnormal
    ax = axes[0]
    colors = plt.cm.RdYlGn_r(stats['mean'])
    bars = ax.bar(range(len(stats)), stats['mean'], yerr=stats['std'],
                  color=colors, edgecolor='black', linewidth=1.5, capsize=5)
    ax.set_xticks(range(len(stats)))
    ax.set_xticklabels(stats['cluster'])
    ax.set_xlabel('Cluster', fontweight='bold', fontsize=12)
    ax.set_ylabel('Mean Fraction Abnormal', fontweight='bold', fontsize=12)
    ax.set_title('Abnormal Cell Fraction by Cluster', fontweight='bold', pad=15)
    ax.axhline(0.5, color='black', linestyle='--', linewidth=2, alpha=0.5, label='Threshold (0.5)')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)

    # Violin plot
    ax = axes[1]
    parts = ax.violinplot([df[df['cluster'] == c]['fraction_abnormal'].values
                           for c in stats['cluster']],
                          positions=range(len(stats)),
                          widths=0.7,
                          showmeans=True,
                          showmedians=True)

    for pc in parts['bodies']:
        pc.set_alpha(0.7)

    ax.set_xticks(range(len(stats)))
    ax.set_xticklabels(stats['cluster'])
    ax.set_xlabel('Cluster', fontweight='bold', fontsize=12)
    ax.set_ylabel('Fraction Abnormal', fontweight='bold', fontsize=12)
    ax.set_title('Distribution of Abnormal Fraction by Cluster', fontweight='bold', pad=15)
    ax.axhline(0.5, color='black', linestyle='--', linewidth=2, alpha=0.5)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  Saved cluster composition plot to {output_path}")
    plt.close()

    return stats

def characterize_clusters(adata, cluster_key, output_path):
    """Characterize each cluster by marker expression."""
    if cluster_key not in adata.obs:
        print(f"  Cluster key '{cluster_key}' not found; skipping characterization.")
        return None

    print(f"\n--- Characterizing Clusters ---")

    # Get expression data
    X = adata.X
    if hasattr(X, 'toarray'):
        X = X.toarray()

    df = pd.DataFrame(X, columns=adata.var_names)
    df['cluster'] = adata.obs[cluster_key].values

    # Calculate mean expression per cluster
    cluster_means = df.groupby('cluster', observed=True).mean()

    # Z-score normalize across clusters for each marker (column-wise)
    from scipy import stats
    zscores = stats.zscore(cluster_means, axis=0, nan_policy='omit')
    cluster_means_zscore = pd.DataFrame(zscores, index=cluster_means.index, columns=cluster_means.columns).fillna(0)

    # Create heatmap
    fig, ax = plt.subplots(figsize=(12, max(8, len(cluster_means) * 0.5)))

    sns.heatmap(cluster_means_zscore, cmap='RdBu_r', center=0,
                cbar_kws={'label': 'Z-score (relative to other clusters)'},
                linewidths=0.5, linecolor='gray',
                annot=False, fmt='.2f', ax=ax, vmin=-2, vmax=2)

    ax.set_xlabel('Cluster', fontweight='bold', fontsize=12)
    ax.set_ylabel('Marker', fontweight='bold', fontsize=12)
    ax.set_title('Cluster Characterization by Marker Expression', fontweight='bold', pad=15)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  Saved cluster characterization heatmap to {output_path}")
    plt.close()

    # Save cluster profiles to CSV
    csv_path = output_path.parent / f'{output_path.stem}_profiles.csv'
    cluster_means.to_csv(csv_path)
    print(f"  Saved cluster profiles to {csv_path}")

    return cluster_means

def create_summary_report(adata, cluster_key, stats, results_df, output_path,
                          phenograph_Q=None, fraction_note=None,
                          sample_type_table=None):
    """Create summary report of clustering analysis."""
    print(f"\n--- Creating Summary Report ---")

    with open(output_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("CLUSTERING ANALYSIS SUMMARY\n")
        f.write("SEACells Metacells\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Total metacells analyzed: {adata.n_obs}\n")
        f.write(f"Number of markers: {adata.n_vars}\n\n")
        if fraction_note:
            f.write(f"Fraction abnormal note: {fraction_note}\n\n")

        if phenograph_Q is not None:
            f.write("PHENOGRAPH CLUSTERING\n")
            f.write("-" * 80 + "\n")
            n_phenograph = len(adata.obs['phenograph'].unique())
            f.write(f"Number of communities: {n_phenograph}\n")
            f.write(f"Modularity Q: {phenograph_Q:.4f}\n\n")

        f.write("LEIDEN CLUSTERING\n")
        f.write("-" * 80 + "\n")
        f.write(f"Resolutions tested: {list(results_df['resolution'])}\n\n")

        if not results_df.empty:
            f.write("Clustering Quality Metrics:\n")
            for _, row in results_df.iterrows():
                f.write(f"\nResolution {row['resolution']}:\n")
                f.write(f"  Number of clusters: {row['n_clusters']}\n")
                f.write(f"  Silhouette score: {row['silhouette']:.4f}\n")
                f.write(f"  Calinski-Harabasz score: {row['calinski_harabasz']:.2f}\n")
                f.write(f"  Davies-Bouldin score: {row['davies_bouldin']:.4f}\n")

            optimal_idx = results_df['silhouette'].idxmax()
            optimal_res = results_df.loc[optimal_idx, 'resolution']
            f.write(f"\nOptimal resolution (by silhouette): {optimal_res}\n")
            f.write(f"  Number of clusters: {results_df.loc[optimal_idx, 'n_clusters']}\n\n")
        else:
            f.write("No valid Leiden clustering metrics were computed.\n\n")

        f.write("CLUSTER COMPOSITION (Optimal Clustering)\n")
        f.write("-" * 80 + "\n")
        if stats is not None and not stats.empty:
            f.write(f"{'Cluster':<10} {'Size':<8} {'Mean Abn':<12} {'Std Abn':<12} {'Classification'}\n")
            f.write("-" * 80 + "\n")

            for _, row in stats.iterrows():
                classification = "Abnormal-Predominant" if row['mean'] >= 0.5 else "Normal-Predominant"
                f.write(f"{str(row['cluster']):<10} {int(row['count']):<8} "
                        f"{row['mean']:<12.4f} {row['std']:<12.4f} {classification}\n")
        else:
            f.write("Fraction abnormal data unavailable; cluster composition statistics omitted.\n")

        if sample_type_table is not None:
            f.write("\nSample type breakdown per cluster:\n")
            f.write(f"{'Cluster':<10}")
            for col in sample_type_table.columns:
                f.write(f"{col:<12}")
            f.write("\n" + "-" * 80 + "\n")
            for cluster, row in sample_type_table.iterrows():
                f.write(f"{str(cluster):<10}")
                for value in row:
                    f.write(f"{int(value):<12}")
                f.write("\n")

        f.write("\n" + "=" * 80 + "\n")

    print(f"  Saved summary report to {output_path}")

def main():
    parser = argparse.ArgumentParser(
        description="PhenoGraph and Leiden clustering analysis of SEACells metacells"
    )
    parser.add_argument(
        "--input",
        type=Path,
        default="seacells_output/batch_corrected/integrated_metacells_harmony.h5ad",
        help="Path to integrated metacells h5ad file or directory with individual metacell files"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default="seacells_output/clustering_analysis",
        help="Output directory for plots and results"
    )
    parser.add_argument(
        "--phenograph-k",
        type=int,
        default=30,
        help="Number of nearest neighbors for PhenoGraph"
    )
    parser.add_argument(
        "--leiden-resolutions",
        type=float,
        nargs='+',
        default=[0.5, 1.0, 1.5, 2.0],
        help="Leiden resolution parameters to test"
    )

    args = parser.parse_args()

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("CLUSTERING ANALYSIS OF SEACELLS METACELLS")
    print("=" * 80)

    # Load data
    adata = load_metacells(args.input)

    fraction_available, fraction_note = ensure_fraction_abnormal(adata)
    if fraction_available:
        print(f"Fraction abnormal information: {fraction_note}")
    else:
        print(f"WARNING: {fraction_note}")

    # Run PhenoGraph
    adata, phenograph_Q = run_phenograph(adata, k=args.phenograph_k)

    if phenograph_Q is not None:
        plot_umap_clusters(adata, 'phenograph',
                          args.output_dir / 'umap_phenograph.png',
                          title='PhenoGraph Clustering')

        if fraction_available:
            plot_cluster_composition(adata, 'phenograph',
                                    args.output_dir / 'phenograph_composition.png')
        else:
            print("  Skipping PhenoGraph composition plot due to missing fraction data.")

        characterize_clusters(adata, 'phenograph',
                            args.output_dir / 'phenograph_characterization.png')

    # Run Leiden clustering
    adata, results_df, optimal_res = run_leiden_clustering(adata, args.leiden_resolutions)

    # Save clustering metrics
    results_df.to_csv(args.output_dir / 'leiden_metrics.csv', index=False)
    print(f"\n  Saved Leiden metrics to {args.output_dir / 'leiden_metrics.csv'}")

    # Plot clustering metrics
    plot_clustering_metrics(results_df, args.output_dir / 'leiden_metrics.png')

    # Plot optimal clustering
    if optimal_res is not None:
        plot_umap_clusters(adata, 'leiden_optimal',
                          args.output_dir / 'umap_leiden_optimal.png',
                          title=f'Leiden Clustering (resolution={optimal_res})')

        stats = None
        if fraction_available:
            stats = plot_cluster_composition(adata, 'leiden_optimal',
                                            args.output_dir / 'leiden_composition.png')
        else:
            print("  Skipping Leiden composition plot due to missing fraction data.")

        characterize_clusters(adata, 'leiden_optimal',
                            args.output_dir / 'leiden_characterization.png')

        sample_type_table = summarize_sample_types(adata, 'leiden_optimal')

        # Create summary report
        create_summary_report(adata, 'leiden_optimal', stats, results_df,
                            args.output_dir / 'clustering_summary.txt',
                            phenograph_Q=phenograph_Q,
                            fraction_note=fraction_note,
                            sample_type_table=sample_type_table)

    # Save updated adata with clustering results
    adata.write(args.output_dir / 'metacells_with_clusters.h5ad')
    print(f"\n  Saved updated AnnData to {args.output_dir / 'metacells_with_clusters.h5ad'}")

    print("\n" + "=" * 80)
    print("CLUSTERING ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"\nAll outputs saved to: {args.output_dir}")

if __name__ == "__main__":
    main()
