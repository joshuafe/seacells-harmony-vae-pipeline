"""
Create comprehensive UMAP visualizations with marker expression overlays
and sample type (normal vs abnormal) comparisons.

Uses beautiful color palettes instead of viridis.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from pathlib import Path
import argparse

# Set style for publication-quality figures
sns.set_style("white")
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['font.size'] = 11
plt.rcParams['font.family'] = 'sans-serif'

def create_marker_umaps(adata, output_dir):
    """
    Create individual UMAP plots for each marker with expression overlay.
    Uses 'magma' colormap which is better than viridis for expression data.
    """
    print("\n" + "=" * 80)
    print("CREATING MARKER EXPRESSION UMAPS")
    print("=" * 80)

    output_dir = Path(output_dir)
    marker_dir = output_dir / "marker_umaps"
    marker_dir.mkdir(parents=True, exist_ok=True)

    markers = adata.var_names
    n_markers = len(markers)

    print(f"\nGenerating UMAPs for {n_markers} markers...")
    print(f"Using 'magma' colormap (dark purple → orange → yellow)")

    # Create individual marker plots
    for i, marker in enumerate(markers, 1):
        print(f"  [{i}/{n_markers}] {marker}")

        fig, ax = plt.subplots(figsize=(10, 8))

        # Get expression values
        expr = adata[:, marker].X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray().flatten()
        else:
            expr = np.array(expr).flatten()

        # Create scatter plot with magma colormap
        scatter = ax.scatter(
            adata.obsm['X_umap'][:, 0],
            adata.obsm['X_umap'][:, 1],
            c=expr,
            s=30,
            alpha=0.7,
            cmap='magma',  # Beautiful purple to yellow
            edgecolors='none',
            rasterized=True
        )

        # Colorbar
        cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
        cbar.set_label(f'{marker} Expression', fontweight='bold', fontsize=12)

        # Styling
        ax.set_xlabel('UMAP 1', fontweight='bold', fontsize=12)
        ax.set_ylabel('UMAP 2', fontweight='bold', fontsize=12)
        ax.set_title(f'{marker} Expression', fontweight='bold', fontsize=14, pad=15)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Remove tick labels for cleaner look
        ax.set_xticks([])
        ax.set_yticks([])

        plt.tight_layout()
        plt.savefig(marker_dir / f'umap_{marker}.png', dpi=300, bbox_inches='tight')
        plt.close()

    print(f"\n✓ Saved {n_markers} marker UMAPs to {marker_dir}/")
    return marker_dir

def create_marker_grid(adata, output_dir, markers_per_page=12):
    """
    Create grid layouts of marker UMAPs for easy comparison.
    """
    print("\n" + "=" * 80)
    print("CREATING MARKER GRID LAYOUTS")
    print("=" * 80)

    output_dir = Path(output_dir)

    markers = list(adata.var_names)
    n_markers = len(markers)
    n_pages = (n_markers + markers_per_page - 1) // markers_per_page

    print(f"\nCreating {n_pages} grid pages with up to {markers_per_page} markers each...")

    for page in range(n_pages):
        start_idx = page * markers_per_page
        end_idx = min(start_idx + markers_per_page, n_markers)
        page_markers = markers[start_idx:end_idx]
        n_markers_page = len(page_markers)

        # Calculate grid dimensions
        ncols = 4
        nrows = (n_markers_page + ncols - 1) // ncols

        fig, axes = plt.subplots(nrows, ncols, figsize=(20, 5 * nrows))
        axes = axes.flatten() if nrows > 1 else [axes] if ncols == 1 else axes

        print(f"\n  Page {page + 1}/{n_pages}: {page_markers}")

        for idx, marker in enumerate(page_markers):
            ax = axes[idx]

            # Get expression
            expr = adata[:, marker].X
            if hasattr(expr, 'toarray'):
                expr = expr.toarray().flatten()
            else:
                expr = np.array(expr).flatten()

            # Plot
            scatter = ax.scatter(
                adata.obsm['X_umap'][:, 0],
                adata.obsm['X_umap'][:, 1],
                c=expr,
                s=15,
                alpha=0.7,
                cmap='magma',
                edgecolors='none',
                rasterized=True
            )

            # Styling
            ax.set_title(marker, fontweight='bold', fontsize=12)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_visible(False)

            # Small colorbar
            plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)

        # Hide unused subplots
        for idx in range(n_markers_page, len(axes)):
            axes[idx].axis('off')

        plt.suptitle(f'Marker Expression Overview (Page {page + 1}/{n_pages})',
                     fontweight='bold', fontsize=16, y=0.995)
        plt.tight_layout()
        plt.savefig(output_dir / f'marker_grid_page{page + 1}.png', dpi=300, bbox_inches='tight')
        plt.close()

    print(f"\n✓ Saved {n_pages} marker grid pages")

def create_sample_type_umap(adata, output_dir):
    """
    Create UMAP showing Normal vs Abnormal samples with beautiful colors.
    """
    print("\n" + "=" * 80)
    print("CREATING SAMPLE TYPE UMAP (Normal vs Abnormal)")
    print("=" * 80)

    output_dir = Path(output_dir)

    # Define beautiful colors
    colors = {
        'Normal': '#2E7D32',    # Forest green
        'Abnormal': '#C62828'   # Deep red
    }

    fig, ax = plt.subplots(figsize=(12, 10))

    # Plot each sample type
    for sample_type in ['Normal', 'Abnormal']:
        mask = adata.obs['sample_type'] == sample_type
        n_samples = mask.sum()

        ax.scatter(
            adata.obsm['X_umap'][mask, 0],
            adata.obsm['X_umap'][mask, 1],
            c=colors[sample_type],
            label=f'{sample_type} (n={n_samples})',
            s=40,
            alpha=0.7,
            edgecolors='white',
            linewidth=0.5,
            rasterized=True
        )

    # Styling
    ax.set_xlabel('UMAP 1', fontweight='bold', fontsize=14)
    ax.set_ylabel('UMAP 2', fontweight='bold', fontsize=14)
    ax.set_title('Metacells by Sample Type', fontweight='bold', fontsize=16, pad=20)
    ax.legend(frameon=True, loc='upper right', fontsize=12, framealpha=0.95)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])

    plt.tight_layout()
    plt.savefig(output_dir / 'umap_sample_type.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n✓ Saved sample type UMAP")

def create_density_comparison(adata, output_dir):
    """
    Create density plots comparing Normal vs Abnormal distribution in UMAP space.
    """
    print("\n" + "=" * 80)
    print("CREATING DENSITY COMPARISON PLOTS")
    print("=" * 80)

    output_dir = Path(output_dir)

    fig, axes = plt.subplots(1, 3, figsize=(24, 7))

    # Get UMAP coordinates
    umap_coords = adata.obsm['X_umap']

    # Plot 1: Normal density
    ax = axes[0]
    mask_normal = adata.obs['sample_type'] == 'Normal'
    ax.hexbin(
        umap_coords[mask_normal, 0],
        umap_coords[mask_normal, 1],
        gridsize=50,
        cmap='Greens',
        mincnt=1,
        alpha=0.8
    )
    ax.set_title('Normal Samples Density', fontweight='bold', fontsize=14)
    ax.set_xlabel('UMAP 1', fontweight='bold')
    ax.set_ylabel('UMAP 2', fontweight='bold')
    ax.set_xticks([])
    ax.set_yticks([])

    # Plot 2: Abnormal density
    ax = axes[1]
    mask_abnormal = adata.obs['sample_type'] == 'Abnormal'
    ax.hexbin(
        umap_coords[mask_abnormal, 0],
        umap_coords[mask_abnormal, 1],
        gridsize=50,
        cmap='Reds',
        mincnt=1,
        alpha=0.8
    )
    ax.set_title('Abnormal Samples Density', fontweight='bold', fontsize=14)
    ax.set_xlabel('UMAP 1', fontweight='bold')
    ax.set_ylabel('UMAP 2', fontweight='bold')
    ax.set_xticks([])
    ax.set_yticks([])

    # Plot 3: Overlay with transparency
    ax = axes[2]
    ax.scatter(
        umap_coords[mask_normal, 0],
        umap_coords[mask_normal, 1],
        c='#2E7D32',
        s=30,
        alpha=0.5,
        label='Normal',
        edgecolors='none'
    )
    ax.scatter(
        umap_coords[mask_abnormal, 0],
        umap_coords[mask_abnormal, 1],
        c='#C62828',
        s=30,
        alpha=0.5,
        label='Abnormal',
        edgecolors='none'
    )
    ax.set_title('Overlay Comparison', fontweight='bold', fontsize=14)
    ax.set_xlabel('UMAP 1', fontweight='bold')
    ax.set_ylabel('UMAP 2', fontweight='bold')
    ax.legend(frameon=True, fontsize=11)
    ax.set_xticks([])
    ax.set_yticks([])

    plt.suptitle('Spatial Distribution: Normal vs Abnormal', fontweight='bold', fontsize=16, y=0.98)
    plt.tight_layout()
    plt.savefig(output_dir / 'density_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n✓ Saved density comparison plots")

def create_clustering_overlays(adata, output_dir):
    """
    Create UMAPs with clustering overlays using distinct colors.
    """
    print("\n" + "=" * 80)
    print("CREATING CLUSTERING OVERLAYS")
    print("=" * 80)

    output_dir = Path(output_dir)

    # Find clustering columns
    cluster_cols = [col for col in adata.obs.columns if 'leiden' in col or 'phenograph' in col]

    if not cluster_cols:
        print("  No clustering results found. Skipping.")
        return

    print(f"\nFound {len(cluster_cols)} clustering results: {cluster_cols}")

    for cluster_col in cluster_cols:
        fig, axes = plt.subplots(1, 2, figsize=(20, 8))

        # Left: Clusters with distinct colors
        ax = axes[0]

        # Get unique clusters and assign colors
        clusters = adata.obs[cluster_col].astype(str)
        unique_clusters = sorted(clusters.unique())
        n_clusters = len(unique_clusters)

        # Use tab20 for many clusters, Set3 for fewer
        if n_clusters <= 10:
            colors = plt.cm.Set3(np.linspace(0, 1, n_clusters))
        else:
            colors = plt.cm.tab20(np.linspace(0, 1, n_clusters))

        for i, cluster in enumerate(unique_clusters):
            mask = clusters == cluster
            ax.scatter(
                adata.obsm['X_umap'][mask, 0],
                adata.obsm['X_umap'][mask, 1],
                c=[colors[i]],
                label=f'Cluster {cluster}',
                s=40,
                alpha=0.7,
                edgecolors='white',
                linewidth=0.3
            )

        ax.set_title(f'{cluster_col.replace("_", " ").title()}', fontweight='bold', fontsize=14)
        ax.set_xlabel('UMAP 1', fontweight='bold')
        ax.set_ylabel('UMAP 2', fontweight='bold')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9, frameon=True)
        ax.set_xticks([])
        ax.set_yticks([])

        # Right: Clusters with sample type overlay
        ax = axes[1]

        # Plot with sample type colors
        for sample_type in ['Normal', 'Abnormal']:
            color = '#2E7D32' if sample_type == 'Normal' else '#C62828'
            mask = adata.obs['sample_type'] == sample_type
            ax.scatter(
                adata.obsm['X_umap'][mask, 0],
                adata.obsm['X_umap'][mask, 1],
                c=color,
                label=sample_type,
                s=40,
                alpha=0.6,
                edgecolors='white',
                linewidth=0.3
            )

        ax.set_title(f'{cluster_col.replace("_", " ").title()} + Sample Type',
                    fontweight='bold', fontsize=14)
        ax.set_xlabel('UMAP 1', fontweight='bold')
        ax.set_ylabel('UMAP 2', fontweight='bold')
        ax.legend(frameon=True, fontsize=11)
        ax.set_xticks([])
        ax.set_yticks([])

        plt.tight_layout()
        plt.savefig(output_dir / f'umap_{cluster_col}.png', dpi=300, bbox_inches='tight')
        plt.close()

        print(f"  ✓ Created {cluster_col} overlay")

    print(f"\n✓ Saved clustering overlay UMAPs")

def main():
    parser = argparse.ArgumentParser(
        description="Create comprehensive UMAP visualizations with marker overlays"
    )
    parser.add_argument(
        "--input",
        type=Path,
        default="seacells_output/clustering_analysis/metacells_with_clusters.h5ad",
        help="Path to clustered metacells h5ad file"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default="seacells_output/clustering_analysis",
        help="Output directory for plots"
    )
    parser.add_argument(
        "--markers-per-page",
        type=int,
        default=12,
        help="Number of markers per grid page"
    )

    args = parser.parse_args()

    print("=" * 80)
    print("COMPREHENSIVE UMAP VISUALIZATION GENERATOR")
    print("=" * 80)
    print(f"\nInput: {args.input}")
    print(f"Output: {args.output_dir}")

    # Load data
    print("\nLoading data...")
    adata = sc.read_h5ad(args.input)
    print(f"  Loaded {adata.n_obs} metacells with {adata.n_vars} markers")
    print(f"  Sample types: {adata.obs['sample_type'].value_counts().to_dict()}")

    # Ensure UMAP exists
    if 'X_umap' not in adata.obsm:
        print("\nERROR: No UMAP coordinates found. Run clustering analysis first.")
        return

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Generate all visualizations
    create_sample_type_umap(adata, args.output_dir)
    create_density_comparison(adata, args.output_dir)
    create_clustering_overlays(adata, args.output_dir)
    create_marker_umaps(adata, args.output_dir)
    create_marker_grid(adata, args.output_dir, args.markers_per_page)

    print("\n" + "=" * 80)
    print("VISUALIZATION COMPLETE!")
    print("=" * 80)
    print(f"\nAll plots saved to: {args.output_dir}")
    print("\nGenerated:")
    print("  • umap_sample_type.png - Normal vs Abnormal comparison")
    print("  • density_comparison.png - Spatial density analysis")
    print("  • umap_[clustering].png - Clustering overlays")
    print("  • marker_umaps/ - Individual marker expression UMAPs")
    print("  • marker_grid_page*.png - Grid layouts for quick comparison")
    print("=" * 80)

if __name__ == "__main__":
    main()
