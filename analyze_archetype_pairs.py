
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import numpy as np
from collections import Counter
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch

def chord_diagram(matrix, names, ax=None, colors=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.set_aspect('equal')

    n = len(matrix)
    if colors is None:
        colors = plt.cm.viridis(np.linspace(0, 1, n))

    # Create the circle
    circle = plt.Circle((0, 0), 1, color='white', ec='black', fill=False)
    ax.add_artist(circle)

    # Add nodes
    node_pos = {}
    for i in range(n):
        angle = (np.pi / 2) - (2 * np.pi * i / n)
        x = np.cos(angle)
        y = np.sin(angle)
        node_pos[names[i]] = (x, y)
        ax.text(x * 1.1, y * 1.1, names[i], ha='center', va='center', fontsize=12)

    # Add edges (chords)
    for i in range(n):
        for j in range(i + 1, n):
            if matrix[i, j] > 0:
                start_node = node_pos[names[i]]
                end_node = node_pos[names[j]]
                
                path_data = [
                    (start_node),
                    ((start_node[0] + end_node[0]) / 2, (start_node[1] + end_node[1]) / 2),
                    (end_node)
                ]
                
                path = MplPath(path_data, [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3])
                patch = PathPatch(path, facecolor='none', ec='black', lw=matrix[i, j] / 100)
                ax.add_patch(patch)

    ax.autoscale_view()
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.axis('off')

def analyze_archetype_pairs(adata, output_dir):
    """
    Analyzes the relationships between pairs of archetypes.

    Args:
        adata: AnnData object with archetype probabilities.
        output_dir: Directory to save the plots.
    """
    print("--- Analyzing archetype pairs and 2D distributions ---")

    archetype_cols = [c for c in adata.obs.columns if 'archetype_' in c and '_prob' in c]
    archetype_probs = adata.obs[archetype_cols]

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

    # Find the top two archetypes for each cell
    top2_archetypes_indices = np.argsort(-archetype_probs.values, axis=1)[:, :2]
    
    # Create a DataFrame to store top 2 archetype probabilities and indices
    top2_probs_df = pd.DataFrame({
        'prob1': [archetype_probs.iloc[i, top2_archetypes_indices[i, 0]] for i in range(len(adata))],
        'prob2': [archetype_probs.iloc[i, top2_archetypes_indices[i, 1]] for i in range(len(adata))],
        'idx1': top2_archetypes_indices[:, 0],
        'idx2': top2_archetypes_indices[:, 1],
        'entropy': adata.obs['archetype_entropy'],
        'cohort': adata.obs['cohort']
    }, index=adata.obs.index)

    # Filter for cells where prob1 + prob2 > 0 to avoid division by zero
    top2_probs_df = top2_probs_df[top2_probs_df['prob1'] + top2_probs_df['prob2'] > 0].copy()
    top2_probs_df['proximity_ratio'] = top2_probs_df['prob1'] / (top2_probs_df['prob1'] + top2_probs_df['prob2'])

    # --- Analyze common pairs ---
    # Sort the archetype indices within each pair to count unique pairs regardless of order
    sorted_pairs = [tuple(sorted((idx1, idx2))) for idx1, idx2 in zip(top2_probs_df['idx1'], top2_probs_df['idx2'])]
    pair_counts = Counter(sorted_pairs)
    common_pairs = pair_counts.most_common(10)

    print("\n  Top 10 most common archetype pairs:")
    for pair, count in common_pairs:
        print(f"    Archetypes {pair[0]} and {pair[1]}: {count} cells")

    # --- Visualize 2D Distribution for each common pair ---
    plots_dir = output_dir / "archetype_pair_distributions"
    plots_dir.mkdir(parents=True, exist_ok=True)

    for pair, _ in common_pairs:
        a1, a2 = pair[0], pair[1]
        
        # Select cells belonging to this pair
        pair_data = top2_probs_df[
            ((top2_probs_df['idx1'] == a1) & (top2_probs_df['idx2'] == a2)) |
            ((top2_probs_df['idx1'] == a2) & (top2_probs_df['idx2'] == a1))
        ].copy()

        # Ensure proximity ratio is always for a1 vs a2
        pair_data['proximity_ratio_a1'] = pair_data.apply(
            lambda row: row['prob1'] / (row['prob1'] + row['prob2']) if row['idx1'] == a1 else row['prob2'] / (row['prob1'] + row['prob2']),
            axis=1
        )

        plt.figure(figsize=(10, 8))
        sns.scatterplot(
            x='proximity_ratio_a1',
            y='entropy',
            hue='cohort',
            data=pair_data,
            s=20,
            alpha=0.6,
            edgecolor=None,
            palette='viridis'
        )
        plt.title(f"Distribution of Cells Between Archetypes {a1} and {a2}", fontsize=14, fontweight='bold')
        plt.xlabel(f"Proximity to Archetype {a1} (1.0) vs {a2} (0.0)", fontsize=12, fontweight='bold')
        plt.ylabel("Archetype Entropy", fontsize=12, fontweight='bold')
        plt.xlim(-0.05, 1.05)
        plt.tight_layout()
        plot_path = plots_dir / f"archetype_pair_{a1}_{a2}_2d_distribution.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved 2D distribution plot for Archetypes {a1} and {a2} to {plot_path}")

    print("\nArchetype pair and 2D distribution analysis complete.")


def main():
    parser = argparse.ArgumentParser(description="Analyze archetype pair relationships.")
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

    analyze_archetype_pairs(adata, args.output_dir)

    print("\nArchetype pair analysis complete.")

if __name__ == "__main__":
    main()
