
import pandas as pd
import anndata as ad
import plotly.graph_objects as go
from pathlib import Path
import argparse

def visualize_cluster_archetype_flow(cluster_adata_path, archetype_adata_path, output_file):
    """
    Creates a Sankey diagram to visualize the flow from clusters to archetypes.

    Args:
        cluster_adata_path: Path to the h5ad file with clustering results.
        archetype_adata_path: Path to the h5ad file with archetype results.
        output_file: Path to save the Sankey diagram HTML file.
    """
    print("--- Visualizing cluster to archetype flow ---")

    # Load the data
    cluster_adata = ad.read_h5ad(cluster_adata_path)
    archetype_adata = ad.read_h5ad(archetype_adata_path)

    # Merge the data
    df = pd.DataFrame({
        'cluster': cluster_adata.obs['phenograph'],
        'archetype': archetype_adata.obs['archetype_id_hard']
    })

    # Calculate the flow from clusters to archetypes
    flow = df.groupby(['cluster', 'archetype']).size().reset_index(name='value')

    # Create the Sankey diagram
    clusters = sorted(df['cluster'].unique())
    archetypes = sorted(df['archetype'].unique())

    nodes = [f"Cluster {c}" for c in clusters] + [f"Archetype {a}" for a in archetypes]
    node_map = {node: i for i, node in enumerate(nodes)}

    source = [node_map[f"Cluster {row['cluster']}"] for index, row in flow.iterrows()]
    target = [node_map[f"Archetype {row['archetype']}"] for index, row in flow.iterrows()]
    value = flow['value'].tolist()

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=nodes,
        ),
        link=dict(
            source=source,
            target=target,
            value=value
        )
    )])

    fig.update_layout(title_text="Flow from Phenograph Clusters to VAE Archetypes", font_size=10)
    fig.write_html(output_file)
    print(f"  Sankey diagram saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Visualize cluster to archetype flow.")
    parser.add_argument(
        "--cluster-input",
        type=Path,
        default="d28_clustering_analysis/metacells_with_clusters.h5ad",
        help="Path to the h5ad file with clustering results."
    )
    parser.add_argument(
        "--archetype-input",
        type=Path,
        default="d28_vae_analysis/metacells_with_archetypes.h5ad",
        help="Path to the h5ad file with archetype results."
    )
    parser.add_argument(
        "--output-file",
        type=Path,
        default="d28_vae_analysis/cluster_archetype_sankey.html",
        help="Output file for the Sankey diagram."
    )
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    args.output_file.parent.mkdir(parents=True, exist_ok=True)

    # Install plotly if not already installed
    try:
        import plotly
    except ImportError:
        print("Installing plotly...")
        import subprocess
        import sys
        subprocess.check_call([sys.executable, "-m", "pip", "install", "plotly"])

    visualize_cluster_archetype_flow(args.cluster_input, args.archetype_input, args.output_file)

    print("\nCluster to archetype flow visualization complete.")

if __name__ == "__main__":
    main()
