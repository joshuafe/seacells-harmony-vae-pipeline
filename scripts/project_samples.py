"""
Project new samples onto trained VAE archetype model

This script takes a trained VAE model and projects new samples into the
learned latent space and assigns them to the established archetypes.

This is crucial for:
1. Evaluating new patient samples against the reference model
2. Comparing samples without retraining
3. Maintaining consistent archetype definitions
"""

import numpy as np
import pandas as pd
import scanpy as sc
import torch
import torch.nn as nn
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import cdist
import warnings
warnings.filterwarnings('ignore')

# Import model from training script
import sys
sys.path.append('scripts')


class VAEWithPriors(nn.Module):
    """Same as in vae_archetypes_with_priors.py"""

    def __init__(self, input_dim, latent_dim=5, hidden_dims=[64, 32], marker_names=None):
        super(VAEWithPriors, self).__init__()

        # Encoder
        encoder_layers = []
        prev_dim = input_dim
        for hidden_dim in hidden_dims:
            encoder_layers.extend([
                nn.Linear(prev_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.LeakyReLU(0.2),
                nn.Dropout(0.1)
            ])
            prev_dim = hidden_dim

        self.encoder = nn.Sequential(*encoder_layers)
        self.fc_mu = nn.Linear(prev_dim, latent_dim)
        self.fc_logvar = nn.Linear(prev_dim, latent_dim)

        # Decoder
        decoder_layers = []
        prev_dim = latent_dim
        for hidden_dim in reversed(hidden_dims):
            decoder_layers.extend([
                nn.Linear(prev_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.LeakyReLU(0.2),
                nn.Dropout(0.1)
            ])
            prev_dim = hidden_dim

        decoder_layers.append(nn.Linear(prev_dim, input_dim))
        self.decoder = nn.Sequential(*decoder_layers)

    def encode(self, x):
        h = self.encoder(x)
        return self.fc_mu(h), self.fc_logvar(h)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z):
        return self.decoder(z)

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        recon_x = self.decode(z)
        return recon_x, mu, logvar


def compute_soft_archetype_assignments(latent_embeddings, archetype_latent, temperature=1.0):
    """Compute soft assignments to archetypes"""
    distances = cdist(latent_embeddings, archetype_latent, metric='euclidean')
    similarities = -distances / temperature
    soft_assignments = np.exp(similarities) / np.exp(similarities).sum(axis=1, keepdims=True)
    entropy = -np.sum(soft_assignments * np.log(soft_assignments + 1e-10), axis=1)
    return soft_assignments, entropy


def load_trained_model(model_path, device='cpu'):
    """Load trained VAE model and associated data"""
    checkpoint = torch.load(model_path, map_location=device)

    # Extract config
    config = checkpoint['config']
    marker_names = checkpoint['marker_names']
    archetype_latent = checkpoint['archetype_latent']

    # Recreate model
    model = VAEWithPriors(
        input_dim=len(marker_names),
        latent_dim=config['latent_dim'],
        hidden_dims=[64, 32],
        marker_names=marker_names
    ).to(device)

    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()

    # Scaler
    scaler = StandardScaler()
    scaler.mean_ = checkpoint['scaler_mean']
    scaler.scale_ = checkpoint['scaler_scale']

    return model, scaler, archetype_latent, config


def project_samples(adata, model, scaler, archetype_latent, device='cpu', temperature=1.0):
    """
    Project samples onto trained model

    Args:
        adata: AnnData object with new samples
        model: trained VAE model
        scaler: fitted StandardScaler
        archetype_latent: archetype centers in latent space
        device: torch device
        temperature: softmax temperature for soft assignments

    Returns:
        adata: annotated with latent embeddings and archetype assignments
    """

    # Standardize
    X_scaled = scaler.transform(adata.X)

    # Convert to tensor
    X_tensor = torch.FloatTensor(X_scaled).to(device)

    # Encode
    with torch.no_grad():
        mu, logvar = model.encode(X_tensor)
        latent = mu.cpu().numpy()

    # Add to adata
    adata.obsm['X_vae'] = latent

    # Compute soft assignments
    soft_assign, entropy = compute_soft_archetype_assignments(
        latent, archetype_latent, temperature=temperature
    )

    # Add assignments to adata
    n_archetypes = archetype_latent.shape[0]
    for i in range(n_archetypes):
        adata.obs[f'archetype_{i}_prob'] = soft_assign[:, i]

    adata.obs['archetype_entropy'] = entropy
    adata.obs['archetype_id_hard'] = np.argmax(soft_assign, axis=1)

    return adata


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Project new samples onto trained VAE model")
    parser.add_argument(
        "--model",
        type=str,
        required=True,
        help="Path to trained model checkpoint (.pt file)"
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Input h5ad file with new samples to project"
    )
    parser.add_argument(
        "--reference-data",
        type=str,
        default=None,
        help="Optional: reference h5ad (for comparison visualization)"
    )
    parser.add_argument(
        "--soft-assignment-temp",
        type=float,
        default=1.0,
        help="Temperature for soft assignments"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="projection_output",
        help="Output directory"
    )

    args = parser.parse_args()

    print("="*80)
    print("PROJECT SAMPLES ONTO TRAINED VAE MODEL")
    print("="*80)
    print(f"Model: {args.model}")
    print(f"Input: {args.input}")
    print(f"Temperature: {args.soft_assignment_temp}")
    print("="*80)

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load trained model
    print("\nLoading trained model...")
    device = torch.device('mps' if torch.backends.mps.is_available() else 'cpu')
    model, scaler, archetype_latent, config = load_trained_model(args.model, device)
    print(f"  ✓ Loaded model (latent_dim={config['latent_dim']}, "
          f"{archetype_latent.shape[0]} archetypes)")

    # Load new data
    print("\nLoading new samples...")
    adata_new = sc.read_h5ad(args.input)
    print(f"  Loaded {adata_new.n_obs} metacells, {adata_new.n_vars} markers")

    # Project
    print("\nProjecting samples...")
    adata_new = project_samples(
        adata_new, model, scaler, archetype_latent,
        device=device, temperature=args.soft_assignment_temp
    )
    print(f"  ✓ Projected {adata_new.n_obs} metacells")
    print(f"  Entropy range: {adata_new.obs['archetype_entropy'].min():.3f} - "
          f"{adata_new.obs['archetype_entropy'].max():.3f}")

    # Compute UMAP
    print("\nComputing UMAP...")
    sc.pp.neighbors(adata_new, n_neighbors=15, use_rep='X_vae')
    sc.tl.umap(adata_new, min_dist=0.3)

    # Save results
    adata_new.write_h5ad(output_dir / "projected_samples.h5ad")
    print(f"  ✓ Saved: {output_dir / 'projected_samples.h5ad'}")

    # Generate visualizations
    print("\nGenerating visualizations...")

    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)

    # 1. Sample type
    if 'sample_type' in adata_new.obs:
        ax1 = fig.add_subplot(gs[0, 0])
        for sample_type in adata_new.obs['sample_type'].unique():
            mask = adata_new.obs['sample_type'] == sample_type
            ax1.scatter(adata_new.obsm['X_umap'][mask, 0],
                       adata_new.obsm['X_umap'][mask, 1],
                       label=sample_type, s=20, alpha=0.6)
        ax1.set_xlabel('UMAP 1')
        ax1.set_ylabel('UMAP 2')
        ax1.set_title('Sample Type')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

    # 2. Archetype assignment
    ax2 = fig.add_subplot(gs[0, 1])
    scatter = ax2.scatter(adata_new.obsm['X_umap'][:, 0],
                         adata_new.obsm['X_umap'][:, 1],
                         c=adata_new.obs['archetype_id_hard'],
                         s=20, alpha=0.6, cmap='tab20')
    plt.colorbar(scatter, ax=ax2, label='Archetype')
    ax2.set_xlabel('UMAP 1')
    ax2.set_ylabel('UMAP 2')
    ax2.set_title('Archetype Assignment')
    ax2.grid(True, alpha=0.3)

    # 3. Entropy
    ax3 = fig.add_subplot(gs[0, 2])
    scatter = ax3.scatter(adata_new.obsm['X_umap'][:, 0],
                         adata_new.obsm['X_umap'][:, 1],
                         c=adata_new.obs['archetype_entropy'],
                         s=20, alpha=0.6, cmap='viridis')
    plt.colorbar(scatter, ax=ax3, label='Entropy')
    ax3.set_xlabel('UMAP 1')
    ax3.set_ylabel('UMAP 2')
    ax3.set_title('Archetype Entropy')
    ax3.grid(True, alpha=0.3)

    # 4. Archetype composition by sample
    if 'sample_id' in adata_new.obs or 'sample_type' in adata_new.obs:
        ax4 = fig.add_subplot(gs[1, :])

        groupby_col = 'sample_id' if 'sample_id' in adata_new.obs else 'sample_type'

        composition = []
        n_archetypes = archetype_latent.shape[0]

        for sample in adata_new.obs[groupby_col].unique():
            mask = adata_new.obs[groupby_col] == sample
            row = {'Sample': sample}
            for i in range(n_archetypes):
                pct = (adata_new.obs.loc[mask, 'archetype_id_hard'] == i).sum() / mask.sum() * 100
                row[f'Archetype_{i}'] = pct
            composition.append(row)

        comp_df = pd.DataFrame(composition).set_index('Sample')

        sns.heatmap(comp_df, cmap='YlOrRd', annot=True, fmt='.1f',
                   cbar_kws={'label': 'Percentage'}, ax=ax4)
        ax4.set_title('Archetype Composition by Sample')
        ax4.set_xlabel('Archetype')
        ax4.set_ylabel('Sample')

        # Save composition
        comp_df.to_csv(output_dir / "archetype_composition.csv")

    plt.savefig(output_dir / "projection_overview.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved: {output_dir / 'projection_overview.png'}")

    # Summary statistics
    print("\n" + "="*80)
    print("PROJECTION SUMMARY")
    print("="*80)
    print(f"Total metacells projected: {adata_new.n_obs}")

    if 'sample_type' in adata_new.obs:
        print(f"\nSample type distribution:")
        print(adata_new.obs['sample_type'].value_counts())

    print(f"\nArchetype assignment distribution:")
    archetype_counts = adata_new.obs['archetype_id_hard'].value_counts().sort_index()
    for arch_id, count in archetype_counts.items():
        pct = count / adata_new.n_obs * 100
        print(f"  Archetype {arch_id}: {count} ({pct:.1f}%)")

    print(f"\nEntropy statistics:")
    print(f"  Mean: {adata_new.obs['archetype_entropy'].mean():.3f}")
    print(f"  Median: {adata_new.obs['archetype_entropy'].median():.3f}")
    print(f"  Range: {adata_new.obs['archetype_entropy'].min():.3f} - "
          f"{adata_new.obs['archetype_entropy'].max():.3f}")

    # High entropy metacells (mixed phenotypes)
    high_entropy_threshold = adata_new.obs['archetype_entropy'].quantile(0.9)
    n_mixed = (adata_new.obs['archetype_entropy'] > high_entropy_threshold).sum()
    print(f"\nMixed phenotype metacells (entropy > {high_entropy_threshold:.2f}): {n_mixed} ({n_mixed/adata_new.n_obs*100:.1f}%)")

    print("="*80)


if __name__ == "__main__":
    main()
