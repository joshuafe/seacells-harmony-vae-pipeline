"""
VAE Archetypes with Biological Priors and Train-Project Framework

Key features:
1. Train on reference samples (Normal) and project others
2. Biological prior loss functions for lineage constraints
3. Soft archetype assignments to quantify mixed phenotypes
4. Hyperparameter tuning support
"""

import numpy as np
import pandas as pd
import scanpy as sc
import torch
import torch.nn as nn
import torch.nn.functional as F
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import cdist
import warnings
warnings.filterwarnings('ignore')

# Set random seeds
np.random.seed(42)
torch.manual_seed(42)


class BiologicalPriorLoss(nn.Module):
    """
    Biological prior loss to enforce known lineage relationships

    Despite merged channels (CD3+CD19, CD4+CD33, etc.), we can still
    enforce soft constraints on expected co-expression patterns
    """

    def __init__(self, marker_names, weights=None):
        super().__init__()
        self.marker_names = list(marker_names)
        self.marker_idx = {name: i for i, name in enumerate(marker_names)}

        # Default weights for each prior type
        if weights is None:
            weights = {
                'tcell_nk_mutual_exclusion': 1.0,
                'cd4_cd8_preference': 0.5,
                'memory_naive_exclusion': 0.3,
            }
        self.weights = weights

    def forward(self, x):
        """
        Compute biological prior loss

        Args:
            x: reconstructed expression (batch_size, n_markers)

        Returns:
            loss: scalar prior loss
        """
        total_loss = 0.0

        # 1. T-cell vs NK mutual exclusion
        # CD3+CD19 should be high XOR CD16+TIGIT/CD56 should be high
        if 'CD3+CD19' in self.marker_idx and 'CD16+TIGIT' in self.marker_idx and 'CD56' in self.marker_idx:
            cd3_cd19 = x[:, self.marker_idx['CD3+CD19']]
            cd16_tigit = x[:, self.marker_idx['CD16+TIGIT']]
            cd56 = x[:, self.marker_idx['CD56']]

            # NK score: high CD16 or CD56
            nk_score = torch.max(cd16_tigit, cd56)

            # Penalize high CD3 AND high NK markers (should be mutually exclusive)
            mutual_exclusion = torch.relu(cd3_cd19) * torch.relu(nk_score)
            total_loss += self.weights['tcell_nk_mutual_exclusion'] * mutual_exclusion.mean()

        # 2. CD4 vs CD8 preference (most T cells should be one or the other)
        # CD4+CD33 vs CD8+CD14
        if 'CD4+CD33' in self.marker_idx and 'CD8+CD14' in self.marker_idx and 'CD3+CD19' in self.marker_idx:
            cd3_cd19 = x[:, self.marker_idx['CD3+CD19']]
            cd4_cd33 = x[:, self.marker_idx['CD4+CD33']]
            cd8_cd14 = x[:, self.marker_idx['CD8+CD14']]

            # For cells with high CD3, penalize intermediate CD4/CD8
            # (they should be clearly CD4+ or CD8+)
            tcell_mask = torch.sigmoid(cd3_cd19)  # soft mask for T cells
            cd4_cd8_product = torch.abs(cd4_cd33) * torch.abs(cd8_cd14)

            # Penalize high product (both high or both low is bad)
            cd4_cd8_loss = tcell_mask * cd4_cd8_product
            total_loss += self.weights['cd4_cd8_preference'] * cd4_cd8_loss.mean()

        # 3. Memory vs Naive exclusion
        # CD45RA (naive) vs CD45RO (memory) should be anti-correlated
        if 'CD45RA' in self.marker_idx and 'CD45RO' in self.marker_idx:
            cd45ra = x[:, self.marker_idx['CD45RA']]
            cd45ro = x[:, self.marker_idx['CD45RO']]

            # Penalize both being high
            memory_naive_conflict = torch.relu(cd45ra) * torch.relu(cd45ro)
            total_loss += self.weights['memory_naive_exclusion'] * memory_naive_conflict.mean()

        return total_loss


class VAEWithPriors(nn.Module):
    """Variational Autoencoder with biological priors"""

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

        # Biological prior loss
        if marker_names is not None:
            self.bio_prior_loss = BiologicalPriorLoss(marker_names)
        else:
            self.bio_prior_loss = None

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


def vae_loss_with_priors(recon_x, x, mu, logvar, model, beta=1.0, gamma=1.0):
    """
    VAE loss with reconstruction, KL divergence, and biological priors

    Args:
        recon_x: reconstructed input
        x: original input
        mu: latent mean
        logvar: latent log variance
        model: VAE model (to access bio_prior_loss)
        beta: weight for KL divergence
        gamma: weight for biological priors
    """
    # Reconstruction loss
    recon_loss = nn.MSELoss()(recon_x, x)

    # KL divergence
    kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
    kl_loss /= x.size(0) * x.size(1)

    # Biological prior loss
    if model.bio_prior_loss is not None and gamma > 0:
        bio_loss = model.bio_prior_loss(recon_x)
    else:
        bio_loss = torch.tensor(0.0, device=recon_x.device)

    total_loss = recon_loss + beta * kl_loss + gamma * bio_loss

    return total_loss, recon_loss, kl_loss, bio_loss


def train_vae(model, data_loader, optimizer, device, beta=1.0, gamma=1.0):
    """Train VAE for one epoch"""
    model.train()
    total_loss = 0
    total_recon = 0
    total_kl = 0
    total_bio = 0

    for batch_tuple in data_loader:
        batch_data = batch_tuple[0].to(device)
        optimizer.zero_grad()

        recon_batch, mu, logvar = model(batch_data)
        loss, recon_loss, kl_loss, bio_loss = vae_loss_with_priors(
            recon_batch, batch_data, mu, logvar, model, beta, gamma
        )

        loss.backward()
        optimizer.step()

        total_loss += loss.item()
        total_recon += recon_loss.item()
        total_kl += kl_loss.item()
        total_bio += bio_loss.item() if isinstance(bio_loss, torch.Tensor) else bio_loss

    n_batches = len(data_loader)
    return (total_loss / n_batches, total_recon / n_batches,
            total_kl / n_batches, total_bio / n_batches)


def compute_soft_archetype_assignments(latent_embeddings, archetype_latent, temperature=1.0):
    """
    Compute soft assignments to archetypes using softmax over distances

    Args:
        latent_embeddings: (n_cells, latent_dim) embeddings
        archetype_latent: (n_archetypes, latent_dim) archetype centers
        temperature: softmax temperature (lower = sharper assignments)

    Returns:
        soft_assignments: (n_cells, n_archetypes) probability distribution
        entropy: (n_cells,) assignment entropy (measures "mixedness")
    """
    # Compute distances
    distances = cdist(latent_embeddings, archetype_latent, metric='euclidean')

    # Convert to similarities (negative distance)
    similarities = -distances / temperature

    # Softmax to get probabilities
    soft_assignments = np.exp(similarities) / np.exp(similarities).sum(axis=1, keepdims=True)

    # Compute entropy (measure of uncertainty/mixedness)
    # High entropy = mixed phenotype, low entropy = pure archetype
    entropy = -np.sum(soft_assignments * np.log(soft_assignments + 1e-10), axis=1)

    return soft_assignments, entropy


def identify_archetypes(latent_embeddings, n_archetypes=10, method='kmeans'):
    """Identify archetype metacells in latent space"""
    if method == 'kmeans':
        from sklearn.cluster import KMeans
        kmeans = KMeans(n_clusters=n_archetypes, random_state=42, n_init=20)
        kmeans.fit(latent_embeddings)
        # Find closest real metacells to cluster centers
        archetypes = []
        for center in kmeans.cluster_centers_:
            distances = np.linalg.norm(latent_embeddings - center, axis=1)
            archetypes.append(np.argmin(distances))
        return archetypes

    elif method == 'extremes':
        archetypes = []
        for dim in range(min(n_archetypes, latent_embeddings.shape[1])):
            min_idx = np.argmin(latent_embeddings[:, dim])
            max_idx = np.argmax(latent_embeddings[:, dim])
            archetypes.extend([min_idx, max_idx])
        return list(set(archetypes))[:n_archetypes]

    else:
        raise ValueError(f"Unknown method: {method}")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="VAE archetypes with biological priors and train-project framework"
    )
    parser.add_argument(
        "--input",
        type=str,
        default="seacells_output/three_groups_phenograph/integrated_metacells_harmony_phenograph.h5ad",
        help="Input h5ad file"
    )
    parser.add_argument(
        "--train-on-reference",
        action="store_true",
        help="Train only on Normal samples, then project others"
    )
    parser.add_argument(
        "--reference-type",
        type=str,
        default="Normal",
        help="Sample type to use as reference (default: Normal)"
    )
    parser.add_argument(
        "--latent-dim",
        type=int,
        default=5,
        help="Latent dimension size"
    )
    parser.add_argument(
        "--n-archetypes",
        type=int,
        default=10,
        help="Number of archetypes to identify"
    )
    parser.add_argument(
        "--archetype-method",
        type=str,
        default="kmeans",
        choices=["extremes", "kmeans"],
        help="Method for identifying archetypes"
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=300,
        help="Number of training epochs"
    )
    parser.add_argument(
        "--beta",
        type=float,
        default=0.01,
        help="Beta parameter for KL divergence weighting"
    )
    parser.add_argument(
        "--gamma",
        type=float,
        default=0.1,
        help="Gamma parameter for biological prior weighting (0 = no priors)"
    )
    parser.add_argument(
        "--beta-warmup",
        type=int,
        default=100,
        help="Number of epochs to warm up beta"
    )
    parser.add_argument(
        "--gamma-cooldown",
        type=int,
        default=150,
        help="Number of epochs to cool down gamma (biological priors). "
             "Gamma starts at full strength and decreases to gamma_min. "
             "0 = no cooldown (constant gamma)"
    )
    parser.add_argument(
        "--gamma-min",
        type=float,
        default=0.0,
        help="Minimum gamma value after cooldown (default: 0.0)"
    )
    parser.add_argument(
        "--soft-assignment-temp",
        type=float,
        default=1.0,
        help="Temperature for soft archetype assignments (lower = sharper)"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="vae_archetypes_with_priors_output",
        help="Output directory"
    )

    args = parser.parse_args()

    print("="*80)
    print("VAE ARCHETYPES WITH BIOLOGICAL PRIORS")
    print("="*80)
    print(f"Input: {args.input}")
    print(f"Train on reference: {args.train_on_reference}")
    if args.train_on_reference:
        print(f"  Reference type: {args.reference_type}")
    print(f"Latent dimensions: {args.latent_dim}")
    print(f"Number of archetypes: {args.n_archetypes}")
    print(f"Biological prior weight (gamma): {args.gamma}")
    if args.gamma_cooldown > 0:
        print(f"  Gamma cooldown: {args.gamma} → {args.gamma_min} over {args.gamma_cooldown} epochs")
    print(f"Soft assignment temperature: {args.soft_assignment_temp}")
    print("="*80)

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    print("\nLoading data...")
    adata = sc.read_h5ad(args.input)
    print(f"  Loaded {adata.n_obs} metacells, {adata.n_vars} markers")
    print(f"  Sample types: {adata.obs['sample_type'].value_counts().to_dict()}")

    # Prepare training and projection sets
    if args.train_on_reference:
        print(f"\n** TRAIN-AND-PROJECT MODE **")
        train_mask = adata.obs['sample_type'] == args.reference_type
        adata_train = adata[train_mask].copy()
        adata_project = adata[~train_mask].copy()

        print(f"  Training set ({args.reference_type}): {adata_train.n_obs} metacells")
        print(f"  Projection set (others): {adata_project.n_obs} metacells")
        print(f"    {adata_project.obs['sample_type'].value_counts().to_dict()}")

        X_train = adata_train.X.copy()
        X_project = adata_project.X.copy()
    else:
        print(f"\n** TRAIN-ON-ALL MODE **")
        adata_train = adata.copy()
        X_train = adata_train.X.copy()
        X_project = None

    # Standardize
    print("\nStandardizing data...")
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)

    if X_project is not None:
        X_project_scaled = scaler.transform(X_project)  # Use same scaler

    # Convert to tensors
    X_train_tensor = torch.FloatTensor(X_train_scaled)
    train_dataset = torch.utils.data.TensorDataset(X_train_tensor)
    train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=64, shuffle=True)

    # Setup device
    device = torch.device('mps' if torch.backends.mps.is_available() else 'cpu')
    print(f"  Using device: {device}")

    # Create VAE model
    print(f"\nCreating VAE model with biological priors...")
    print(f"  Input dim: {adata.n_vars}, Latent dim: {args.latent_dim}")
    print(f"  Biological prior weight (gamma): {args.gamma}")

    model = VAEWithPriors(
        input_dim=adata.n_vars,
        latent_dim=args.latent_dim,
        hidden_dims=[64, 32],
        marker_names=adata.var_names
    ).to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=5e-4, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='min', factor=0.5, patience=20
    )

    # Training
    print(f"\nTraining VAE for {args.epochs} epochs...")
    print(f"  Beta warmup: 0 → {args.beta} over {args.beta_warmup} epochs")
    if args.gamma_cooldown > 0:
        print(f"  Gamma cooldown: {args.gamma} → {args.gamma_min} over {args.gamma_cooldown} epochs")
    train_losses = []
    recon_losses = []
    kl_losses = []
    bio_losses = []
    gammas = []

    for epoch in range(args.epochs):
        # Beta warmup schedule
        if epoch < args.beta_warmup:
            current_beta = args.beta * (epoch / args.beta_warmup)
        else:
            current_beta = args.beta

        # Gamma cooldown schedule
        # Start high (strong priors), gradually decrease (relax constraints)
        if args.gamma_cooldown > 0:
            if epoch < args.gamma_cooldown:
                # Linear cooldown from gamma to gamma_min
                progress = epoch / args.gamma_cooldown
                current_gamma = args.gamma * (1 - progress) + args.gamma_min * progress
            else:
                current_gamma = args.gamma_min
        else:
            # No cooldown, constant gamma
            current_gamma = args.gamma

        gammas.append(current_gamma)

        loss, recon, kl, bio = train_vae(
            model, train_loader, optimizer, device,
            beta=current_beta, gamma=current_gamma
        )

        train_losses.append(loss)
        recon_losses.append(recon)
        kl_losses.append(kl)
        bio_losses.append(bio)

        scheduler.step(loss)

        if (epoch + 1) % 50 == 0:
            print(f"  Epoch {epoch+1}/{args.epochs}: Loss={loss:.4f}, "
                  f"Recon={recon:.4f}, KL={kl:.4f}, Bio={bio:.4f}, Gamma={current_gamma:.4f}")

    print("  ✓ Training complete")

    # Get latent embeddings for training data
    print("\nComputing latent embeddings...")
    model.eval()
    with torch.no_grad():
        X_train_device = X_train_tensor.to(device)
        mu_train, _ = model.encode(X_train_device)
        latent_train = mu_train.cpu().numpy()

    adata_train.obsm['X_vae'] = latent_train

    # Identify archetypes on training data
    print(f"\nIdentifying {args.n_archetypes} archetypes from training data...")
    archetype_indices = identify_archetypes(
        latent_train,
        n_archetypes=args.n_archetypes,
        method=args.archetype_method
    )

    # Get archetype centers
    archetype_latent = latent_train[archetype_indices]

    # Compute soft assignments for training data
    print("  Computing soft archetype assignments...")
    soft_assign_train, entropy_train = compute_soft_archetype_assignments(
        latent_train, archetype_latent, temperature=args.soft_assignment_temp
    )

    # Add to adata
    for i in range(args.n_archetypes):
        adata_train.obs[f'archetype_{i}_prob'] = soft_assign_train[:, i]
    adata_train.obs['archetype_entropy'] = entropy_train
    adata_train.obs['archetype_id_hard'] = np.argmax(soft_assign_train, axis=1)

    # Mark archetype centers
    adata_train.obs['is_archetype'] = False
    adata_train.obs.loc[adata_train.obs.index[archetype_indices], 'is_archetype'] = True

    print(f"  ✓ Identified {len(archetype_indices)} archetypes")
    print(f"  ✓ Computed soft assignments (entropy range: {entropy_train.min():.3f} - {entropy_train.max():.3f})")

    # Project other samples if in train-and-project mode
    if args.train_on_reference and X_project is not None:
        print(f"\nProjecting {adata_project.n_obs} non-reference metacells...")
        X_project_tensor = torch.FloatTensor(X_project_scaled)

        with torch.no_grad():
            X_project_device = X_project_tensor.to(device)
            mu_project, _ = model.encode(X_project_device)
            latent_project = mu_project.cpu().numpy()

        adata_project.obsm['X_vae'] = latent_project

        # Compute soft assignments for projected data
        soft_assign_project, entropy_project = compute_soft_archetype_assignments(
            latent_project, archetype_latent, temperature=args.soft_assignment_temp
        )

        for i in range(args.n_archetypes):
            adata_project.obs[f'archetype_{i}_prob'] = soft_assign_project[:, i]
        adata_project.obs['archetype_entropy'] = entropy_project
        adata_project.obs['archetype_id_hard'] = np.argmax(soft_assign_project, axis=1)
        adata_project.obs['is_archetype'] = False

        print(f"  ✓ Projected samples")
        print(f"  ✓ Entropy range: {entropy_project.min():.3f} - {entropy_project.max():.3f}")

        # Combine for visualization
        adata_combined = sc.concat([adata_train, adata_project], join='outer')
        adata_combined.obs['is_reference'] = adata_combined.obs['sample_type'] == args.reference_type
    else:
        adata_combined = adata_train
        adata_combined.obs['is_reference'] = True

    # UMAP
    print("\nComputing UMAP...")
    sc.pp.neighbors(adata_combined, n_neighbors=15, use_rep='X_vae')
    sc.tl.umap(adata_combined, min_dist=0.3)

    # Save model
    print("\nSaving model...")
    torch.save({
        'model_state_dict': model.state_dict(),
        'scaler_mean': scaler.mean_,
        'scaler_scale': scaler.scale_,
        'archetype_latent': archetype_latent,
        'archetype_indices': archetype_indices,
        'marker_names': list(adata.var_names),
        'config': vars(args)
    }, output_dir / "vae_model.pt")
    print("  ✓ Saved model checkpoint")

    # Save results
    adata_combined.write_h5ad(output_dir / "metacells_with_archetypes.h5ad")
    print("  ✓ Saved annotated data")

    # Export archetype profiles
    archetype_data = adata_train[adata_train.obs['is_archetype']].copy()
    archetype_profiles = pd.DataFrame(
        archetype_data.X,
        index=[f"Archetype_{i}" for i in range(len(archetype_indices))],
        columns=adata.var_names
    )
    archetype_profiles.to_csv(output_dir / "archetype_profiles.csv")

    # Generate visualizations
    print("\nGenerating visualizations...")
    generate_visualizations(
        adata_combined, archetype_indices, archetype_profiles,
        train_losses, recon_losses, kl_losses, bio_losses, gammas,
        output_dir, args
    )

    print("\n" + "="*80)
    print("COMPLETE")
    print("="*80)


def generate_visualizations(adata, archetype_indices, archetype_profiles,
                           train_losses, recon_losses, kl_losses, bio_losses, gammas,
                           output_dir, args):
    """Generate comprehensive visualizations"""

    # 1. Training curves
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    axes[0, 0].plot(train_losses)
    axes[0, 0].set_title('Total Loss')
    axes[0, 0].set_xlabel('Epoch')
    axes[0, 0].grid(True, alpha=0.3)

    axes[0, 1].plot(recon_losses)
    axes[0, 1].set_title('Reconstruction Loss')
    axes[0, 1].set_xlabel('Epoch')
    axes[0, 1].grid(True, alpha=0.3)

    axes[0, 2].plot(kl_losses)
    axes[0, 2].set_title('KL Divergence')
    axes[0, 2].set_xlabel('Epoch')
    axes[0, 2].grid(True, alpha=0.3)

    axes[1, 0].plot(bio_losses)
    axes[1, 0].set_title('Biological Prior Loss')
    axes[1, 0].set_xlabel('Epoch')
    axes[1, 0].grid(True, alpha=0.3)

    # Gamma schedule plot
    axes[1, 1].plot(gammas, color='purple', linewidth=2)
    axes[1, 1].set_title('Gamma Schedule\n(Biological Prior Weight)')
    axes[1, 1].set_xlabel('Epoch')
    axes[1, 1].set_ylabel('Gamma')
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].axhline(y=args.gamma, color='purple', linestyle='--', alpha=0.3, label=f'Initial γ={args.gamma}')
    if args.gamma_cooldown > 0:
        axes[1, 1].axhline(y=args.gamma_min, color='purple', linestyle='--', alpha=0.3, label=f'Final γ={args.gamma_min}')
        axes[1, 1].axvline(x=args.gamma_cooldown, color='red', linestyle=':', alpha=0.3, label=f'Cooldown end')
    axes[1, 1].legend(fontsize=8)

    # Combined loss components
    epochs = np.arange(len(train_losses))
    axes[1, 2].plot(epochs, recon_losses, label='Reconstruction', alpha=0.7)
    axes[1, 2].plot(epochs, kl_losses, label='KL Divergence', alpha=0.7)
    axes[1, 2].plot(epochs, bio_losses, label='Biological Prior', alpha=0.7)
    axes[1, 2].set_title('Loss Components')
    axes[1, 2].set_xlabel('Epoch')
    axes[1, 2].set_ylabel('Loss')
    axes[1, 2].legend()
    axes[1, 2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / "training_curves.png", dpi=300, bbox_inches='tight')
    plt.close()

    # 2. UMAP visualizations
    fig = plt.figure(figsize=(20, 12))
    gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)

    # Sample type
    ax1 = fig.add_subplot(gs[0, 0])
    for sample_type, color in zip(['Normal', 'Abnormal', 'PTCy'],
                                   ['#1f77b4', '#ff7f0e', '#2ca02c']):
        mask = adata.obs['sample_type'] == sample_type
        if mask.sum() > 0:
            ax1.scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
                       c=color, label=sample_type, s=20, alpha=0.6, edgecolors='none')

    # Highlight archetypes
    if 'is_archetype' in adata.obs:
        arch_mask = adata.obs['is_archetype']
        ax1.scatter(adata.obsm['X_umap'][arch_mask, 0], adata.obsm['X_umap'][arch_mask, 1],
                   c='red', s=200, marker='*', edgecolors='black', linewidths=2,
                   label='Archetypes', zorder=10)
    ax1.set_xlabel('UMAP 1')
    ax1.set_ylabel('UMAP 2')
    ax1.set_title('Sample Type')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Archetype entropy (mixedness)
    ax2 = fig.add_subplot(gs[0, 1])
    scatter = ax2.scatter(adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1],
                         c=adata.obs['archetype_entropy'], s=20, alpha=0.6,
                         cmap='viridis', edgecolors='none')
    plt.colorbar(scatter, ax=ax2, label='Entropy (mixedness)')
    ax2.set_xlabel('UMAP 1')
    ax2.set_ylabel('UMAP 2')
    ax2.set_title('Archetype Entropy\n(high = mixed phenotype)')
    ax2.grid(True, alpha=0.3)

    # Hard archetype assignment
    ax3 = fig.add_subplot(gs[0, 2])
    scatter = ax3.scatter(adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1],
                         c=adata.obs['archetype_id_hard'], s=20, alpha=0.6,
                         cmap='tab20', edgecolors='none')
    plt.colorbar(scatter, ax=ax3, label='Archetype ID')
    ax3.set_xlabel('UMAP 1')
    ax3.set_ylabel('UMAP 2')
    ax3.set_title('Hard Archetype Assignment')
    ax3.grid(True, alpha=0.3)

    # Archetype heatmap
    ax4 = fig.add_subplot(gs[1, :])
    sns.heatmap(archetype_profiles, cmap='RdBu_r', center=0,
                xticklabels=True, yticklabels=True,
                cbar_kws={'label': 'Expression'}, ax=ax4)
    ax4.set_title('Archetype Marker Profiles')
    ax4.set_xlabel('Markers')
    ax4.set_ylabel('Archetypes')

    plt.savefig(output_dir / "vae_archetypes_overview.png", dpi=300, bbox_inches='tight')
    plt.close()

    print("  ✓ Saved visualizations")


if __name__ == "__main__":
    main()
