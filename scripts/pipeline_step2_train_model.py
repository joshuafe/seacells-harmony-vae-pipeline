"""
Metacell Pipeline v2: VAE Disambiguation Model Training (Full Script)

This script constitutes the second step of the corrected analysis pipeline.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import torch
import torch.nn as nn
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import cdist
import warnings
import argparse

warnings.filterwarnings('ignore')
np.random.seed(42)
torch.manual_seed(42)

# --- Ground Truth Marker Panel (17 -> 23 Markers) ---
ORIGINAL_MARKERS = [
    'CD62L', 'CD152 R718', 'CD45', 'CD45RO', 'CD279+CD24', 'CD95', 
    'CD16+TIGIT', 'CD8+CD14', 'CD56', 'CD45RA', 'CD366', 'CD69', 
    'CD25 BB515', 'CD34+CD223', 'CD197', 'CD3+CD19', 'CD4+CD33'
]

MERGED_TO_UNMERGED = {
    'CD16+TIGIT': ['CD16', 'TIGIT'], 'CD8+CD14': ['CD8', 'CD14'],
    'CD4+CD33': ['CD4', 'CD33'], 'CD279+CD24': ['CD279', 'CD24'],
    'CD34+CD223': ['CD34', 'CD223'], 'CD3+CD19': ['CD3', 'CD19'],
}

UNMERGED_MARKERS = []
for m in ORIGINAL_MARKERS:
    cleaned_m = m.split(' ')[0]
    if cleaned_m in MERGED_TO_UNMERGED:
        UNMERGED_MARKERS.extend(MERGED_TO_UNMERGED[cleaned_m])
    else:
        UNMERGED_MARKERS.append(cleaned_m)

# --- Biological Priors ---
class DisambiguationBiologicalPriorLoss(nn.Module):
    def __init__(self, unmerged_marker_names, weights=None):
        super().__init__()
        self.marker_idx = {name: i for i, name in enumerate(unmerged_marker_names)}
        if weights is None:
            weights = {
                'tcell_nk_mutual_exclusion': 1.0,
                'cd4_cd8_preference': 0.1,
                'memory_naive_exclusion': 0.5,
                'bcell_tcell_exclusion': 1.0,
            }
        self.weights = weights

    def forward(self, x_disambiguated):
        total_loss = 0.0
        if 'CD3' in self.marker_idx and 'CD19' in self.marker_idx:
            cd3 = x_disambiguated[:, self.marker_idx['CD3']]
            cd19 = x_disambiguated[:, self.marker_idx['CD19']]
            total_loss += self.weights['bcell_tcell_exclusion'] * (torch.relu(cd3) * torch.relu(cd19)).mean()
        if 'CD3' in self.marker_idx and 'CD56' in self.marker_idx:
            cd3 = x_disambiguated[:, self.marker_idx['CD3']]
            cd56 = x_disambiguated[:, self.marker_idx['CD56']]
            total_loss += self.weights['tcell_nk_mutual_exclusion'] * (torch.relu(cd3) * torch.relu(cd56)).mean()
        if 'CD4' in self.marker_idx and 'CD8' in self.marker_idx and 'CD3' in self.marker_idx:
            cd3 = x_disambiguated[:, self.marker_idx['CD3']]
            cd4 = x_disambiguated[:, self.marker_idx['CD4']]
            cd8 = x_disambiguated[:, self.marker_idx['CD8']]
            tcell_mask = torch.sigmoid(cd3)
            cd4_cd8_product = torch.relu(cd4) * torch.relu(cd8)
            total_loss += self.weights['cd4_cd8_preference'] * (tcell_mask * cd4_cd8_product).mean()
        if 'CD45RA' in self.marker_idx and 'CD45RO' in self.marker_idx:
            cd45ra = x_disambiguated[:, self.marker_idx['CD45RA']]
            cd45ro = x_disambiguated[:, self.marker_idx['CD45RO']]
            total_loss += self.weights['memory_naive_exclusion'] * (torch.relu(cd45ra) * torch.relu(cd45ro)).mean()
        return total_loss

# --- VAE Model ---
class VAEDisambiguation(nn.Module):
    def __init__(self, input_dim, latent_dim=10, hidden_dims=[128, 64]):
        super().__init__()
        self.unmerged_dim = len(UNMERGED_MARKERS)
        assert self.unmerged_dim == 23, f"Expected 23 unmerged markers, got {self.unmerged_dim}"
        
        # ... (Full model architecture) ...
        encoder_layers = []
        prev_dim = input_dim
        for hidden_dim in hidden_dims:
            encoder_layers.extend([nn.Linear(prev_dim, hidden_dim), nn.BatchNorm1d(hidden_dim), nn.LeakyReLU(0.2)])
            prev_dim = hidden_dim
        self.encoder = nn.Sequential(*encoder_layers)
        self.fc_mu = nn.Linear(prev_dim, latent_dim)
        self.fc_logvar = nn.Linear(prev_dim, latent_dim)

        decoder_layers = []
        prev_dim = latent_dim
        for hidden_dim in reversed(hidden_dims):
            decoder_layers.extend([nn.Linear(prev_dim, hidden_dim), nn.BatchNorm1d(hidden_dim), nn.LeakyReLU(0.2)])
            prev_dim = hidden_dim
        self.decoder_base = nn.Sequential(*decoder_layers)
        self.disambiguation_layer = nn.Linear(prev_dim, self.unmerged_dim)
        self.bio_prior_loss = DisambiguationBiologicalPriorLoss(UNMERGED_MARKERS)

    def encode(self, x):
        return self.fc_mu(self.encoder(x)), self.fc_logvar(self.encoder(x))

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        return mu + torch.randn_like(std) * std

    def decode(self, z):
        x_disambiguated = self.disambiguation_layer(self.decoder_base(z))
        return self.remerge(x_disambiguated), x_disambiguated

    def remerge(self, x_disambiguated):
        recon_list = []
        unmerged_idx_map = {name: i for i, name in enumerate(UNMERGED_MARKERS)}
        for marker in ORIGINAL_MARKERS:
            cleaned_m = marker.split(' ')[0]
            if cleaned_m in MERGED_TO_UNMERGED:
                m1_name, m2_name = MERGED_TO_UNMERGED[cleaned_m]
                recon_list.append(x_disambiguated[:, unmerged_idx_map[m1_name]] + x_disambiguated[:, unmerged_idx_map[m2_name]])
            else:
                recon_list.append(x_disambiguated[:, unmerged_idx_map[cleaned_m]])
        return torch.stack(recon_list, dim=1)

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        return self.decode(z) + (mu, logvar)

# --- Helper Functions ---
def vae_loss_disambiguation(recon_x, x, x_disambiguated, mu, logvar, model, beta=1.0, gamma=1.0):
    recon_loss = nn.MSELoss()(recon_x, x)
    kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp()) / (x.size(0) * x.size(1))
    bio_loss = model.bio_prior_loss(x_disambiguated) if model.bio_prior_loss is not None and gamma > 0 else torch.tensor(0.0, device=recon_x.device)
    return recon_loss + beta * kl_loss + gamma * bio_loss, recon_loss, kl_loss, bio_loss

def train_vae_disambiguation(model, data_loader, optimizer, device, beta=1.0, gamma=1.0):
    model.train()
    total_loss, total_recon, total_kl, total_bio = 0, 0, 0, 0
    for batch_data, in data_loader:
        batch_data = batch_data.to(device)
        optimizer.zero_grad()
        recon_batch, x_disambiguated, mu, logvar = model(batch_data)
        loss, recon, kl, bio = vae_loss_disambiguation(recon_batch, batch_data, x_disambiguated, mu, logvar, model, beta, gamma)
        loss.backward()
        optimizer.step()
        total_loss += loss.item()
        total_recon += recon.item()
        total_kl += kl.item()
        total_bio += bio.item()
    n_batches = len(data_loader)
    return (total_loss/n_batches, total_recon/n_batches, total_kl/n_batches, total_bio/n_batches)

def compute_soft_archetype_assignments(latent_embeddings, archetype_latent, temperature=1.0):
    distances = cdist(latent_embeddings, archetype_latent, metric='euclidean')
    similarities = -distances / temperature
    soft_assignments = np.exp(similarities) / np.exp(similarities).sum(axis=1, keepdims=True)
    entropy = -np.sum(soft_assignments * np.log(soft_assignments + 1e-10), axis=1)
    return soft_assignments, entropy

def identify_archetypes(latent_embeddings, n_archetypes=10):
    from sklearn.cluster import KMeans
    kmeans = KMeans(n_clusters=n_archetypes, random_state=42, n_init=20).fit(latent_embeddings)
    return [np.argmin(np.linalg.norm(latent_embeddings - center, axis=1)) for center in kmeans.cluster_centers_]

def main():
    parser = argparse.ArgumentParser(description="Metacell Pipeline v2: VAE Model Training")
    parser.add_argument("--input", type=str, default="seacells_pipeline_output/integrated_metacells_17_markers.h5ad")
    parser.add_argument("--output-dir", type=str, default="vae_output")
    parser.add_argument("--n-archetypes", type=int, default=12)
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--latent-dim", type=int, default=10)
    parser.add_argument("--beta", type=float, default=0.01)
    parser.add_argument("--gamma", type=float, default=0.5)
    parser.add_argument("--beta-warmup", type=int, default=100)
    parser.add_argument("--gamma-cooldown", type=int, default=150)
    parser.add_argument("--gamma-min", type=float, default=0.0)
    parser.add_argument("--soft-assignment-temp", type=float, default=1.0)
    args = parser.parse_args()

    print("="*80 + f"\nVAE WITH DISAMBIGUATION LAYER (v6)\nUsing 17 available markers -> 23 un-merged markers\n" + "="*80)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.input)
    adata = adata[:, ORIGINAL_MARKERS].copy()

    X_train = adata.X.copy()

    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    train_loader = torch.utils.data.DataLoader(torch.utils.data.TensorDataset(torch.FloatTensor(X_train_scaled)), batch_size=64, shuffle=True)

    device = torch.device('mps' if torch.backends.mps.is_available() else 'cpu')
    model = VAEDisambiguation(input_dim=adata.n_vars, latent_dim=args.latent_dim).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=5e-4, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', patience=20)

    print(f"\nTraining VAE for {args.epochs} epochs on {device}...")
    for epoch in range(args.epochs):
        current_beta = args.beta * min(1.0, epoch / args.beta_warmup if args.beta_warmup > 0 else 1.0)
        if args.gamma_cooldown > 0:
            progress = min(1.0, epoch / args.gamma_cooldown)
            current_gamma = args.gamma * (1 - progress) + args.gamma_min * progress
        else:
            current_gamma = args.gamma

        loss, recon, kl, bio = train_vae_disambiguation(model, train_loader, optimizer, device, beta=current_beta, gamma=current_gamma)
        scheduler.step(loss)
        if (epoch + 1) % 50 == 0: print(f"  Epoch {epoch+1}/{args.epochs}: Loss={loss:.4f}, Recon={recon:.4f}, KL={kl:.4f}, Bio={bio:.4f}")

    model.eval()
    with torch.no_grad():
        mu_train, _ = model.encode(torch.FloatTensor(X_train_scaled).to(device))
        latent_train = mu_train.cpu().numpy()
    adata.obsm['X_vae'] = latent_train

    archetype_indices = identify_archetypes(latent_train, n_archetypes=args.n_archetypes)
    archetype_latent = latent_train[archetype_indices]

    soft_assign_train, entropy_train = compute_soft_archetype_assignments(latent_train, archetype_latent, temperature=args.soft_assignment_temp)
    for i in range(args.n_archetypes): adata.obs[f'archetype_{i}_prob'] = soft_assign_train[:, i]

    adata.obs['archetype_entropy'] = entropy_train
    train_labels = np.argmax(soft_assign_train, axis=1)
    adata.obs['archetype_id_hard'] = pd.Series(train_labels, index=adata.obs.index).astype('category')

    is_archetype_series = pd.Series('False', index=adata.obs.index)
    is_archetype_series.iloc[archetype_indices] = 'True'
    adata.obs['is_archetype'] = is_archetype_series.astype('category')

    sc.pp.neighbors(adata, use_rep='X_vae')
    sc.tl.umap(adata)

    torch.save({'model_state_dict': model.state_dict(), 'scaler': scaler, 'archetype_latent': archetype_latent}, output_dir / "vae_model.pt")
    adata.write_h5ad(output_dir / "metacells_with_archetypes.h5ad")

    with torch.no_grad():
        _, x_disambiguated_archetypes = model.decode(torch.FloatTensor(archetype_latent).to(device))
        pd.DataFrame(x_disambiguated_archetypes.cpu().numpy(), index=[f"Archetype_{i}" for i in range(args.n_archetypes)], columns=UNMERGED_MARKERS).to_csv(output_dir / "archetype_profiles_disambiguated.csv")

    print(f"\n✓ Analysis complete. Outputs saved to: {args.output_dir}")

if __name__ == "__main__":
    main()