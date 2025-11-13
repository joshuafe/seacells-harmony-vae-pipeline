import argparse
import json
import os
import time
import optuna
import joblib
import glob
from dataclasses import dataclass, asdict
from pathlib import Path
import pandas as pd
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset
from tqdm.auto import tqdm
import scanpy as sc
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply
from dotenv import load_dotenv

load_dotenv(dotenv_path=Path(__file__).parent.parent / '.env')

def set_seed(seed: int):
    """Sets the random seed for reproducibility."""
    import random
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    print(f"Random seed set to {seed} for reproducibility.")

class HeatDiffusionFast:
    """Heat diffusion using PyNNDescent."""
    def __init__(self, k=15, t=0.5, use_pynndescent=True):
        self.k = k
        self.t = t
        self.use_pynndescent = use_pynndescent
        self.L = None
        
    def fit(self, X):
        if self.k == 0:
            print("k=0: Skipping heat diffusion")
            return self
            
        t0 = time.time()
        print(f"\n🔥 Building {self.k}-NN graph for {X.shape[0]:,} data points...")
        
        if self.use_pynndescent:
            from pynndescent import NNDescent
            from scipy.sparse import csr_matrix
            
            index = NNDescent(X, n_neighbors=self.k+1, metric='euclidean', n_jobs=-1)
            indices, distances = index.neighbor_graph
            indices, distances = indices[:, 1:], distances[:, 1:]
            
            n = X.shape[0]
            row_ind = np.repeat(np.arange(n), self.k)
            col_ind = indices.flatten()
            data = np.ones_like(col_ind) # unweighted graph
            
            A = csr_matrix((data, (row_ind, col_ind)), shape=(n, n))
        else:
            from sklearn.neighbors import kneighbors_graph
            A = kneighbors_graph(X, self.k, mode='connectivity', n_jobs=-1)
        
        t1 = time.time()
        print(f"   → k-NN graph built in {t1-t0:.1f}s")
        
        W = (A + A.T) / 2
        
        print(f"   → Computing Laplacian...")
        t2 = time.time()
        
        D = np.asarray(W.sum(axis=1)).flatten()
        D_inv_sqrt = diags(1.0 / np.sqrt(D + 1e-10))
        I = diags(np.ones(X.shape[0]))
        self.L = I - D_inv_sqrt @ W @ D_inv_sqrt
        
        t3 = time.time()
        print(f"   → Laplacian computed in {t3-t2:.1f}s")
        return self
    
    def transform(self, X):
        if self.k == 0 or self.L is None:
            return X
        
        t0 = time.time()
        print(f"   → Applying heat diffusion (t={self.t})...")
        X_smooth = expm_multiply(-self.t * self.L, X)
        t1 = time.time()
        print(f"✅ Diffusion applied in {t1-t0:.1f}s")
        return X_smooth
    
    def fit_transform(self, X):
        return self.fit(X).transform(X)

def _approx_split(composite: np.ndarray, bias: float = 0.5) -> Tuple[np.ndarray, np.ndarray]:
    """Approximates splitting of composite channels."""
    bias = np.clip(bias, 0.05, 0.95)
    return composite * bias, composite * (1.0 - bias)

def get_marker_data(adata: sc.AnnData, marker_name: str) -> np.ndarray:
    """Safely gets marker data from an AnnData object, handling composite markers."""
    if marker_name in adata.var_names:
        return adata[:, marker_name].X.flatten()
    
    # Handle composite markers
    if marker_name == 'CD3':
        return _approx_split(adata[:, 'CD3+CD19'].X.flatten(), 0.55)[0]
    if marker_name == 'CD19':
        return _approx_split(adata[:, 'CD3+CD19'].X.flatten(), 0.55)[1]
    if marker_name == 'CD8':
        return _approx_split(adata[:, 'CD8+CD14'].X.flatten(), 0.4)[0]
    if marker_name == 'CD4':
        return _approx_split(adata[:, 'CD4+CD33'].X.flatten(), 0.35)[0]
    
    raise ValueError(f"Marker {marker_name} not found in adata.var_names or composite definitions.")

def compute_metacell_aux_labels(adata: sc.AnnData) -> Dict[str, np.ndarray]:
    """Computes auxiliary labels for metacells from the AnnData object."""
    cd4_data = get_marker_data(adata, 'CD4')
    cd8_data = get_marker_data(adata, 'CD8')
    cd69_data = get_marker_data(adata, 'CD69')
    
    return {
        'cd4': (cd4_data > np.percentile(cd4_data, 60)).astype(np.float32).reshape(-1, 1),
        'cd8': (cd8_data > np.percentile(cd8_data, 60)).astype(np.float32).reshape(-1, 1),
        'activation': (cd69_data > np.percentile(cd69_data, 70)).astype(np.float32).reshape(-1, 1)
    }

class VAEDataset(Dataset):
    def __init__(self, data: np.ndarray, aux_targets: Optional[np.ndarray] = None):
        self.data = torch.tensor(data, dtype=torch.float32)
        self.aux = torch.tensor(aux_targets, dtype=torch.float32) if aux_targets is not None else None
    def __len__(self): return len(self.data)
    def __getitem__(self, idx):
        return (self.data[idx], self.aux[idx]) if self.aux is not None else self.data[idx]

class ArchetypeLayer(nn.Module):
    def __init__(self, input_dim: int, n_archetypes: int, archetype_dim: int, tau: float = 0.5):
        super().__init__()
        self.encoder = nn.Linear(input_dim, n_archetypes)
        self.tau = tau
        self.register_parameter("archetypes", nn.Parameter(torch.randn(n_archetypes, archetype_dim) * 0.1))

    def forward(self, h: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        logits = self.encoder(h)
        logits = logits - logits.max(dim=-1, keepdim=True).values
        alpha = F.softmax(logits / self.tau, dim=-1)
        archetypes = F.normalize(self.archetypes, p=2, dim=1)
        return alpha @ archetypes, alpha

class ArchetypeVAE(nn.Module):
    def __init__(self, input_dim: int, latent_dim: int, hidden_dim: int, n_archetypes: int, archetype_dim: int, tau: float):
        super().__init__()
        self.encoder = nn.Sequential(nn.Linear(input_dim, hidden_dim), nn.BatchNorm1d(hidden_dim), nn.LeakyReLU(0.2))
        self.fc_mu, self.fc_logvar = nn.Linear(hidden_dim, latent_dim), nn.Linear(hidden_dim, latent_dim)
        self.archetype_layer = ArchetypeLayer(latent_dim, n_archetypes, archetype_dim, tau=tau)
        decoder_dim = max(hidden_dim // 2, archetype_dim)
        self.decoder = nn.Sequential(nn.Linear(archetype_dim, decoder_dim), nn.LeakyReLU(0.2), nn.Linear(decoder_dim, input_dim))
        self.aux_heads = nn.ModuleDict({'cd4': nn.Linear(archetype_dim, 1), 'cd8': nn.Linear(archetype_dim, 1), 'activation': nn.Linear(archetype_dim, 1)})

    def encode(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        h = self.encoder(x)
        return self.fc_mu(h), self.fc_logvar(h).clamp(-6.0, 6.0), h

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        return mu + torch.randn_like(std) * std

    def forward(self, x: torch.Tensor):
        mu, logvar, h = self.encode(x)
        z = self.reparameterize(mu, logvar)
        archetype_latent, alpha = self.archetype_layer(z)
        recon = self.decoder(archetype_latent)
        aux_logits = {name: head(archetype_latent) for name, head in self.aux_heads.items()}
        return recon, mu, logvar, alpha, aux_logits

@dataclass
class VAEParams:
    latent_dim: int = 8
    hidden_dim: int = 64
    n_archetypes: int = 12
    archetype_dim: int = 12
    tau_init: float = 0.8
    tau_final: float = 0.2
    tau_decay_epochs: int = 150
    lr: float = 1e-3
    kl_target: float = 0.01
    kl_warmup_epochs: int = 50
    entropy_weight: float = 0.02
    diversity_weight: float = 0.2
    uniform_weight: float = 0.05
    aux_weight: float = 0.3
    epochs: int = 200
    batch_size: int = 64
    early_stop_patience: int = 30
    early_stop_delta: float = 1e-5
    use_heat_diffusion: bool = False
    heat_diffusion_k: int = 5
    heat_diffusion_t: float = 1.0

def _kl_gaussian(mu, logvar):
    return -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=1)

def archetype_diversity_loss(archetypes: torch.Tensor) -> torch.Tensor:
    norm = F.normalize(archetypes, p=2, dim=1)
    cosine_sim = norm @ norm.t()
    return (cosine_sim - torch.eye(norm.size(0), device=norm.device)).pow(2).mean()

def train_vae(model: ArchetypeVAE, data: np.ndarray, params: VAEParams, aux_targets: Dict[str, np.ndarray], device: str = "cpu") -> Dict:
    aux_matrix = np.concatenate([aux_targets[name] for name in model.aux_heads.keys()], axis=1)
    # Handle small batch sizes
    batch_size = min(params.batch_size, data.shape[0] // 2)
    if batch_size <= 0: batch_size = data.shape[0]
    
    loader = DataLoader(VAEDataset(data, aux_matrix), batch_size=batch_size, shuffle=True, drop_last=True)
    model.to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=params.lr)
    
    best_loss, wait = float('inf'), 0

    for epoch in tqdm(range(params.epochs), desc="VAE Training", leave=False):
        model.train()
        model.archetype_layer.tau = max(params.tau_init + (params.tau_final - params.tau_init) * (epoch / max(1, params.tau_decay_epochs)), params.tau_final)
        kl_weight = params.kl_target * min(1.0, epoch / max(1, params.kl_warmup_epochs))
        
        epoch_loss = 0.0
        for batch_data, aux_batch in loader:
            batch = batch_data.to(device)
            aux_batch = aux_batch.to(device)
            
            recon, mu, logvar, alpha, aux_logits = model(batch)
            
            recon_loss = F.mse_loss(recon, batch)
            kl_loss = _kl_gaussian(mu, logvar).mean()
            entropy = -(alpha * torch.log(alpha + 1e-10)).sum(dim=1).mean()
            div_loss = archetype_diversity_loss(model.archetype_layer.archetypes)
            uniform_loss = F.mse_loss(alpha.mean(dim=0), torch.full((alpha.size(1),), 1.0 / alpha.size(1), device=alpha.device))
            aux_loss = sum(F.binary_cross_entropy_with_logits(logits, aux_batch[:, i:i+1]) for i, logits in enumerate(aux_logits.values()))
            
            loss = (recon_loss + 
                    kl_weight * kl_loss + 
                    params.entropy_weight * entropy + 
                    params.diversity_weight * div_loss + 
                    params.uniform_weight * uniform_loss + 
                    params.aux_weight * aux_loss)
            
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()

        avg_epoch_loss = epoch_loss / len(loader)
        if (epoch + 1) % 10 == 0:
            print(f"  Epoch {epoch+1}: loss={avg_epoch_loss:.4f}")

        if avg_epoch_loss < best_loss - params.early_stop_delta:
            best_loss, wait = avg_epoch_loss, 0
        else:
            wait += 1
        if params.early_stop_patience > 0 and wait >= params.early_stop_patience:
            print(f"Early stopping at epoch {epoch+1}.")
            break
            
    model.eval()
    with torch.no_grad():
        full_data_tensor = torch.tensor(data, dtype=torch.float32).to(device)
        recon, mu, logvar, alpha, _ = model(full_data_tensor)
        metrics = {
            'final_recon': F.mse_loss(recon, full_data_tensor).item(),
            'final_kl': _kl_gaussian(mu, logvar).mean().item(),
            'mean_purity': alpha.max(dim=1).values.mean().item(),
            'archetype_harmony_profiles': model.decoder(F.normalize(model.archetype_layer.archetypes, p=2, dim=1)).cpu().numpy(),
            'alpha': alpha.cpu().numpy(),
        }
    return metrics

def get_expression_profiles(adata: sc.AnnData, harmony_profiles: np.ndarray) -> np.ndarray:
    """Maps low-dimensional harmony profiles back to high-dimensional expression space."""
    pca_comps = adata.varm['PCs']
    n_harmony_dims = harmony_profiles.shape[1]
    
    if n_harmony_dims > pca_comps.shape[1]:
        raise ValueError(f"Harmony profiles have more dimensions ({n_harmony_dims}) than available PCA components ({pca_comps.shape[1]})")
    
    scaled_expression = harmony_profiles @ pca_comps[:, :n_harmony_dims].T
    
    mean = adata.var['mean'].to_numpy()
    std = adata.var['std'].to_numpy()
    
    expression_profiles = scaled_expression * (std + 1e-8) + mean
    return expression_profiles

def save_visual_summary(metrics: Dict, adata: sc.AnnData, output_dir: Path):
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    output_dir.mkdir(parents=True, exist_ok=True)

    harmony_profiles = metrics['archetype_harmony_profiles']
    expression_profiles = get_expression_profiles(adata, harmony_profiles)
    channel_names = adata.var_names.tolist()
    alpha = metrics['alpha']
    
    archetype_names = [f"Archetype {i+1}" for i in range(expression_profiles.shape[0])]

    plt.figure(figsize=(14, 8))
    sns.heatmap(expression_profiles, yticklabels=archetype_names, xticklabels=channel_names, cmap='viridis', center=0)
    plt.title('Archetype Marker Expression Profiles (Reconstructed)')
    plt.tight_layout()
    plt.savefig(output_dir / 'archetype_profiles_heatmap.png')
    plt.close()

    embedding = adata.obsm['X_umap']
    cell_archetype_labels = np.argmax(alpha, axis=1)

    fig, ax = plt.subplots(figsize=(12, 10))
    scatter = ax.scatter(embedding[:, 0], embedding[:, 1], c=cell_archetype_labels, cmap='tab20', s=15, alpha=0.7)
    ax.set_title('UMAP of Metacells Colored by Archetype Assignment')
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    legend_elements = scatter.legend_elements(prop='colors', num=len(archetype_names))
    ax.legend(legend_elements[0], archetype_names, bbox_to_anchor=(1.02, 1), loc='upper left', title="Archetypes")
    fig.tight_layout()
    plt.savefig(output_dir / 'metacell_distribution_umap.png')
    plt.close(fig)

    archetype_usage = np.mean(alpha, axis=0)
    usage_df = pd.DataFrame({'archetype': archetype_names, 'proportion': archetype_usage})
    usage_df = usage_df.sort_values('proportion', ascending=False)

    plt.figure(figsize=(10, 8))
    sns.stripplot(data=usage_df, y='archetype', x='proportion', orient='h', size=15, order=usage_df['archetype'])
    plt.title('Archetype Proportions (Usage on Metacells)')
    plt.xlabel('Proportion')
    plt.ylabel('')
    plt.grid(axis='x', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_dir / 'archetype_usage_dotplot.png')
    plt.close()
    
    print(f"Visual summary saved to {output_dir}")

@dataclass
class ExperimentConfig:
    vae: VAEParams

def run_experiment(adata: sc.AnnData, config: ExperimentConfig, output_dir: Path, device: str = "cpu") -> Dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    
    vae_data = adata.obsm['X_harmony'].copy()
    params = config.vae

    if params.use_heat_diffusion:
        print(f"Applying Heat Diffusion (k={params.heat_diffusion_k}, t={params.heat_diffusion_t})")
        diffuser = HeatDiffusionFast(k=params.heat_diffusion_k, t=params.heat_diffusion_t)
        vae_data = diffuser.fit_transform(vae_data)

    aux_targets = compute_metacell_aux_labels(adata)

    model = ArchetypeVAE(
        input_dim=vae_data.shape[1],
        latent_dim=params.latent_dim,
        hidden_dim=params.hidden_dim,
        n_archetypes=params.n_archetypes,
        archetype_dim=params.archetype_dim,
        tau=params.tau_init
    )

    metrics = train_vae(model, vae_data, params, aux_targets=aux_targets, device=device)
    
    torch.save(model.state_dict(), output_dir / 'metacell_vae_model.pt')
    metrics['input_dimensions'] = vae_data.shape[1]
    
    return metrics

def _to_serialisable(obj):
    if isinstance(obj, (np.ndarray, torch.Tensor)): return obj.tolist()
    if isinstance(obj, dict): return {k: _to_serialisable(v) for k, v in obj.items()}
    if isinstance(obj, list): return [_to_serialisable(v) for v in obj]
    return obj

def objective(trial: optuna.trial.Trial, adata: sc.AnnData, device: str, output_dir: Path):
    trial_dir = output_dir / f"trial_{trial.number:03d}"

    try:
        params = VAEParams(
            lr=trial.suggest_float("lr", 1e-4, 5e-3, log=True),
            latent_dim=trial.suggest_categorical("latent_dim", [8, 12, 16]),
            hidden_dim=trial.suggest_categorical("hidden_dim", [32, 64, 96]),
            n_archetypes=trial.suggest_int("n_archetypes", 6, 20),
            archetype_dim=trial.suggest_categorical("archetype_dim", [8, 12, 16]),
            kl_target=trial.suggest_float("kl_target", 1e-3, 5e-2, log=True),
            diversity_weight=trial.suggest_float("diversity_weight", 0.1, 1.0),
            batch_size=trial.suggest_categorical("batch_size", [32, 64]),
            epochs=250,
            use_heat_diffusion=trial.suggest_categorical("use_heat_diffusion", [True, False]),
            heat_diffusion_k=trial.suggest_int("heat_diffusion_k", 3, 10) if trial.params['use_heat_diffusion'] else 0,
            heat_diffusion_t=trial.suggest_float("heat_diffusion_t", 0.5, 2.0) if trial.params['use_heat_diffusion'] else 0,
        )
        config = ExperimentConfig(vae=params)

        trial_dir.mkdir(parents=True, exist_ok=True)
        with open(trial_dir / 'config.json', 'w') as f:
            json.dump(asdict(config), f, indent=2, default=_to_serialisable)

        metrics = run_experiment(adata, config, trial_dir, device=device)

        metrics_to_save = {k: _to_serialisable(v) for k, v in metrics.items() if k not in ['archetype_harmony_profiles', 'alpha']}
        with open(trial_dir / 'metrics.json', 'w') as f:
            json.dump(metrics_to_save, f, indent=2)

        np.save(trial_dir / 'alpha.npy', metrics['alpha'])
        np.save(trial_dir / 'archetype_harmony_profiles.npy', metrics['archetype_harmony_profiles'])

        save_visual_summary(metrics, adata, trial_dir)

        return metrics.get('mean_purity', 0.0), metrics.get('final_recon', float('inf'))

    except Exception as e:
        import traceback
        print(f"Trial {trial.number} FAILED: {e}")
        traceback.print_exc()
        return 0.0, float('inf')

def main():
    parser = argparse.ArgumentParser(description="Archetype VAE on Harmony-corrected Metacells")
    parser.add_argument("--mode", default="optuna", choices=["single", "optuna"], help="Execution mode.")
    parser.add_argument("--input-h5ad", type=Path, default="seacells_output/batch_corrected/integrated_metacells_harmony.h5ad", help="Path to the integrated h5ad file.")
    parser.add_argument("--output", type=Path, default="./metacell_vae_output", help="Directory for outputs.")
    parser.add_argument("--n-trials", type=int, default=50, help="Number of Optuna trials.")
    parser.add_argument("--device", default="cpu", help="Device to use for training (e.g., 'cpu', 'cuda').")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility.")
    args = parser.parse_args()

    set_seed(args.seed)
    
    if not args.input_h5ad.exists():
        raise FileNotFoundError(f"Input file not found: {args.input_h5ad}. Please run the seacells pipeline first.")

    print("Loading integrated metacell data...")
    adata = sc.read_h5ad(args.input_h5ad)
    
    required_keys = ['X_harmony', 'PCs', 'mean', 'std', 'X_umap']
    for key in required_keys:
        if (key in adata.obsm) or (key in adata.varm) or (key in adata.var):
            continue
        raise ValueError(f"Required key '{key}' not found in the AnnData object. Please ensure the preprocessing pipeline was run completely.")

    if args.mode == 'optuna':
        study_name = "metacell_archetype_optimization"
        args.output.mkdir(parents=True, exist_ok=True)
        storage_path = f"sqlite:///{args.output / 'optuna_metacell_study.db'}"
        study = optuna.create_study(
            study_name=study_name,
            storage=storage_path,
            directions=["maximize", "minimize"], # Maximize purity, Minimize reconstruction error
            load_if_exists=True
        )

        objective_func = lambda trial: objective(trial, adata, args.device, args.output)
        study.optimize(objective_func, n_trials=args.n_trials)
        
        print("\nOptuna sweep complete.")
        print(f"Best trials: {study.best_trials}")

    elif args.mode == 'single':
        print("Running a single experiment with default parameters...")
        params = VAEParams(use_heat_diffusion=True) # Example: turn on heat diffusion
        config = ExperimentConfig(vae=params)
        output_dir = args.output / "single_run"
        
        with open(output_dir / 'config.json', 'w') as f:
            json.dump(asdict(config), f, indent=2, default=_to_serialisable)
            
        metrics = run_experiment(adata, config, output_dir, device=args.device)
        save_visual_summary(metrics, adata, output_dir)
        print(f"Single run finished. Results in {output_dir}")

if __name__ == "__main__":
    main()
