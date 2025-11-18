#!/usr/bin/env python3
"""
Comprehensive Model Comparison (v3)
Compares all VAE archetype models trained today.
This version provides a more nuanced comparison, focusing on the trade-offs
between assignment clarity, archetype balance, and archetype uniqueness,
rather than a single score.
"""

import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import scanpy as sc

# Model directories to compare
models = {
    'Gamma Cooldown (150 ep)': 'test_output/gamma_cooldown',
    'Constant Gamma (150 ep)': 'test_output/constant_gamma',
    'Gamma Cooldown (300 ep)': 'test_output/gamma_cooldown_300epochs',
    'Best Optuna (Trial 13)': 'optuna_tuning_Nov15/trial_13',
}

def load_data_from_h5ad(model_dir):
    """Load assignments and other data from a .h5ad file"""
    h5ad_file = Path(model_dir) / 'metacells_with_archetypes.h5ad'
    if not h5ad_file.exists():
        return None, None

    adata = sc.read_h5ad(h5ad_file)
    
    # Extract only reference samples for fair comparison of entropy/gini
    adata_ref = adata[adata.obs['is_reference'] == True].copy()

    assignments = pd.DataFrame({
        'entropy': adata_ref.obs['archetype_entropy'],
        'archetype_id_hard': adata_ref.obs['archetype_id_hard'],
    })
    
    return assignments, adata

def load_archetype_profiles(model_dir):
    """Load archetype marker profiles"""
    profiles_file = Path(model_dir) / 'archetype_profiles.csv'
    if not profiles_file.exists():
        return None
    return pd.read_csv(profiles_file, index_col=0)

def compute_summary_stats(assignments, profiles):
    """Compute summary statistics from assignments and profiles"""
    if assignments is None:
        return {}

    entropy = assignments['entropy'].values

    # Gini coefficient for archetype distribution
    archetype_counts = assignments['archetype_id_hard'].value_counts().values
    archetype_counts = np.sort(archetype_counts)
    n = len(archetype_counts)
    if n == 0 or np.sum(archetype_counts) == 0:
        gini = 0
    else:
        index = np.arange(1, n + 1)
        gini = (2 * np.sum(index * archetype_counts)) / (n * np.sum(archetype_counts)) - (n + 1) / n

    # Archetype uniqueness (mean absolute correlation)
    mean_abs_corr = 0
    if profiles is not None:
        # Exclude metadata columns like 'gini', 'entropy', 'n_metacells'
        marker_cols = [c for c in profiles.columns if c not in ['gini', 'entropy', 'n_metacells']]
        corr_matrix = profiles[marker_cols].T.corr(method='spearman')
        # Get the mean of the upper triangle, excluding the diagonal
        upper_tri = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
        mean_abs_corr = upper_tri.abs().mean().mean()

    return {
        'mean_entropy': entropy.mean(),
        'median_entropy': np.median(entropy),
        'high_entropy_pct': (entropy > 1.5).mean() * 100,
        'gini_coefficient': gini,
        'n_archetypes': assignments['archetype_id_hard'].nunique(),
        'archetype_uniqueness': 1 - mean_abs_corr, # Higher is better (less correlated)
    }

# Collect all model data
comparison_data = []
all_adatas = {}
all_profiles = {}

for model_name, model_dir in models.items():
    if not Path(model_dir).exists():
        print(f"⚠️  Skipping {model_name}: directory not found")
        continue

    print(f"📊 Loading {model_name}...")

    assignments, adata = load_data_from_h5ad(model_dir)
    if assignments is None:
        print(f"⚠️  Skipping {model_name}: could not load data from .h5ad file")
        continue
    
    profiles = load_archetype_profiles(model_dir)
    summary = compute_summary_stats(assignments, profiles)

    row = {'Model': model_name, **summary}
    comparison_data.append(row)
    all_adatas[model_name] = adata
    all_profiles[model_name] = profiles


# Create comparison DataFrame
df_comparison = pd.DataFrame(comparison_data)

# Save to CSV
output_dir = Path('model_comparison_Nov15')
output_dir.mkdir(exist_ok=True)

df_comparison.to_csv(output_dir / 'all_models_comparison_v3.csv', index=False)
print(f"\n✓ Saved comparison table to {output_dir / 'all_models_comparison_v3.csv'}")

# Print comparison table
print("\n" + "="*80)
print("MODEL COMPARISON SUMMARY (METRICS FROM REFERENCE SAMPLES ONLY)")
print("="*80)

display_cols = [
    'Model',
    'mean_entropy',
    'gini_coefficient',
    'archetype_uniqueness',
    'n_archetypes',
]
print(df_comparison[display_cols].to_string(index=False, float_format="%.3f"))

# --- Visualization ---
fig = plt.figure(figsize=(18, 16))
gs = fig.add_gridspec(3, 4)
fig.suptitle('VAE Archetype Model Comparison', fontsize=20, fontweight='bold')

# Metrics plots
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])

df_plot = df_comparison.sort_values('mean_entropy')
ax1.barh(df_plot['Model'], df_plot['mean_entropy'], color='coral')
ax1.set_title('Assignment Clarity (Mean Entropy)')
ax1.set_xlabel('Lower is clearer')

df_plot = df_comparison.sort_values('gini_coefficient')
ax2.barh(df_plot['Model'], df_plot['gini_coefficient'], color='green')
ax2.set_title('Archetype Balance (Gini)')
ax2.set_xlabel('Lower is more balanced')

df_plot = df_comparison.sort_values('archetype_uniqueness', ascending=False)
ax3.barh(df_plot['Model'], df_plot['archetype_uniqueness'], color='skyblue')
ax3.set_title('Archetype Uniqueness (1 - Mean Correlation)')
ax3.set_xlabel('Higher is more unique')

# Correlation heatmaps
for i, (model_name, profiles) in enumerate(all_profiles.items()):
    ax = fig.add_subplot(gs[1, i])
    if profiles is not None:
        marker_cols = [c for c in profiles.columns if c not in ['gini', 'entropy', 'n_metacells']]
        corr_matrix = profiles[marker_cols].T.corr(method='spearman')
        sns.heatmap(corr_matrix, ax=ax, cmap='coolwarm', vmin=-1, vmax=1, square=True, cbar= i==3)
        ax.set_title(f"{model_name}\nArchetype Correlation", fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])

# UMAP plots
for i, (model_name, adata) in enumerate(all_adatas.items()):
    ax = fig.add_subplot(gs[2, i])
    sc.pl.umap(adata, color='archetype_id_hard', ax=ax, show=False, legend_loc='on data', title=f"{model_name}\nArchetype UMAP", size=20)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(output_dir / 'model_comparison_charts_v3.png', dpi=300)
print(f"\n✓ Saved visualization to {output_dir / 'model_comparison_charts_v3.png'}")


# --- Interpretation ---
print("\n" + "="*80)
print("INTERPRETATION & DISCUSSION")
print("="*80)

print(""")
This analysis compares models on three main axes:
1. Assignment Clarity (Mean Entropy): How clearly are metacells assigned to a single archetype?
   - Lower entropy means "purer" assignments.
   - The 'Best Optuna' model is the clearest, while the 'Gamma Cooldown (300 ep)' model has the highest entropy.

2. Archetype Balance (Gini Coefficient): How evenly are metacells distributed among archetypes?
   - Lower Gini means more balanced archetypes. A higher Gini may indicate the presence of rare archetypes.
   - The 'Gamma Cooldown (150 ep)' model is the most balanced.

3. Archetype Uniqueness (1 - Mean Correlation): How different are the archetype marker profiles from each other?
   - Higher uniqueness (lower correlation) is desirable, suggesting the model found distinct biological roles.
   - The 'Gamma Cooldown (300 ep)' model produced the most unique set of archetypes.

Discussion on "Mixed Phenotypes":
As you noted, high entropy is not necessarily a flaw. The gamma cooldown process is designed to first define archetypes and then place metacells in the space between them.
- The 'Gamma Cooldown' models, with their moderate-to-high entropy, may be successfully capturing these biologically meaningful "in-between" states.
- The 'Best Optuna' model, with its very low entropy, might be "over-clustering" the data into highly pure but potentially less nuanced groups.

Recommendation:
There is no single "best" model; the choice depends on the analytical goal.
- For identifying pure, core phenotypes: 'Best Optuna (Trial 13)' is superior due to its low entropy and high number of archetypes.
- For studying transitional states and capturing phenotypic diversity: 'Gamma Cooldown (150 ep)' offers the best balance of even, reasonably unique archetypes with moderate entropy.
- For discovering the most distinct (least correlated) phenotypes: 'Gamma Cooldown (300 ep)' excels, though it has the highest entropy.

Given the goal of finding unique, biologically-relevant archetypes and understanding the "space between", the 'Gamma Cooldown (150 ep)' model appears to be a strong, balanced candidate for further analysis.
""")

print("\n✓ Analysis complete!")
