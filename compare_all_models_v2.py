#!/usr/bin/env python3
"""
Comprehensive Model Comparison (v2)
Compares all VAE archetype models trained today.
This version is adapted to work with the output files that are actually generated,
reading data from .h5ad files instead of missing .csv and .json files.
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

def load_assignments_from_h5ad(model_dir):
    """Load archetype assignments from a .h5ad file"""
    h5ad_file = Path(model_dir) / 'metacells_with_archetypes.h5ad'
    if not h5ad_file.exists():
        # Disable ignores to check if the file is just being ignored
        if not Path(h5ad_file).is_file():
             return None

    adata = sc.read_h5ad(h5ad_file)
    
    # Extract only reference samples for fair comparison of entropy/gini
    adata_ref = adata[adata.obs['is_reference'] == True]

    df = pd.DataFrame({
        'entropy': adata_ref.obs['archetype_entropy'],
        'archetype_id_hard': adata_ref.obs['archetype_id_hard'],
    })
    return df

def compute_summary_stats(assignments):
    """Compute summary statistics from assignments"""
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

    return {
        'mean_entropy': entropy.mean(),
        'median_entropy': np.median(entropy),
        'high_entropy_pct': (entropy > 1.5).mean() * 100,
        'gini_coefficient': gini,
        'n_archetypes': assignments['archetype_id_hard'].nunique(),
    }

# Collect all model data
comparison_data = []

for model_name, model_dir in models.items():
    if not Path(model_dir).exists():
        print(f"⚠️  Skipping {model_name}: directory not found")
        continue

    print(f"📊 Loading {model_name}...")

    assignments = load_assignments_from_h5ad(model_dir)
    if assignments is None:
        print(f"⚠️  Skipping {model_name}: could not load assignments from .h5ad file")
        continue
        
    summary = compute_summary_stats(assignments)

    row = {
        'Model': model_name,
        **summary,
    }

    comparison_data.append(row)

# Create comparison DataFrame
df_comparison = pd.DataFrame(comparison_data)

# Save to CSV
output_dir = Path('model_comparison_Nov15')
output_dir.mkdir(exist_ok=True)

df_comparison.to_csv(output_dir / 'all_models_comparison_v2.csv', index=False)
print(f"\n✓ Saved comparison table to {output_dir / 'all_models_comparison_v2.csv'}")

# Print comparison table
print("\n" + "="*80)
print("MODEL COMPARISON SUMMARY (METRICS FROM REFERENCE SAMPLES ONLY)")
print("="*80)

# Select key columns for display
display_cols = [
    'Model',
    'mean_entropy',
    'high_entropy_pct',
    'gini_coefficient',
    'n_archetypes',
]

print(df_comparison[display_cols].to_string(index=False))

print("\n" + "="*80)
print("INTERPRETATION")
print("="*80)

# Find best models by different criteria
best_entropy = df_comparison.loc[df_comparison['mean_entropy'].idxmin()]
best_gini = df_comparison.loc[df_comparison['gini_coefficient'].idxmin()]

print(f"\n🎯 Clearest Assignments: {best_entropy['Model']} (entropy = {best_entropy['mean_entropy']:.4f})")
print(f"   - This model produces the most distinct archetypes on the reference dataset.")
print(f"   - {best_entropy['high_entropy_pct']:.1f}% high-entropy metacells in reference set.")

print(f"\n📊 Most Balanced: {best_gini['Model']} (Gini = {best_gini['gini_coefficient']:.3f})")
print(f"   - This model distributes reference metacells most evenly across archetypes.")
print(f"   - Avoids archetype collapse or duplication.")

# Create visualization
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('VAE Archetype Model Comparison (Assignment-Based Metrics)', fontsize=16, fontweight='bold')

# 1. Entropy comparison
ax = axes[0, 0]
df_plot = df_comparison.sort_values('mean_entropy')
ax.barh(df_plot['Model'], df_plot['mean_entropy'], color='coral')
ax.set_xlabel('Mean Entropy (lower = clearer assignments)')
ax.set_title('Assignment Clarity (Reference Samples)')
ax.invert_yaxis()
for i, (model, ent) in enumerate(zip(df_plot['Model'], df_plot['mean_entropy'])):
    ax.text(ent + 0.03, i, f'{ent:.3f}', va='center')

# 2. Gini coefficient
ax = axes[0, 1]
df_plot = df_comparison.sort_values('gini_coefficient')
ax.barh(df_plot['Model'], df_plot['gini_coefficient'], color='green')
ax.set_xlabel('Gini Coefficient (lower = more balanced)')
ax.set_title('Archetype Distribution Balance (Reference Samples)')
ax.invert_yaxis()
for i, (model, gini) in enumerate(zip(df_plot['Model'], df_plot['gini_coefficient'])):
    ax.text(gini + 0.005, i, f'{gini:.3f}', va='center')

# 3. High-entropy percentage
ax = axes[1, 0]
df_plot = df_comparison.sort_values('high_entropy_pct', ascending=False)
ax.barh(df_plot['Model'], df_plot['high_entropy_pct'], color='purple')
ax.set_xlabel('% High-Entropy Metacells (>1.5)')
ax.set_title('Mixed Phenotype Prevalence (Reference Samples)')
ax.invert_yaxis()
for i, (model, pct) in enumerate(zip(df_plot['Model'], df_plot['high_entropy_pct'])):
    ax.text(pct + 1, i, f'{pct:.1f}%', va='center')
    
# 4. Number of archetypes
ax = axes[1, 1]
df_plot = df_comparison.sort_values('n_archetypes')
ax.barh(df_plot['Model'], df_plot['n_archetypes'], color='orange')
ax.set_xlabel('Number of unique archetypes')
ax.set_title('Number of Archetypes (Reference Samples)')
ax.invert_yaxis()
ax.set_xlim(0, df_plot['n_archetypes'].max() + 1)
for i, (model, n) in enumerate(zip(df_plot['Model'], df_plot['n_archetypes'])):
    ax.text(n + 0.1, i, f'{n}', va='center')


plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(output_dir / 'model_comparison_charts_v2.png', dpi=300, bbox_inches='tight')
print(f"\n✓ Saved visualization to {output_dir / 'model_comparison_charts_v2.png'}")

# Recommendation
print("\n" + "="*80)
print("RECOMMENDATION")
print("="*80)

# Score each model (lower is better)
# Removed loss from score, re-weighting entropy and gini
df_comparison['score'] = (
    df_comparison['mean_entropy'] / df_comparison['mean_entropy'].min() * 0.6 +  # 60% weight
    df_comparison['gini_coefficient'] / df_comparison['gini_coefficient'].min() * 0.4  # 40% weight
)

best_overall = df_comparison.loc[df_comparison['score'].idxmin()]

print(f"\n🌟 RECOMMENDED MODEL: {best_overall['Model']}")
print(f"\n   Based on a composite score of assignment clarity (entropy) and archetype balance (Gini).")
print(f"   Composite score: {best_overall['score']:.3f} (lower is better)")
print(f"\n   Key metrics (from reference samples):")
print(f"   - Mean entropy: {best_overall['mean_entropy']:.4f}")
print(f"   - Gini coefficient: {best_overall['gini_coefficient']:.3f}")
print(f"   - High-entropy metacells: {best_overall['high_entropy_pct']:.1f}%")
print(f"   - Number of archetypes: {best_overall['n_archetypes']:.0f}")

print("\n" + "="*80)
print("NEXT STEPS")
print("="*80)
print(f"\n1. Use the recommended model ('{best_overall['Model']}') for production analysis.")
print(f"2. Apply to full dataset for sample composition analysis using scripts like 'analyze_sample_composition.py'.")
print(f"3. Proceed with external validation to address batch effects (see VALIDATION_ROADMAP.md).")

print("\n✓ Analysis complete!")
