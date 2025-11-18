#!/usr/bin/env python3
"""
Comprehensive Model Comparison
Compares all VAE archetype models trained today
"""

import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Model directories to compare
models = {
    'Gamma Cooldown (150 ep)': 'test_output/gamma_cooldown',
    'Constant Gamma (150 ep)': 'test_output/constant_gamma',
    'Gamma Cooldown (300 ep)': 'test_output/gamma_cooldown_300epochs',
    'Best Optuna (Trial 13)': 'optuna_tuning_Nov15/trial_13',
}

def load_metrics(model_dir):
    """Load training metrics from a model directory"""
    metrics_file = Path(model_dir) / 'training_metrics.csv'
    if not metrics_file.exists():
        return None

    df = pd.read_csv(metrics_file)
    final = df.iloc[-1]

    return {
        'total_loss': final['total_loss'],
        'recon_loss': final['recon_loss'],
        'kl_div': final['kl_div'],
        'prior_loss': final.get('prior_loss', 0.0),
        'final_epoch': final['epoch'],
    }

def load_config(model_dir):
    """Load model configuration"""
    config_file = Path(model_dir) / 'config.json'
    if not config_file.exists():
        return {}

    with open(config_file, 'r') as f:
        return json.load(f)

def load_assignments(model_dir):
    """Load archetype assignments"""
    assign_file = Path(model_dir) / 'archetype_assignments.csv'
    if not assign_file.exists():
        return None

    return pd.read_csv(assign_file)

def compute_summary_stats(assignments):
    """Compute summary statistics from assignments"""
    if assignments is None:
        return {}

    entropy = assignments['entropy'].values

    # Gini coefficient for archetype distribution
    archetype_counts = assignments['archetype_id_hard'].value_counts().values
    archetype_counts = np.sort(archetype_counts)
    n = len(archetype_counts)
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

    metrics = load_metrics(model_dir)
    config = load_config(model_dir)
    assignments = load_assignments(model_dir)
    summary = compute_summary_stats(assignments)

    row = {
        'Model': model_name,
        **metrics,
        **summary,
        **{f'config_{k}': v for k, v in config.items() if k in ['n_archetypes', 'latent_dim', 'beta', 'gamma', 'gamma_cooldown', 'gamma_min', 'soft_temp']},
    }

    comparison_data.append(row)

# Create comparison DataFrame
df_comparison = pd.DataFrame(comparison_data)

# Save to CSV
output_dir = Path('model_comparison_Nov15')
output_dir.mkdir(exist_ok=True)

df_comparison.to_csv(output_dir / 'all_models_comparison.csv', index=False)
print(f"\n✓ Saved comparison table to {output_dir / 'all_models_comparison.csv'}")

# Print comparison table
print("\n" + "="*80)
print("MODEL COMPARISON SUMMARY")
print("="*80)

# Select key columns for display
display_cols = [
    'Model',
    'total_loss',
    'mean_entropy',
    'high_entropy_pct',
    'gini_coefficient',
    'n_archetypes',
    'config_gamma',
    'config_gamma_cooldown',
]

print(df_comparison[display_cols].to_string(index=False))

print("\n" + "="*80)
print("INTERPRETATION")
print("="*80)

# Find best models by different criteria
best_loss = df_comparison.loc[df_comparison['total_loss'].idxmin()]
best_entropy = df_comparison.loc[df_comparison['mean_entropy'].idxmin()]
best_gini = df_comparison.loc[df_comparison['gini_coefficient'].idxmin()]

print(f"\n🏆 Lowest Loss: {best_loss['Model']} (loss = {best_loss['total_loss']:.4f})")
print(f"   - This model fits the training data best")
print(f"   - Risk: May be overfitting")

print(f"\n🎯 Clearest Assignments: {best_entropy['Model']} (entropy = {best_entropy['mean_entropy']:.4f})")
print(f"   - This model produces the most distinct archetypes")
print(f"   - {best_entropy['high_entropy_pct']:.1f}% high-entropy metacells")

print(f"\n📊 Most Balanced: {best_gini['Model']} (Gini = {best_gini['gini_coefficient']:.3f})")
print(f"   - This model distributes metacells most evenly across archetypes")
print(f"   - Avoids archetype collapse or duplication")

# Create visualization
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('VAE Archetype Model Comparison', fontsize=16, fontweight='bold')

# 1. Loss comparison
ax = axes[0, 0]
df_plot = df_comparison.sort_values('total_loss')
ax.barh(df_plot['Model'], df_plot['total_loss'], color='steelblue')
ax.set_xlabel('Total Loss (lower = better fit)')
ax.set_title('Training Loss')
ax.invert_yaxis()
for i, (model, loss) in enumerate(zip(df_plot['Model'], df_plot['total_loss'])):
    ax.text(loss + 0.005, i, f'{loss:.4f}', va='center')

# 2. Entropy comparison
ax = axes[0, 1]
df_plot = df_comparison.sort_values('mean_entropy')
ax.barh(df_plot['Model'], df_plot['mean_entropy'], color='coral')
ax.set_xlabel('Mean Entropy (lower = clearer assignments)')
ax.set_title('Assignment Clarity')
ax.invert_yaxis()
for i, (model, ent) in enumerate(zip(df_plot['Model'], df_plot['mean_entropy'])):
    ax.text(ent + 0.03, i, f'{ent:.3f}', va='center')

# 3. Gini coefficient
ax = axes[1, 0]
df_plot = df_comparison.sort_values('gini_coefficient')
ax.barh(df_plot['Model'], df_plot['gini_coefficient'], color='green')
ax.set_xlabel('Gini Coefficient (lower = more balanced)')
ax.set_title('Archetype Distribution Balance')
ax.invert_yaxis()
for i, (model, gini) in enumerate(zip(df_plot['Model'], df_plot['gini_coefficient'])):
    ax.text(gini + 0.005, i, f'{gini:.3f}', va='center')

# 4. High-entropy percentage
ax = axes[1, 1]
df_plot = df_comparison.sort_values('high_entropy_pct', ascending=False)
ax.barh(df_plot['Model'], df_plot['high_entropy_pct'], color='purple')
ax.set_xlabel('% High-Entropy Metacells (>1.5)')
ax.set_title('Mixed Phenotype Prevalence')
ax.invert_yaxis()
for i, (model, pct) in enumerate(zip(df_plot['Model'], df_plot['high_entropy_pct'])):
    ax.text(pct + 1, i, f'{pct:.1f}%', va='center')

plt.tight_layout()
plt.savefig(output_dir / 'model_comparison_charts.png', dpi=300, bbox_inches='tight')
print(f"\n✓ Saved visualization to {output_dir / 'model_comparison_charts.png'}")

# Recommendation
print("\n" + "="*80)
print("RECOMMENDATION")
print("="*80)

# Score each model (lower is better)
df_comparison['score'] = (
    df_comparison['total_loss'] / df_comparison['total_loss'].min() * 0.3 +  # 30% weight
    df_comparison['mean_entropy'] / df_comparison['mean_entropy'].min() * 0.4 +  # 40% weight
    df_comparison['gini_coefficient'] / df_comparison['gini_coefficient'].min() * 0.3  # 30% weight
)

best_overall = df_comparison.loc[df_comparison['score'].idxmin()]

print(f"\n🌟 RECOMMENDED MODEL: {best_overall['Model']}")
print(f"\n   Composite score: {best_overall['score']:.3f} (lower is better)")
print(f"\n   Key metrics:")
print(f"   - Total loss: {best_overall['total_loss']:.4f}")
print(f"   - Mean entropy: {best_overall['mean_entropy']:.4f}")
print(f"   - Gini coefficient: {best_overall['gini_coefficient']:.3f}")
print(f"   - High-entropy metacells: {best_overall['high_entropy_pct']:.1f}%")
print(f"   - Number of archetypes: {best_overall['n_archetypes']:.0f}")
print(f"\n   Configuration:")
print(f"   - Gamma: {best_overall.get('config_gamma', 'N/A')}")
print(f"   - Gamma cooldown: {best_overall.get('config_gamma_cooldown', 'N/A')} epochs")
print(f"   - Gamma min: {best_overall.get('config_gamma_min', 'N/A')}")
print(f"   - Beta: {best_overall.get('config_beta', 'N/A')}")
print(f"   - Latent dim: {best_overall.get('config_latent_dim', 'N/A')}")

print("\n" + "="*80)
print("NEXT STEPS")
print("="*80)
print(f"\n1. Use the recommended model for production analysis")
print(f"2. Apply to full dataset for sample composition analysis")
print(f"3. Proceed with external validation (see VALIDATION_ROADMAP.md)")

print("\n✓ Analysis complete!")
