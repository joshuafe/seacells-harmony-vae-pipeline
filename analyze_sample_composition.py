"""
Analyze archetype composition differences between sample types

Compares Normal (reference) vs Abnormal vs PTCy archetype distributions
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path

def analyze_sample_composition(adata_path, output_dir):
    """Analyze and visualize archetype composition by sample type"""

    print("="*80)
    print("ARCHETYPE COMPOSITION BY SAMPLE TYPE")
    print("="*80)

    # Load data
    print(f"\nLoading: {adata_path}")
    adata = sc.read_h5ad(adata_path)

    # Get sample types
    sample_types = adata.obs['sample_type'].unique()
    print(f"Sample types: {list(sample_types)}")
    print(f"  Total metacells: {adata.n_obs}")
    for st in sample_types:
        n = (adata.obs['sample_type'] == st).sum()
        print(f"    {st}: {n} ({n/adata.n_obs*100:.1f}%)")

    # 1. Archetype composition by sample type
    print("\n" + "="*80)
    print("1. ARCHETYPE DISTRIBUTION BY SAMPLE TYPE")
    print("="*80)

    composition = []
    n_archetypes = adata.obs['archetype_id_hard'].nunique()

    for sample_type in sample_types:
        mask = adata.obs['sample_type'] == sample_type
        row = {'Sample_Type': sample_type, 'N_Metacells': mask.sum()}

        for arch_id in range(n_archetypes):
            count = (adata.obs.loc[mask, 'archetype_id_hard'] == arch_id).sum()
            pct = count / mask.sum() * 100
            row[f'Arch_{arch_id}_count'] = count
            row[f'Arch_{arch_id}_pct'] = pct

        composition.append(row)

    comp_df = pd.DataFrame(composition)

    # Print percentage composition
    print("\nPercentage of each archetype by sample type:")
    pct_cols = [c for c in comp_df.columns if c.endswith('_pct')]
    pct_df = comp_df[['Sample_Type'] + pct_cols].set_index('Sample_Type')
    pct_df.columns = [c.replace('_pct', '') for c in pct_df.columns]
    print(pct_df.to_string())

    # 2. Statistical comparison (chi-square test)
    print("\n" + "="*80)
    print("2. STATISTICAL COMPARISON (Chi-square test)")
    print("="*80)

    # Create contingency table
    contingency = []
    for sample_type in sample_types:
        mask = adata.obs['sample_type'] == sample_type
        counts = [
            (adata.obs.loc[mask, 'archetype_id_hard'] == i).sum()
            for i in range(n_archetypes)
        ]
        contingency.append(counts)

    chi2, p_value, dof, expected = stats.chi2_contingency(contingency)

    print(f"Chi-square statistic: {chi2:.2f}")
    print(f"P-value: {p_value:.2e}")
    print(f"Degrees of freedom: {dof}")

    if p_value < 0.001:
        print("\n*** HIGHLY SIGNIFICANT DIFFERENCE between sample types (p < 0.001) ***")
    elif p_value < 0.05:
        print("\n** Significant difference between sample types (p < 0.05) **")
    else:
        print("\nNo significant difference between sample types")

    # 3. Identify enriched/depleted archetypes
    print("\n" + "="*80)
    print("3. ENRICHED/DEPLETED ARCHETYPES (vs Normal)")
    print("="*80)

    normal_pcts = pct_df.loc['Normal'].values

    for sample_type in ['Abnormal', 'PTCy']:
        if sample_type in pct_df.index:
            print(f"\n{sample_type} vs Normal:")
            test_pcts = pct_df.loc[sample_type].values

            for i, (norm_pct, test_pct) in enumerate(zip(normal_pcts, test_pcts)):
                fold_change = test_pct / (norm_pct + 1e-10)  # avoid div by zero
                diff = test_pct - norm_pct

                if abs(diff) > 2.0:  # >2% difference
                    direction = "↑ ENRICHED" if diff > 0 else "↓ DEPLETED"
                    print(f"  Archetype {i}: {direction} ({test_pct:.1f}% vs {norm_pct:.1f}%, {diff:+.1f}%, {fold_change:.2f}x)")

    # 4. Entropy comparison
    print("\n" + "="*80)
    print("4. ENTROPY COMPARISON (mixedness)")
    print("="*80)

    for sample_type in sample_types:
        mask = adata.obs['sample_type'] == sample_type
        entropy = adata.obs.loc[mask, 'archetype_entropy']

        print(f"\n{sample_type}:")
        print(f"  Mean entropy: {entropy.mean():.3f}")
        print(f"  Median entropy: {entropy.median():.3f}")
        print(f"  Std: {entropy.std():.3f}")
        print(f"  High-entropy (>1.5): {(entropy > 1.5).sum()} ({(entropy > 1.5).sum()/len(entropy)*100:.1f}%)")

    # Statistical test for entropy differences
    normal_entropy = adata.obs.loc[adata.obs['sample_type'] == 'Normal', 'archetype_entropy']

    for sample_type in ['Abnormal', 'PTCy']:
        if sample_type in sample_types:
            test_entropy = adata.obs.loc[adata.obs['sample_type'] == sample_type, 'archetype_entropy']
            t_stat, p_val = stats.ttest_ind(normal_entropy, test_entropy)
            print(f"\n{sample_type} vs Normal t-test: t={t_stat:.2f}, p={p_val:.2e}")
            if p_val < 0.001:
                direction = "higher" if t_stat < 0 else "lower"
                print(f"  *** {sample_type} has SIGNIFICANTLY {direction} entropy than Normal ***")

    # 5. Visualizations
    print("\n" + "="*80)
    print("5. CREATING VISUALIZATIONS")
    print("="*80)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(20, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

    # Heatmap of archetype percentages
    ax1 = fig.add_subplot(gs[0, :])
    sns.heatmap(pct_df.T, annot=True, fmt='.1f', cmap='YlOrRd',
                cbar_kws={'label': 'Percentage'}, ax=ax1)
    ax1.set_title('Archetype Composition by Sample Type (%)', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Sample Type')
    ax1.set_ylabel('Archetype')

    # Stacked bar chart
    ax2 = fig.add_subplot(gs[1, 0])
    arch_cols = [f'Arch_{i}' for i in range(n_archetypes)]
    pct_df[arch_cols].T.plot(kind='bar', stacked=True, ax=ax2, legend=False)
    ax2.set_title('Stacked Archetype Distribution')
    ax2.set_xlabel('Archetype')
    ax2.set_ylabel('Percentage')
    ax2.set_xticklabels(range(n_archetypes), rotation=0)
    ax2.legend(sample_types, loc='upper right', fontsize=8)

    # Grouped bar chart
    ax3 = fig.add_subplot(gs[1, 1])
    pct_df[arch_cols].plot(kind='bar', ax=ax3, width=0.8, alpha=0.8)
    ax3.set_title('Grouped Archetype Distribution')
    ax3.set_xlabel('Sample Type')
    ax3.set_ylabel('Percentage')
    ax3.set_xticklabels(ax3.get_xticklabels(), rotation=0)
    ax3.legend(arch_cols, ncol=2, fontsize=7, title='Archetype')
    ax3.grid(True, alpha=0.3, axis='y')

    # Entropy distributions
    ax4 = fig.add_subplot(gs[1, 2])
    for sample_type in sample_types:
        mask = adata.obs['sample_type'] == sample_type
        entropy = adata.obs.loc[mask, 'archetype_entropy']
        ax4.hist(entropy, bins=30, alpha=0.5, label=sample_type, density=True)
    ax4.set_xlabel('Archetype Entropy')
    ax4.set_ylabel('Density')
    ax4.set_title('Entropy Distribution by Sample Type')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    # Fold change heatmap (vs Normal)
    ax5 = fig.add_subplot(gs[2, :2])
    fold_changes = []
    for sample_type in ['Abnormal', 'PTCy']:
        if sample_type in pct_df.index:
            fc = pct_df.loc[sample_type, arch_cols].values / (pct_df.loc['Normal', arch_cols].values + 1e-10)
            fold_changes.append(fc)

    if fold_changes:
        fc_df = pd.DataFrame(
            fold_changes,
            index=[st for st in ['Abnormal', 'PTCy'] if st in pct_df.index],
            columns=arch_cols
        )

        sns.heatmap(fc_df, annot=True, fmt='.2f', cmap='RdBu_r', center=1.0,
                   vmin=0, vmax=3, cbar_kws={'label': 'Fold Change vs Normal'}, ax=ax5)
        ax5.set_title('Archetype Enrichment (Fold Change vs Normal)', fontsize=12, fontweight='bold')
        ax5.set_xlabel('Archetype')
        ax5.set_ylabel('Sample Type')

    # Scatter: entropy vs archetype diversity
    ax6 = fig.add_subplot(gs[2, 2])
    for sample_type in sample_types:
        mask = adata.obs['sample_type'] == sample_type
        entropy = adata.obs.loc[mask, 'archetype_entropy']

        # Compute archetype diversity (number of archetypes with >5% probability)
        diversity = []
        for idx in adata.obs.index[mask]:
            probs = [adata.obs.loc[idx, f'archetype_{i}_prob'] for i in range(n_archetypes)]
            n_diverse = sum([p > 0.05 for p in probs])
            diversity.append(n_diverse)

        ax6.scatter(diversity, entropy, alpha=0.3, s=10, label=sample_type)

    ax6.set_xlabel('Archetype Diversity (# with >5% prob)')
    ax6.set_ylabel('Entropy')
    ax6.set_title('Entropy vs Archetype Diversity')
    ax6.legend()
    ax6.grid(True, alpha=0.3)

    plt.savefig(output_dir / 'sample_type_composition.png', dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_dir / 'sample_type_composition.png'}")

    # Save composition table
    comp_df.to_csv(output_dir / 'archetype_composition_by_sample.csv', index=False)
    print(f"  ✓ Saved: {output_dir / 'archetype_composition_by_sample.csv'}")

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Chi-square p-value: {p_value:.2e} - Sample types have {'SIGNIFICANTLY' if p_value < 0.001 else 'somewhat'} different archetype distributions")

    print("\nKey findings will be in visualizations and CSV output.")
    print("="*80)

if __name__ == "__main__":
    # Use the 150 epoch model (better interpretability)
    analyze_sample_composition(
        'test_output/gamma_cooldown/metacells_with_archetypes.h5ad',
        'test_output/sample_composition_analysis'
    )
