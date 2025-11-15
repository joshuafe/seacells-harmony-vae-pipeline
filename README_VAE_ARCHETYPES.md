# VAE Archetypes Framework

## Overview

This framework provides a comprehensive approach to identifying and analyzing lymphocyte archetypes using Variational Autoencoders (VAE) with biological priors.

### Key Features

1. **Train-and-Project Framework**: Train on reference samples (e.g., Normal) and project others (Abnormal, PTCy)
2. **Biological Prior Loss Functions**: Enforce known lineage relationships despite merged flow cytometry channels
3. **Soft Archetype Assignments**: Quantify "mixed" phenotypes using probability distributions and entropy
4. **Hyperparameter Tuning**: Optuna-based optimization to find optimal model configuration
5. **Future Sample Projection**: Project new samples onto established archetype model

---

## Scripts

### 1. `vae_archetypes_with_priors.py`

Main training script with biological priors.

**Usage:**

```bash
# Train on all samples
python scripts/vae_archetypes_with_priors.py \
  --input seacells_output/three_groups_phenograph/integrated_metacells_harmony_phenograph.h5ad \
  --latent-dim 5 \
  --n-archetypes 10 \
  --beta 0.01 \
  --gamma 0.1 \
  --epochs 300 \
  --output-dir vae_output_test

# Train on Normal samples only, project others
python scripts/vae_archetypes_with_priors.py \
  --input seacells_output/three_groups_phenograph/integrated_metacells_harmony_phenograph.h5ad \
  --train-on-reference \
  --reference-type Normal \
  --latent-dim 5 \
  --n-archetypes 10 \
  --gamma 0.5 \
  --output-dir vae_reference_model
```

**Key Parameters:**

- `--train-on-reference`: Train only on reference samples
- `--reference-type`: Which sample type to use as reference (default: Normal)
- `--latent-dim`: Size of latent space (4-8 recommended)
- `--n-archetypes`: Number of archetypes to identify (8-15 recommended)
- `--beta`: KL divergence weight (0.001-0.1, lower = less regularization)
- `--gamma`: Biological prior weight - INITIAL value (0-1, higher = stronger constraints)
- `--gamma-cooldown`: Epochs to cool down gamma (default: 150, 0 = no cooldown)
- `--gamma-min`: Minimum gamma after cooldown (default: 0.0)
- `--soft-assignment-temp`: Temperature for soft assignments (0.5-2.0, lower = sharper)

**Outputs:**

- `metacells_with_archetypes.h5ad`: Annotated AnnData object
- `vae_model.pt`: Trained model checkpoint
- `archetype_profiles.csv`: Marker expression profiles of archetypes
- `training_curves.png`: Training loss curves
- `vae_archetypes_overview.png`: UMAP visualizations

---

### 2. `tune_vae_archetypes.py`

Hyperparameter optimization using Optuna.

**Usage:**

```bash
python scripts/tune_vae_archetypes.py \
  --input seacells_output/three_groups_phenograph/integrated_metacells_harmony_phenograph.h5ad \
  --train-on-reference \
  --n-trials 50 \
  --epochs 200 \
  --output-dir optuna_tuning
```

**Optimized Metrics:**

1. **Assignment clarity**: Average entropy (lower = more distinct assignments)
2. **Latent space quality**: Silhouette score (higher = better separation)

**Outputs:**

- `tuning_results.csv`: All trial results
- `best_config.json`: Best hyperparameters
- `trial_*/`: Individual trial outputs

---

### 3. `project_samples.py`

Project new samples onto a trained model.

**Usage:**

```bash
# After training a reference model
python scripts/project_samples.py \
  --model vae_reference_model/vae_model.pt \
  --input new_samples.h5ad \
  --output-dir projection_results
```

**Use Cases:**

- Evaluate new patient samples against reference
- Compare samples without retraining
- Maintain consistent archetype definitions over time

**Outputs:**

- `projected_samples.h5ad`: Annotated samples
- `projection_overview.png`: Visualizations
- `archetype_composition.csv`: Sample composition summary

---

## Biological Priors

### Gamma Cooldown Schedule

**NEW**: The framework now supports a gamma cooldown schedule that enforces biological priors strongly at the beginning of training and gradually relaxes them:

- **Rationale**: Guide the model toward biologically plausible regions early (when parameters are random), then allow it to find nuanced patterns later
- **Default behavior**: Linear decay from `--gamma` to `--gamma-min` (default 0.0) over 150 epochs
- **Control**:
  - `--gamma 0.5`: Initial strong biological constraints
  - `--gamma-min 0.0`: No constraints at end of training
  - `--gamma-cooldown 150`: Decay over 150 epochs
  - `--gamma-cooldown 0`: Constant gamma (no cooldown)

**Example schedules:**

```bash
# Strong initial guidance, then full freedom
--gamma 0.5 --gamma-cooldown 150 --gamma-min 0.0

# Moderate initial guidance, maintain weak constraints
--gamma 0.3 --gamma-cooldown 100 --gamma-min 0.05

# Constant moderate constraints
--gamma 0.2 --gamma-cooldown 0
```

### Biological Constraints

The framework implements soft constraints for known biological relationships:

### 1. T-cell vs NK Mutual Exclusion

- **Constraint**: CD3+ (T-cells) and CD16+/CD56+ (NK cells) should be mutually exclusive
- **Implementation**: Penalize high co-expression of CD3+CD19 AND (CD16+TIGIT OR CD56)

### 2. CD4 vs CD8 Preference

- **Constraint**: T-cells should be predominantly CD4+ OR CD8+, not both or neither
- **Implementation**: For high CD3 cells, penalize intermediate CD4/CD8 expression

### 3. Memory vs Naive Exclusion

- **Constraint**: CD45RA (naive) and CD45RO (memory) should be anti-correlated
- **Implementation**: Penalize both being high simultaneously

**Note**: These are soft constraints that can be violated if there's strong signal, controlled by the `--gamma` parameter.

---

## Soft Archetype Assignments

Unlike hard assignments (each cell → one archetype), soft assignments provide:

1. **Probability distribution**: Each metacell has probabilities for all archetypes
2. **Entropy**: Measures "mixedness" (high entropy = mixed phenotype)
3. **Biological interpretation**: Metacells transitioning between states or with dual markers

**Accessing soft assignments:**

```python
import scanpy as sc

adata = sc.read_h5ad("vae_output/metacells_with_archetypes.h5ad")

# Hard assignment
archetype_id = adata.obs['archetype_id_hard']

# Soft assignments (probabilities)
prob_archetype_0 = adata.obs['archetype_0_prob']
prob_archetype_1 = adata.obs['archetype_1_prob']
# ...

# Entropy (mixedness)
entropy = adata.obs['archetype_entropy']

# Find mixed phenotype metacells (high entropy)
mixed = adata[adata.obs['archetype_entropy'] > 1.5]
```

---

## Recommended Workflow

### Stage 1: Initial Exploration

```bash
# Quick run with default parameters
python scripts/vae_archetypes_with_priors.py \
  --input <your_data.h5ad> \
  --train-on-reference \
  --epochs 200 \
  --output-dir initial_test
```

Review outputs, check if archetypes are biologically meaningful.

### Stage 2: Hyperparameter Tuning

```bash
# Optimize hyperparameters
python scripts/tune_vae_archetypes.py \
  --input <your_data.h5ad> \
  --train-on-reference \
  --n-trials 50 \
  --output-dir tuning_results
```

Check `best_config.json` for optimal parameters.

### Stage 3: Final Training

```bash
# Train with optimized parameters
python scripts/vae_archetypes_with_priors.py \
  --input <your_data.h5ad> \
  --train-on-reference \
  --latent-dim <from_best_config> \
  --n-archetypes <from_best_config> \
  --beta <from_best_config> \
  --gamma <from_best_config> \
  --epochs 500 \
  --output-dir final_model
```

### Stage 4: Analysis

Analyze archetypes, mixed phenotypes, sample comparisons.

### Stage 5: Future Projections

```bash
# Project new samples
python scripts/project_samples.py \
  --model final_model/vae_model.pt \
  --input new_samples.h5ad \
  --output-dir new_projection
```

---

## Analysis Examples

### 1. Identify Mixed Phenotype Metacells

```python
import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("vae_output/metacells_with_archetypes.h5ad")

# High entropy = mixed phenotypes
threshold = adata.obs['archetype_entropy'].quantile(0.9)
mixed = adata[adata.obs['archetype_entropy'] > threshold]

print(f"Mixed phenotype metacells: {mixed.n_obs}")

# Look at their top 2 archetype associations
for i in mixed.obs.index[:10]:
    probs = [adata.obs.loc[i, f'archetype_{j}_prob'] for j in range(10)]
    top2 = sorted(enumerate(probs), key=lambda x: x[1], reverse=True)[:2]
    print(f"Metacell {i}: Archetype {top2[0][0]} ({top2[0][1]:.2f}), Archetype {top2[1][0]} ({top2[1][1]:.2f})")
```

### 2. Compare Sample Types

```python
import matplotlib.pyplot as plt
import seaborn as sns

# Archetype composition by sample type
composition = []
for sample_type in ['Normal', 'Abnormal', 'PTCy']:
    mask = adata.obs['sample_type'] == sample_type
    row = {'Sample_Type': sample_type}
    for i in range(10):
        pct = (adata.obs.loc[mask, 'archetype_id_hard'] == i).sum() / mask.sum() * 100
        row[f'Archetype_{i}'] = pct
    composition.append(row)

comp_df = pd.DataFrame(composition).set_index('Sample_Type')

plt.figure(figsize=(12, 6))
sns.heatmap(comp_df, cmap='YlOrRd', annot=True, fmt='.1f')
plt.title('Archetype Composition by Sample Type')
plt.savefig('composition_heatmap.png', dpi=300, bbox_inches='tight')
```

### 3. Measure Deviation from Normal

```python
# For each Abnormal/PTCy metacell, compute distance to nearest Normal metacell in latent space
from scipy.spatial.distance import cdist

normal_mask = adata.obs['sample_type'] == 'Normal'
abnormal_mask = adata.obs['sample_type'] == 'Abnormal'

normal_latent = adata[normal_mask].obsm['X_vae']
abnormal_latent = adata[abnormal_mask].obsm['X_vae']

# Compute distances
distances = cdist(abnormal_latent, normal_latent, metric='euclidean')
min_distances = distances.min(axis=1)

adata.obs.loc[abnormal_mask, 'distance_to_normal'] = min_distances

# High distance = divergent from normal
divergent = adata[abnormal_mask][min_distances > min_distances.quantile(0.9)]
print(f"Divergent abnormal metacells: {divergent.n_obs}")
```

---

## Questions Addressed

### Q1: Should I train on all samples or a subset?

**Answer**: Train on Normal samples (reference), project others.

- **Pros**: Normals define the "baseline", deviations are measured relative to reference
- **Cons**: Might miss archetypes unique to Abnormal/PTCy
- **Solution**: Use `--train-on-reference` flag

### Q2: How do I quantify "mixed" phenotypes?

**Answer**: Use soft archetype assignments and entropy.

- High entropy = metacell has significant probability across multiple archetypes
- Access via `adata.obs['archetype_entropy']` and `adata.obs['archetype_*_prob']`

### Q3: How do I enforce biological priors?

**Answer**: Use the `--gamma` parameter.

- `--gamma 0`: No biological constraints (pure VAE)
- `--gamma 0.1-0.5`: Soft guidance toward known biology
- `--gamma 1.0+`: Strong enforcement (may be too restrictive)

### Q4: Do I still need a two-stage model (lymphocyte gating → archetypes)?

**Answer**: Depends on your data.

- If using SEACells metacells already enriched for lymphocytes: **Single stage is fine**
- If working with raw data including monocytes, granulocytes, etc.: **Two-stage recommended**

The current framework assumes lymphocyte-enriched data (CD45+ cells).

### Q5: How do I project future samples?

**Answer**: Use `project_samples.py` with the trained model.

This maintains consistent archetype definitions and enables longitudinal comparisons.

---

## Next Steps

1. **Run initial test** with `vae_archetypes_with_priors.py`
2. **Hyperparameter tuning** with `tune_vae_archetypes.py`
3. **Train final model** with optimized parameters
4. **Analyze archetypes**: Are they biologically meaningful? Do they align with expected markers?
5. **Compare samples**: Normal vs Abnormal vs PTCy
6. **Quantify mixed phenotypes**: Which metacells have high entropy? What do they represent?
7. **Iterate on biological priors**: Adjust `--gamma` based on domain knowledge

---

## Citation

If you use this framework, please cite:
- VAE architecture inspired by Kingma & Welling (2013)
- SEACells metacell generation (Persad et al., 2022)
- Harmony batch correction (Korsunsky et al., 2019)
