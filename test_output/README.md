# Experiment Outputs

This directory contains results from gamma cooldown experiments conducted on November 15, 2025.

## Experiments

### Experiment 1: Gamma Cooldown
**Directory**: `gamma_cooldown/`

**Configuration**:
- Train on reference: Normal samples only (n=605)
- Latent dimensions: 5
- Number of archetypes: 10
- Beta (KL weight): 0.01
- **Gamma schedule**: 0.5 → 0.0 over 100 epochs (cooldown)
- Total epochs: 150

**Results**:
- Final loss: 0.1577
- Reconstruction loss: 0.1445
- Mean entropy (reference): 1.194
- Gini coefficient: 0.182 (more even distribution)

**Key Finding**: Better optimization, cleaner archetypes, fewer duplicates

---

### Experiment 2: Constant Gamma
**Directory**: `constant_gamma/`

**Configuration**:
- Train on reference: Normal samples only (n=605)
- Latent dimensions: 5
- Number of archetypes: 10
- Beta (KL weight): 0.01
- **Gamma**: 0.1 (constant, no cooldown)
- Total epochs: 150

**Results**:
- Final loss: 0.1756
- Reconstruction loss: 0.1485
- Mean entropy (reference): 1.129
- Gini coefficient: 0.232 (less even distribution)

**Key Finding**: Creates duplicate memory archetypes (Arch 2 & 4 nearly identical)

---

## Comparison

**File**: `experiment_comparison.png`

Side-by-side comparison showing:
1. Entropy distributions
2. Archetype distribution evenness
3. Archetype profile correlation matrix
4. Mean entropy by sample type

---

## Verdict

✅ **Gamma cooldown is superior**

**Advantages**:
- 13% lower total loss
- Better reconstruction
- 27% more even archetype distribution
- Fewer duplicate archetypes
- Biologically guided initialization + data-driven refinement

**Recommendation**: Use gamma cooldown (0.5→0.0) as default for production models

---

## File Structure

Each experiment directory contains:
- `metacells_with_archetypes.h5ad` - Annotated AnnData object with latent embeddings and archetype assignments
- `archetype_profiles.csv` - Marker expression profiles for each archetype
- `vae_model.pt` - Saved PyTorch model (can be loaded for future projections)
- `training_curves.png` - Training loss curves including gamma schedule
- `vae_archetypes_overview.png` - UMAP visualizations and archetype heatmap

---

## Reproducibility

To reproduce these experiments:

```bash
# Experiment 1 (Gamma Cooldown)
python scripts/vae_archetypes_with_priors.py \
  --input seacells_output/three_groups_phenograph/integrated_metacells_harmony_phenograph.h5ad \
  --train-on-reference \
  --latent-dim 5 \
  --n-archetypes 10 \
  --beta 0.01 \
  --gamma 0.5 \
  --gamma-cooldown 100 \
  --gamma-min 0.0 \
  --epochs 150 \
  --output-dir test_output/gamma_cooldown

# Experiment 2 (Constant Gamma)
python scripts/vae_archetypes_with_priors.py \
  --input seacells_output/three_groups_phenograph/integrated_metacells_harmony_phenograph.h5ad \
  --train-on-reference \
  --latent-dim 5 \
  --n-archetypes 10 \
  --beta 0.01 \
  --gamma 0.1 \
  --gamma-cooldown 0 \
  --epochs 150 \
  --output-dir test_output/constant_gamma

# Comparison analysis
python compare_experiments.py
```

Or run both with:
```bash
./test_gamma_cooldown.sh
```

---

**Date**: November 15, 2025
**Git Branch**: `nov15_gamma_cooldown_experiments`
**Documentation**: See `Nov15.md` for detailed analysis
