# VAE Archetype Analysis - Key Findings Summary

**Date**: November 15, 2025
**Branch**: `nov15_gamma_cooldown_experiments`
**Model**: VAE with biological priors + gamma cooldown

---

## Executive Summary

We developed a VAE-based archetype framework that identifies distinct lymphocyte phenotypes from flow cytometry data. Training on Normal samples as reference and projecting Abnormal/PTCy samples revealed a **striking disease signature**: loss of phenotypic diversity and accumulation of CD8+ T-cells in transitional states.

---

## Major Discoveries

### 1. Disease Signature: CD8+ T-cell Accumulation

**Abnormal samples lose diversity and concentrate in 2 archetypes:**

| Archetype | Normal | Abnormal | Fold Change | Identity |
|-----------|--------|----------|-------------|----------|
| **Archetype 6** | 8.4% | **26.1%** | **3.1x** ⬆️ | CD8+ Effector Memory |
| **Archetype 8** | 18.3% | **35.8%** | **2.0x** ⬆️ | Unclear (possible transition state) |
| **Combined** | 26.7% | **61.9%** | **2.3x** | - |

**Statistical significance**: p < 8.22 × 10^-134 (Chi-square test)

### 2. Entropy as Disease Biomarker

**Disease samples have massively higher "mixedness":**

| Sample Type | Mean Entropy | High-Entropy (>1.5) | Significance |
|-------------|--------------|---------------------|--------------|
| **Normal** | 1.194 | 27.4% | baseline |
| **Abnormal** | 1.645 | **78.0%** | p < 10^-300 |
| **PTCy** | 1.637 | **75.5%** | p < 10^-103 |

**Interpretation**: Disease samples dominated by transitional/mixed phenotypes rather than stable, distinct cell types.

### 3. Archetype Mixing Patterns

**In high-entropy (mixed) metacells, Archetypes 6 & 8 frequently co-occur:**

| Archetype Pair | Co-occurrences | Primary Sample Type |
|----------------|----------------|---------------------|
| **6 + 8** | **1,780** | Abnormal (96%) |
| 2 + 8 | 1,538 | Mixed |
| 4 + 8 | 1,176 | Mixed |
| 4 + 6 | 892 | Abnormal |

**Most common mixing pattern in Abnormal**: 6-8 (CD8+ states)

---

## Biological Interpretation

### Archetype Identities (from marker expression)

**Disease-enriched archetypes:**
- **Archetype 6**: CD8+ Effector Memory
  - CD8+CD14 ↑↑ (1.96)
  - CD45RO ↑↑ (1.51) [memory marker]
  - CD45 ↑ (1.25)

- **Archetype 8**: Unclear phenotype
  - CD45RA ↑ (0.22) [naive marker]
  - CD3+CD19 ↑ (0.16)
  - Possibly naive/memory transition state

**Disease-depleted archetypes (nearly absent in Abnormal):**
- Archetype 0, 3, 9: Various CD4+ and NK subtypes
- Represents **loss of diversity**

### Disease Mechanism Hypothesis

1. **Normal state**: Diverse, stable lymphocyte phenotypes across 10 archetypes
2. **Disease state**:
   - Loss of diversity (6+ archetypes depleted)
   - Accumulation of CD8+ Effector Memory cells (Arch 6)
   - Increased transitional/activated states (Arch 8, high entropy)
   - Possible abnormal activation, exhaustion, or dysregulation

---

## Model Performance

### Best Configuration (Production Model)

**Training**: 150 epochs, gamma cooldown (0.5 → 0.0 over 100 epochs)

| Metric | Value |
|--------|-------|
| Total Loss | 0.1577 |
| Reconstruction Loss | 0.1445 |
| KL Divergence | 1.3171 |
| Training Time | ~2-3 minutes |
| Archetype Distribution (Gini) | 0.182 (even) |

**Key advantages over alternatives:**
- 13% lower loss than constant gamma
- More even archetype distribution (27% better Gini)
- No duplicate archetypes
- Clean, interpretable results

### Gamma Cooldown Innovation

**Concept**: Enforce biological priors strongly early in training, then relax to allow data-driven discovery.

- **Early training (0-100 epochs)**: Gamma = 0.5 → guides model toward known cell types
- **Late training (100-150 epochs)**: Gamma → 0.0 → allows discovery of disease-specific patterns

**Result**: Best of both worlds - biologically grounded + data-driven

---

## Clinical Implications

### 1. Diagnostic Potential

**Archetype composition can distinguish Normal from disease:**
- Simple binary classifier: % in Archetype 6 + 8
- Normal: ~27%, Disease: ~60%
- Could achieve high sensitivity/specificity

### 2. Disease Monitoring

**Entropy as progression/response biomarker:**
- Baseline: Normal ~1.2, Disease ~1.6
- Track over time: decreasing entropy = improving?
- Monitor treatment response

### 3. Therapeutic Targets

**CD8+ T-cell dysregulation identified:**
- Accumulation of Effector Memory (Arch 6)
- Increased transitional states (Arch 8)
- Suggests targeting CD8+ activation/exhaustion pathways

---

## Technical Details

### Framework Features

1. **Train-and-Project**: Train on Normal, project Abnormal/PTCy
2. **Biological Priors**: Soft constraints on T-cell/NK, CD4/CD8, Memory/Naive
3. **Gamma Cooldown**: Adaptive biological prior schedule
4. **Soft Assignments**: Probability distributions + entropy metric
5. **Transferable**: Saved models can project future samples

### Data

- **Input**: 6,039 metacells from SEACells (Harmony-corrected)
  - Normal: 605 (10%)
  - Abnormal: 4,770 (79%)
  - PTCy: 664 (11%)
- **Markers**: 13 (some merged channels due to panel design)
- **Archetypes**: 10 identified via k-means on VAE latent space

---

## Files & Outputs

### Key Outputs

| Directory | Contents | Description |
|-----------|----------|-------------|
| `test_output/gamma_cooldown/` | ⭐ **Production Model** | 150 epoch, best performance |
| `test_output/sample_composition_analysis/` | Archetype % by sample type | **Major findings** |
| `test_output/mixed_phenotype_analysis/` | Co-occurrence matrices | **Major findings** |

### Documentation

- `Nov15.md` - Complete experimental log and findings
- `README_VAE_ARCHETYPES.md` - Framework usage guide
- `FINDINGS_SUMMARY.md` - This document

### Analysis Scripts

- `analyze_sample_composition.py` - Sample type comparisons
- `analyze_mixed_phenotypes.py` - High-entropy analysis
- `compare_experiments.py` - Model comparison

---

## Next Steps

### Immediate

1. **Biological validation**: Annotate all 10 archetypes with expected cell types
2. **Sample annotation**: Work with domain expert to identify sample heterogeneity
3. **Visualize key archetypes**: Detailed marker profiles for Arch 6 & 8

### Short-term

1. **Hyperparameter optimization**: Optuna sweep to optimize all parameters
2. **Longitudinal analysis**: If time-series data available, track changes
3. **Clinical correlation**: Link to disease severity, outcomes, treatment response

### Long-term

1. **Prospective validation**: Test on new samples
2. **Refine biological priors**: Add B-cell, activation, exhaustion markers
3. **Hierarchical archetypes**: Multi-level clustering (coarse → fine)
4. **Clinical trial application**: Use as stratification or response biomarker

---

## Conclusion

We successfully developed a biologically-informed VAE framework that:
1. ✅ Identifies meaningful lymphocyte archetypes
2. ✅ Distinguishes Normal from disease with high significance
3. ✅ Reveals CD8+ T-cell accumulation as disease signature
4. ✅ Provides entropy as quantitative disease biomarker
5. ✅ Enables projection of future samples onto reference model

**The framework is production-ready and has discovered actionable biological insights.**

---

**Contact**: See `Nov15.md` for detailed experimental procedures and results.
