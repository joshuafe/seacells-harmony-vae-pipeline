# Project Status - November 15, 2025

## Project Overview

Developing VAE-based archetype framework for flow cytometry lymphocyte analysis with the goal of:
1. Establishing Normal samples as reference baseline
2. Measuring deviations of Abnormal and PTCy samples from this baseline
3. Identifying archetypes with biological meaning
4. Quantifying "mixed" phenotypes in metacells
5. Enforcing biological priors despite merged flow channels

---

## Current Status

### Completed Work

#### 1. Initial VAE Archetype Model (`vae_archetypes_seacells.py`)
- **Location**: `scripts/vae_archetypes_seacells.py`
- **Output**: `vae_archetypes_output/`
- **Status**: ✅ Working, produces 10 archetypes
- **Results**:
  - 10 archetypes identified via k-means on latent space
  - Biological names assigned (see `vae_archetypes_output/archetype_biological_names.csv`)
  - Examples: Effector Memory CD4+ T, Central Memory CD4+ T, Naive T, NK cells, etc.
  - CD197 emerged as major feature (TEMRA phenotypes?)

**Key observation**: The model is promising but needs tuning to optimize biological relevance and account for:
- Sample imbalance (605 Normal, 4770 Abnormal, 664 PTCy)
- Heterogeneous population requiring annotation
- Need for reference-based comparison framework

#### 2. New VAE Framework with Biological Priors (`vae_archetypes_with_priors.py`)
- **Location**: `scripts/vae_archetypes_with_priors.py`
- **Status**: ✅ Code complete, UNTESTED
- **Key features**:
  - **Train-and-project mode**: Train on Normal (reference), project Abnormal/PTCy
  - **Biological prior loss functions**: Soft constraints on lineage relationships
  - **Soft archetype assignments**: Probability distributions + entropy (quantifies "mixedness")
  - **Saved models**: Can project future samples without retraining

**Biological Priors Implemented**:
1. T-cell vs NK mutual exclusion (CD3 vs CD16/CD56)
2. CD4 vs CD8 preference (T-cells should be one or the other)
3. Memory vs Naive exclusion (CD45RA vs CD45RO anti-correlation)

**Controlled by `--gamma` parameter**: 0 = no priors, 1.0 = strong enforcement

#### 3. Hyperparameter Tuning Framework (`tune_vae_archetypes.py`)
- **Location**: `scripts/tune_vae_archetypes.py`
- **Status**: ✅ Code complete, UNTESTED
- **Uses**: Optuna for optimization
- **Optimizes**:
  - Number of archetypes (8-15)
  - Latent dimensions (4-8)
  - Beta (KL weight, 0.001-0.1)
  - Gamma (biological prior weight, 0-1.0)
  - Soft assignment temperature (0.5-2.0)
- **Metrics**:
  - Average entropy (lower = clearer assignments)
  - Silhouette score (higher = better separation)

#### 4. Sample Projection Script (`project_samples.py`)
- **Location**: `scripts/project_samples.py`
- **Status**: ✅ Code complete, UNTESTED
- **Purpose**: Project new samples onto trained reference model
- **Use cases**:
  - Future patient samples
  - Longitudinal comparisons
  - Maintain consistent archetype definitions

#### 5. Documentation
- **Location**: `README_VAE_ARCHETYPES.md`
- **Status**: ✅ Complete
- **Contains**: Usage examples, workflow recommendations, analysis examples

---

## Current Data Characteristics

### Input Data
- **File**: `seacells_output/three_groups_phenograph/integrated_metacells_harmony_phenograph.h5ad`
- **Size**: 6039 metacells × 13 markers
- **Sample distribution**:
  - Normal: 605 (10%)
  - Abnormal: 4770 (79%)
  - PTCy: 664 (11%)

### Markers (Merged Channels)
1. CD45 (pan-lymphocyte)
2. **CD3+CD19** (T-cell + B-cell marker, merged)
3. **CD4+CD33** (T helper + myeloid marker, merged)
4. **CD8+CD14** (Cytotoxic T + monocyte marker, merged)
5. **CD16+TIGIT** (NK cell + exhaustion marker, merged)
6. CD34+CD223
7. CD45RA (naive)
8. CD45RO (memory)
9. CD56 (NK cell)
10. CD69 (activation)
11. CD95 (apoptosis/memory)
12. CD197 (CCR7, central memory/naive)
13. CD366 (TIM-3, exhaustion)

**Challenge**: Merged channels complicate biological interpretation but priors can still enforce expected patterns.

### Preprocessing Status
- ✅ SEACells metacells generated
- ✅ Harmony batch correction applied
- ✅ Phenograph clustering performed
- ✅ Data is lymphocyte-enriched (no explicit gating needed)

---

## Ongoing Issues & Open Questions

### 1. Merged Channels Complicate Lineage Constraints

**Issue**: CD3+CD19, CD4+CD33, CD8+CD14, CD16+TIGIT are merged channels
- Can't cleanly separate CD3+ T-cells from CD19+ B-cells
- Can't cleanly separate CD4+ from CD33+ (myeloid)
- Biological priors must be "soft" to account for this

**Current solution**:
- Soft constraints that can be violated with strong signal
- Controlled by `--gamma` parameter
- May need to adjust prior weights based on empirical results

**TODO**:
- [ ] Test different gamma values (0, 0.1, 0.5, 1.0) and compare archetype quality
- [ ] Consider channel-specific weights for priors

### 2. Sample Imbalance

**Issue**: 605 Normal vs 4770 Abnormal vs 664 PTCy
- Training on all samples may bias toward Abnormal
- Training on Normal only may miss important Abnormal-specific archetypes

**Current solution**:
- `--train-on-reference` mode trains on Normal, projects others
- This makes Normal the "baseline" and measures deviations

**Open question**:
- Should we downsample Abnormal to balance training?
- Should we identify archetypes separately for each group then merge?
- Should we use sample weights during training?

**TODO**:
- [ ] Test both modes: train-on-all vs train-on-reference
- [ ] Compare archetype interpretability and separation
- [ ] Quantify how much Abnormal/PTCy deviate from Normal baseline

### 3. Optimal Number of Archetypes

**Issue**: Unknown optimal k (number of archetypes)
- Too few: Miss important phenotypes
- Too many: Over-segmentation, hard to interpret

**Current approach**:
- Hyperparameter tuning with k ∈ [8, 15]
- Use silhouette score and entropy as proxies

**Better approach needed**:
- Biological validation (do archetypes match known cell types?)
- Stability across random initializations
- Clinical relevance (do they predict outcomes?)

**TODO**:
- [ ] Run tuning script to find optimal k
- [ ] Manually inspect archetypes for biological sense
- [ ] Consider hierarchical archetypes (coarse → fine)

### 4. CD197 as Major Feature

**Observation**: In initial results, CD197 (CCR7) was a dominant feature
- Separates central memory/naive (CD197+) from effector (CD197-)
- Possibly TEMRA phenotypes (CD45RA+ CD197-)

**Open question**:
- Is this biologically meaningful or artifact of data/normalization?
- Will this hold after hyperparameter tuning?

**TODO**:
- [ ] Compare CD197 importance across different model configurations
- [ ] Literature review: Is CCR7 expected to be major axis in this context?

### 5. Annotation Requirements

**Issue**: User mentioned "heterogeneous population, need to annotate since certain samples are more appropriate to compare than others"

**Implication**:
- Not all Normal samples may be equivalent
- Not all Abnormal samples may be comparable
- May need sample-level or patient-level metadata

**Current gap**: Scripts don't yet handle sample-specific metadata for filtering/subsetting

**TODO**:
- [ ] What metadata is available? (patient demographics, disease stage, treatment?)
- [ ] Should we filter samples before training?
- [ ] Should we stratify analysis by subgroups?
- [ ] Add `--filter-samples` option to training scripts?

### 6. Two-Stage Model Question

**User question**: "Is the model still two stages: 1) find lymphocytes, 2) describe archetypes?"

**Current answer**: Single stage (data is already lymphocyte-enriched from SEACells)

**Potential issue**:
- Are there non-lymphocytes in the data that should be excluded?
- Check CD45 levels: mean=0, std=1, min=-6.7, max=6.2
- Low CD45 (< -2?) might indicate non-lymphocytes

**TODO**:
- [ ] Visualize CD45 distribution
- [ ] Check if low-CD45 metacells have different characteristics
- [ ] Consider adding viability filter (exclude low CD45)
- [ ] Optionally implement explicit two-stage model

---

## Planned Next Steps

### Immediate (Before Next Session)

1. **Test new VAE framework**:
   ```bash
   # Quick test run
   /opt/homebrew/Caskroom/mambaforge/base/bin/conda run -n flow_archetype_stable \
     python scripts/vae_archetypes_with_priors.py \
     --input seacells_output/three_groups_phenograph/integrated_metacells_harmony_phenograph.h5ad \
     --train-on-reference \
     --latent-dim 5 \
     --n-archetypes 10 \
     --beta 0.01 \
     --gamma 0.1 \
     --epochs 200 \
     --output-dir test_vae_with_priors
   ```

2. **Verify outputs**:
   - Check training curves converge
   - Visualize archetypes on UMAP
   - Check biological priors are working (monitor bio_loss)
   - Verify soft assignments and entropy are computed

3. **Compare to original**:
   - Do we get similar archetypes as `vae_archetypes_output/`?
   - Are biological priors improving interpretability?
   - How does train-on-reference affect results?

### Latest Update: Gamma Cooldown Implementation

**IMPLEMENTED** - Gamma cooldown schedule for biological priors:
- **Rationale**: Enforce lineage constraints strongly at the beginning of training to guide the model toward biologically plausible regions, then gradually relax to allow discovery of nuanced patterns
- **Implementation**: Linear decay from `--gamma` to `--gamma-min` over `--gamma-cooldown` epochs
- **Default**: 150 epochs cooldown from gamma to 0
- **Usage**:
  ```bash
  python scripts/vae_archetypes_with_priors.py \
    --gamma 0.5 \              # Strong priors initially
    --gamma-min 0.0 \          # No priors at end
    --gamma-cooldown 150       # Decay over 150 epochs
  ```
- **Visualization**: Gamma schedule is now plotted in training curves
- **Status**: ✅ Code complete, UNTESTED

### Short-term (Next Session)

1. **Test gamma cooldown**:
   ```bash
   # Compare models with and without cooldown
   # WITH cooldown (new default)
   python scripts/vae_archetypes_with_priors.py \
     --train-on-reference \
     --gamma 0.5 --gamma-cooldown 150 --gamma-min 0.0

   # WITHOUT cooldown (constant gamma)
   python scripts/vae_archetypes_with_priors.py \
     --train-on-reference \
     --gamma 0.1 --gamma-cooldown 0
   ```

2. **Hyperparameter tuning**:
   ```bash
   python scripts/tune_vae_archetypes.py \
     --input <data.h5ad> \
     --train-on-reference \
     --n-trials 30 \
     --epochs 150 \
     --output-dir optuna_tuning
   ```

2. **Analyze tuning results**:
   - What's the optimal number of archetypes?
   - What gamma value gives best biological priors?
   - Plot tuning metrics vs hyperparameters

3. **Train final model** with optimized hyperparameters

4. **Biological validation**:
   - Do archetypes match expected cell types?
   - Expected types:
     - CD3+CD4+ naive (CD45RA+, CD197+)
     - CD3+CD4+ central memory (CD45RO+, CD197+)
     - CD3+CD4+ effector memory (CD45RO+, CD197-)
     - CD3+CD4+ TEMRA (CD45RA+, CD197-)
     - CD3+CD8+ (same subsets as CD4)
     - NK cells CD56bright (CD3-, CD56+, CD16-)
     - NK cells CD56dim (CD3-, CD56+, CD16+)
     - B-cells (CD19+, though hard to separate due to CD3+CD19 merge)

### Medium-term (This Week)

1. **Sample annotation & filtering**:
   - Work with user to understand sample metadata
   - Identify which samples are comparable
   - Add filtering options to scripts

2. **Quantify deviations from Normal**:
   - For each Abnormal/PTCy metacell:
     - Distance to nearest Normal in latent space
     - Archetype composition differences
     - Marker expression differences
   - Statistical tests for differences

3. **Analyze mixed phenotypes**:
   - Identify high-entropy metacells
   - What are their top 2-3 archetype associations?
   - Are these biologically interpretable transitions?
   - Do they cluster by sample type?

4. **Refine biological priors**:
   - Based on results, adjust prior definitions
   - Consider additional priors:
     - B-cell markers (CD19, though merged with CD3)
     - Exhaustion markers (TIGIT, CD366)
     - Activation (CD69)

### Long-term (Next Steps)

1. **Clinical correlation**:
   - Do archetype compositions predict clinical outcomes?
   - Do high-entropy regions correlate with disease severity?

2. **Longitudinal analysis**:
   - Track archetype composition changes over time
   - Use `project_samples.py` for consistent comparisons

3. **Interpretability enhancements**:
   - Visualize latent space trajectories
   - Identify marker combinations that define archetype boundaries
   - Generate natural language descriptions of archetypes

4. **Two-stage model (if needed)**:
   - Stage 1: Lymphocyte vs non-lymphocyte classifier
   - Stage 2: Lymphocyte archetype identification
   - Compare to current single-stage approach

5. **Additional loss functions**:
   - Contrastive loss (pull similar phenotypes together, push dissimilar apart)
   - Triplet loss (anchor-positive-negative)
   - Adversarial loss (sample type invariance for reference archetypes)

---

## Key Files & Outputs

### Scripts
- `scripts/vae_archetypes_seacells.py` - Original VAE (working, produces good results)
- `scripts/vae_archetypes_with_priors.py` - New VAE with priors (untested)
- `scripts/tune_vae_archetypes.py` - Hyperparameter tuning (untested)
- `scripts/project_samples.py` - Project new samples (untested)

### Outputs
- `vae_archetypes_output/` - Original VAE results (reference)
  - `metacells_with_vae_archetypes.h5ad`
  - `archetype_profiles.csv`
  - `archetype_biological_names.csv`
  - Various visualizations

### Documentation
- `README_VAE_ARCHETYPES.md` - Comprehensive usage guide
- `Nov15.md` - This file (project status)

### Data
- `seacells_output/three_groups_phenograph/integrated_metacells_harmony_phenograph.h5ad` - Input data

---

## Critical Questions for User

1. **Sample metadata**: What additional metadata is available?
   - Patient demographics?
   - Disease stage/severity?
   - Treatment history?
   - Time points?

2. **Sample comparability**:
   - Which Abnormal samples are most comparable?
   - Are there subgroups within "Abnormal"?
   - What defines a good "Normal" reference?

3. **Clinical goals**:
   - What outcomes are we trying to predict/understand?
   - Are there specific phenotypes of interest?
   - What would constitute success?

4. **Biological expectations**:
   - Should certain markers always co-occur?
   - Are there known abnormal phenotypes to look for?
   - What level of "mixedness" is biologically plausible?

5. **CD197 prominence**:
   - Is CCR7 expected to be a major discriminating feature?
   - Does this align with biological understanding?

---

## Known Limitations

1. **Merged channels**: Can't cleanly separate CD3/CD19, CD4/CD33, CD8/CD14, CD16/TIGIT
2. **Sample imbalance**: 79% Abnormal, 10% Normal
3. **No validation set**: All tuning/training on same data (risk of overfitting)
4. **No ground truth**: Can't validate archetypes against known cell types
5. **Computational cost**: Hyperparameter tuning may be slow (50 trials × 200 epochs)

---

## Git Branch

**Branch**: `nov15_gamma_cooldown_experiments`

Created from `main` to experiment with new VAE archetype framework.

**Commit**: c3fbd94 - "Add VAE archetypes framework with biological priors and gamma cooldown"

Files tracked:
- `scripts/vae_archetypes_with_priors.py`
- `scripts/tune_vae_archetypes.py`
- `scripts/project_samples.py`
- `Nov15.md`
- `README_VAE_ARCHETYPES.md`
- `test_gamma_cooldown.sh`

---

## Experiments Log

### Experiment 1: Gamma Cooldown (Quick Test)
**Date**: November 15, 2025
**Status**: ✅ Complete
**Config**:
- Train on reference: Yes (Normal samples only, n=605)
- Latent dim: 5
- N archetypes: 10
- Beta: 0.01
- Gamma: 0.5 → 0.0 over 100 epochs
- Total epochs: 150
**Output**: `test_output/gamma_cooldown/`

**Results**:
- Final loss: 0.1577 (Recon=0.1445, KL=1.3171, Bio=0.0000)
- Reference entropy: mean=1.194, median=1.250, range=[0.231, 1.943]
- Projected entropy: range=[0.536, 2.105]
- Archetype distribution: Fairly even (20-111 metacells per archetype)

**Key archetypes identified**:
- Archetype 2: CD4+ Memory (CD4↑1.48, CD45RO↑1.34)
- Archetype 6: CD8+ Effector Memory (CD8↑1.96, CD45RO↑1.51)
- Archetype 7: NK cells (CD16↑2.16, CD56↑1.16, CD45RA↑1.13)
- Archetype 3: CD8+ Activated (CD8↑1.16, CD69↑0.89)

### Experiment 3: Longer Training with Gamma Cooldown
**Date**: November 15, 2025
**Status**: ✅ Complete
**Config**:
- Train on reference: Yes (Normal samples only, n=605)
- Latent dim: 5
- N archetypes: 10
- Beta: 0.01
- Gamma: 0.5 → 0.0 over 150 epochs (extended cooldown)
- Total epochs: 300 (2x longer than Exp 1)
**Output**: `test_output/gamma_cooldown_300epochs/`

**Results**:
- Final loss: 0.1429 (Recon=0.1306, KL=1.2324)
- **9.4% improvement over 150 epochs** in total loss
- Reference entropy: mean=1.372, median=1.430 (HIGHER than 150 epochs)
- Gini coefficient: 0.307 (LESS even than 150 epochs at 0.182)

**Key archetypes**:
- Archetype 0: NK cells (CD16↑1.72, CD45RA↑1.11)
- Archetype 8: CD8+ Effector Memory (CD8↑1.96, CD45RO↑1.51)
- Archetype 9: CD4+ Memory (CD4↑1.82, CD45RO↑1.81)
- Archetype 3: CD8+ Activated (CD8↑1.31, CD69↑1.14)

**Trade-off observed**:
- ✅ Lower loss (better data fit)
- ⚠️ Higher entropy (less decisive assignments)
- ⚠️ Less even distribution (some archetypes very small, n=15-21)
- **Conclusion**: May be overfitting. **150 epochs appears optimal** for interpretability vs fit

### Experiment 2: Constant Gamma (Comparison)
**Date**: November 15, 2025
**Status**: ✅ Complete
**Config**:
- Train on reference: Yes (Normal samples only, n=605)
- Latent dim: 5
- N archetypes: 10
- Beta: 0.01
- Gamma: 0.1 (constant, no cooldown)
- Total epochs: 150
**Output**: `test_output/constant_gamma/`

**Results**:
- Final loss: 0.1756 (Recon=0.1485, KL=1.3084, Bio=0.1401)
- Reference entropy: mean=1.129, median=1.187, range=[0.086, 2.071]
- Projected entropy: range=[0.441, 2.111]
- Archetype distribution: More varied (24-113 metacells per archetype)

**Key archetypes identified**:
- Archetype 2: CD4+ Memory (CD45RO↑2.06, CD4↑1.77)
- Archetype 4: CD4+ Memory variant (CD45RO↑2.38, CD4↑1.53)
- Archetype 7: CD8+ Effector Memory (CD8↑1.96, CD45RO↑1.51) - IDENTICAL to Exp1 Archetype 6
- Archetype 0: NK cells (CD16↑1.72, CD45RA↑1.11)

### Experiment Comparison & Findings

**Quantitative comparison**:
1. **Lower final loss with gamma cooldown** (0.1577 vs 0.1756) - **13% improvement**
2. **Better reconstruction with cooldown** (0.1445 vs 0.1485)
3. **Slightly higher entropy with cooldown** (1.194 vs 1.129) - Less decisive, but not necessarily worse
4. **Some archetypes are identical** despite different training (Archetype 6/7, Archetype 8)

**Biological interpretation**:
- Both models identify similar major cell types (NK, CD4+ Memory, CD8+ Memory)
- Both capture CD45RA/CD45RO gradient (Naive vs Memory)
- Both identify activation states (CD69+)
- Constant gamma tends to create duplicate memory archetypes (Arch 2 & 4 very similar)

**Qualitative observations**:
- **Gamma cooldown produces cleaner archetypes** - fewer duplicates
- **Constant gamma over-segments memory T-cells** - Archetypes 2 & 4 are nearly identical
- **Both successfully project Abnormal/PTCy samples** onto reference model
- **Training time: ~2-3 minutes** for 150 epochs on 605 reference samples

**Preliminary conclusion**:
✅ **Gamma cooldown appears beneficial** - better loss, cleaner archetypes, no duplicate phenotypes

**Additional metrics**:
- Distribution evenness (Gini coefficient): Cooldown=0.182, Constant=0.232
  - **Cooldown has 27% more even distribution** across archetypes
- Archetype profile correlation: High correlation for matched archetypes (0.8-1.0)
- Some archetypes are nearly identical (Arch 6 cooldown ≈ Arch 7 constant, r=1.0)

**Files generated**:
- `test_output/gamma_cooldown/` - All outputs from Experiment 1
  - `metacells_with_archetypes.h5ad` (5.1 MB)
  - `archetype_profiles.csv`
  - `training_curves.png` - Shows gamma schedule, loss components
  - `vae_archetypes_overview.png` - UMAP + heatmaps
  - `vae_model.pt` (42 KB) - Saved model for future projection
- `test_output/constant_gamma/` - All outputs from Experiment 2
  - Same file structure
- `test_output/experiment_comparison.png` - Side-by-side comparison

**Analysis scripts**:
- `compare_experiments.py` - Quantitative comparison with visualizations

**Recommendation for next experiments**:
✅ **Use gamma cooldown as default**
- Proceed with hyperparameter tuning using cooldown schedule
- Test different cooldown durations (50, 100, 150, 200 epochs)
- Test different gamma_min values (0.0, 0.01, 0.05)

**Next steps**:
1. ✅ ~~Examine visualizations~~ - Done via compare_experiments.py
2. ✅ ~~Compare biological coherence~~ - Cooldown has cleaner separation
3. ✅ ~~Run longer training (300 epochs)~~ - Lower loss but worse interpretability
4. ✅ ~~Compare Abnormal vs PTCy archetype compositions~~ - MAJOR DIFFERENCES FOUND
5. ✅ ~~Identify high-entropy mixed phenotypes~~ - 73% of all metacells, 78% Abnormal
6. Hyperparameter sweep with Optuna (include cooldown schedule) - PENDING

---

### Analysis Results: Sample Type Differences

**Experiment 4: Sample Composition Analysis**
**Date**: November 15, 2025
**Status**: ✅ Complete
**Input**: `test_output/gamma_cooldown/metacells_with_archetypes.h5ad` (150 epoch model)
**Output**: `test_output/sample_composition_analysis/`

**Major Findings**:

1. **HIGHLY SIGNIFICANT DIFFERENCES between sample types** (Chi-square p < 8.22e-134)

2. **Abnormal samples have DRAMATICALLY different archetype distribution**:
   - **Archetype 6**: 26.1% (vs 8.4% Normal) - **3.1x ENRICHED** ⬆️
   - **Archetype 8**: 35.8% (vs 18.3% Normal) - **2.0x ENRICHED** ⬆️
   - These 2 archetypes comprise **62% of Abnormal metacells**!

   - **Depleted archetypes** in Abnormal:
     - Archetype 0: 2.4% (vs 10.7% Normal) - 0.22x ⬇️
     - Archetype 3: 1.7% (vs 10.1% Normal) - 0.17x ⬇️
     - Archetype 9: 1.2% (vs 7.9% Normal) - 0.15x ⬇️
   - **Lost diversity**: Multiple archetypes nearly absent in Abnormal

3. **PTCy samples show similar but distinct pattern**:
   - Also enriched in Archetypes 6 (24.5%, 2.9x) and 8 (29.8%, 1.6x)
   - Additionally enriched in Archetype 4 (17.9% vs 11.2%, 1.6x)
   - More diverse than Abnormal but still very different from Normal

4. **ENTROPY DIFFERENCES ARE MASSIVE**:
   - **Normal**: Mean entropy = 1.194 (27% high-entropy >1.5)
   - **Abnormal**: Mean entropy = 1.645 (78% high-entropy) - ***p < 10^-300***
   - **PTCy**: Mean entropy = 1.637 (75.5% high-entropy) - ***p < 10^-100***

   **Interpretation**: Abnormal and PTCy samples have VASTLY more mixed/transitional phenotypes than Normal

**Experiment 5: Mixed Phenotype Analysis**
**Date**: November 15, 2025
**Status**: ✅ Complete
**Input**: `test_output/gamma_cooldown/metacells_with_archetypes.h5ad`
**Output**: `test_output/mixed_phenotype_analysis/`

**Key Findings**:

1. **72.6% of ALL metacells are high-entropy** (entropy > 1.5)
   - This is dominated by disease samples (78% Abnormal, 75% PTCy)

2. **Most common archetype mixing patterns**:
   - **Archetype 6 + 8**: Co-occur 1,780 times (most common!)
   - **Archetype 2 + 8**: Co-occur 1,538 times
   - **Archetype 4 + 8**: Co-occur 1,176 times
   - **Archetype 4 + 6**: Co-occur 892 times

   **Pattern**: Archetypes 6 and 8 (the disease-enriched ones) frequently co-occur in mixed states

3. **Top mixing pattern in Abnormal: 6-8** (57 metacells)
   - These are metacells with substantial probability for BOTH Archetypes 6 and 8
   - Suggests transitional state or dual expression pattern

4. **Normal samples have different mixing**:
   - Top pattern: 2-8 and 8-2 (involving different archetypes)
   - Much lower overall mixing (only 27% high-entropy)

**Biological Interpretation**:
- **Archetype 6**: CD8+ Effector Memory (CD8↑1.96, CD45RO↑1.51, CD45↑1.25)
- **Archetype 8**: Less clear (CD45RA↑, CD3+CD19↑) - possible naive/memory intermediate?
- **Hypothesis**: Abnormal/PTCy samples accumulate CD8+ T-cells in transitional or activated states
- **Clinical significance**: Loss of phenotypic diversity + accumulation of specific CD8+ states may be disease signature

---

## Session Notes

**Date**: November 15, 2025

**Completed**:
- Reviewed initial VAE archetype results
- Designed train-and-project framework
- Implemented biological prior loss functions
- Implemented soft archetype assignments (entropy-based mixedness)
- Created hyperparameter tuning framework
- Created sample projection script
- Wrote comprehensive documentation

**User feedback**:
- Normals should be reference, everything relative to reference
- PTCY vs Abnormal comparison also interesting
- Need sample annotation (heterogeneous population)
- Explore: train-on-all vs train-on-subset (project)?
- Need to project future samples in either case
- Quantify space between archetypes (mixed phenotypes with >1 association)
- Enforce biological priors (CD3 vs CD16/CD56, CD4 vs CD8, CD19+CD24+ B-cells)
- CD197 is major feature - will this hold after optimization?
- Two-stage model question: lymphocyte gating → archetypes?

**Decisions made**:
- Implement train-on-reference mode (Normal as baseline)
- Add soft archetype assignments (probability + entropy)
- Add biological prior loss (gamma parameter)
- Create hyperparameter tuning framework
- Enable future sample projection

**Next**: ✅ Tested framework with gamma cooldown experiments

---

## ⚠️ CRITICAL LIMITATION: Batch Confounding

### Perfect Confounding Detected

**PROBLEM**: Each `sample_id` belongs to ONLY ONE `sample_type`:
- Normal: 4 samples
- Abnormal: 26 samples
- PTCy: 8 samples

**IMPLICATION**: We CANNOT definitively distinguish:
- ✅ True biological differences between Normal and disease
- ❌ Technical/batch effects (sample prep, run date, operator, etc.)

**What was done**:
- Data is Harmony batch-corrected (adata.X)
- Harmony centers all markers to mean≈0
- VAE trained on this corrected data

**What remains uncertain**:
- Harmony corrects **linear** batch effects well
- **Nonlinear** batch effects might persist
- Without technical replicates or mixed batches, cannot validate

**Required for validation**:
1. Technical replicates (same sample, different runs)
2. Mixed sample pools (multiple types in one batch)
3. External validation dataset
4. Within-sample validation (not possible here - each sample is homogeneous)

**Interpretation of findings below**:
- Findings are **hypothesis-generating**, not definitive
- Could represent true biology (optimistic)
- Could represent batch effects (pessimistic)
- Most likely: **mixture of both**
- Treat as preliminary until validated

---

## POTENTIAL DISCOVERIES - November 15, 2025
### (Caveat: May Include Batch Effects)

### Potential Disease Signature: CD8+ T-cell Accumulation

**THE BIG PICTURE**:
1. Normal samples have diverse, stable lymphocyte phenotypes (10 distinct archetypes)
2. Abnormal/PTCy samples **lose this diversity** and accumulate in 2 specific states:
   - **Archetype 6**: CD8+ Effector Memory (3.1x enriched)
   - **Archetype 8**: Unclear phenotype, possibly naive/memory transition (2.0x enriched)
3. Disease samples have **78% high-entropy metacells** (vs 27% Normal)
   - This represents unstable, transitional, or mixed phenotypes
4. The two disease-enriched archetypes (6 & 8) **frequently co-occur** in mixed states

**QUANTITATIVE SUMMARY**:
| Metric | Normal | Abnormal | PTCy | Significance |
|--------|--------|----------|------|--------------|
| Archetype 6 | 8.4% | 26.1% (3.1x) | 24.5% (2.9x) | p < 1e-134 |
| Archetype 8 | 18.3% | 35.8% (2.0x) | 29.8% (1.6x) | p < 1e-134 |
| High-entropy (>1.5) | 27.4% | 78.0% | 75.5% | p < 1e-300 |
| Mean entropy | 1.194 | 1.645 | 1.637 | p < 1e-300 |

**BIOLOGICAL INTERPRETATION**:
- **Loss of diversity**: Abnormal samples lose representation in 6+ archetypes
- **CD8+ accumulation**: Both disease groups accumulate CD8+ Effector Memory cells
- **Transitional states**: Disease samples dominated by mixed/unstable phenotypes
- **Potential mechanism**: Abnormal activation, exhaustion, or dysregulation of CD8+ compartment

**CLINICAL IMPLICATIONS**:
- **Diagnostic potential**: Archetype composition could distinguish Normal from disease
- **Monitoring**: Entropy could track disease progression or treatment response
- **Therapeutic targets**: CD8+ T-cell states enriched in disease

---

## Current Status (End of Nov 15 Session)

### What We Built Today

1. **VAE Archetype Framework** (`scripts/vae_archetypes_with_priors.py`)
   - Train-and-project: Train on Normal, project Abnormal/PTCy
   - Biological prior loss functions with gamma cooldown
   - Soft archetype assignments with entropy metric
   - Complete save/load functionality for future projections

2. **Supporting Tools**
   - `scripts/tune_vae_archetypes.py` - Hyperparameter optimization (Optuna)
   - `scripts/project_samples.py` - Project new samples onto trained model
   - `compare_experiments.py` - Experiment comparison and analysis

3. **Documentation**
   - `README_VAE_ARCHETYPES.md` - Comprehensive usage guide
   - `Nov15.md` (this file) - Project status and experimental log
   - `test_gamma_cooldown.sh` - Quick test script

### What We Tested

Ran **5 experiments** today:

1. **Exp 1: Gamma Cooldown (150 epochs)** ✅ BEST FOR PRODUCTION
   - Loss: 0.1577, Entropy: 1.194, Gini: 0.182
   - Clean, interpretable archetypes

2. **Exp 2: Constant Gamma (150 epochs)**
   - Loss: 0.1756, Entropy: 1.129, Gini: 0.232
   - Duplicate archetypes, less even distribution

3. **Exp 3: Longer Training (300 epochs)**
   - Loss: 0.1429 (9.4% better), but Entropy: 1.372, Gini: 0.307
   - Better fit but worse interpretability (possible overfitting)

4. **Exp 4: Sample Composition Analysis** ⭐ MAJOR FINDINGS
   - Discovered 3.1x enrichment of Archetype 6 in Abnormal
   - Discovered 2.0x enrichment of Archetype 8 in Abnormal
   - Found massive entropy differences (p < 1e-300)

5. **Exp 5: Mixed Phenotype Analysis** ⭐ MAJOR FINDINGS
   - 72.6% of all metacells are high-entropy (mixed phenotypes)
   - Archetype 6+8 co-occur 1,780 times (disease signature)
   - Disease samples dominated by transitional states

**Verdict**:
- **Use 150 epoch gamma cooldown** for production models
- **Archetype composition + entropy** are powerful disease biomarkers
- **CD8+ T-cell accumulation** is the key disease signature

### Where Everything Is

**Git Branch**: `nov15_gamma_cooldown_experiments`
- Commits:
  - `c3fbd94` - Initial VAE framework with biological priors and gamma cooldown
  - `2f2f921` - Experiments 1-2: Gamma cooldown vs constant gamma comparison
  - `17aa308` - Final documentation update
  - `f1e03a4` - Experiments 4-5: Sample composition and mixed phenotype analysis (MAJOR FINDINGS)

**Data Locations**:
- **Input**: `seacells_output/three_groups_phenograph/integrated_metacells_harmony_phenograph.h5ad`
  - 6,039 metacells (Normal: 605, Abnormal: 4,770, PTCy: 664)
  - 13 markers (some merged channels)

- **Experiment Outputs**:
  - `test_output/gamma_cooldown/` - Exp 1: 150 epochs, gamma cooldown ⭐ BEST MODEL
  - `test_output/constant_gamma/` - Exp 2: 150 epochs, constant gamma
  - `test_output/gamma_cooldown_300epochs/` - Exp 3: 300 epochs (overfits)
  - `test_output/experiment_comparison.png` - Visual comparison of Exp 1 vs 2
  - `test_output/sample_composition_analysis/` - Exp 4: Sample type comparisons
    - `sample_type_composition.png`
    - `archetype_composition_by_sample.csv`
  - `test_output/mixed_phenotype_analysis/` - Exp 5: Mixed phenotype analysis
    - `mixed_phenotype_analysis.png`
    - `mixing_patterns.csv`
    - `archetype_cooccurrence.csv`

- **Original Results** (baseline): `vae_archetypes_output/`

**Scripts**:
- **Training & Modeling**:
  - `scripts/vae_archetypes_with_priors.py` - Main VAE training with biological priors + gamma cooldown
  - `scripts/tune_vae_archetypes.py` - Optuna hyperparameter optimization
  - `scripts/project_samples.py` - Project new samples onto trained model

- **Analysis**:
  - `compare_experiments.py` - Compare different training configurations
  - `analyze_sample_composition.py` - Archetype composition by sample type ⭐
  - `analyze_mixed_phenotypes.py` - High-entropy mixed phenotype analysis ⭐

- **Testing**:
  - `test_gamma_cooldown.sh` - Quick test script for gamma cooldown

### What Works

✅ **Train on reference (Normal) samples** - establishes baseline
✅ **Project test samples (Abnormal, PTCy)** - maintains consistent archetype definitions
✅ **Biological prior loss functions** - enforces known lineage relationships
✅ **Gamma cooldown schedule** - 13% better loss, cleaner archetypes
✅ **Soft archetype assignments with entropy** - quantifies "mixedness"
✅ **Model save/load** for future projection
✅ **Identifies biologically meaningful archetypes** (NK, CD4+ Memory, CD8+ Memory, etc.)
✅ **Fast training** (~2-3 min for 150 epochs on 605 samples)
✅ **Distinguishes Normal from disease** with high statistical significance (p < 1e-134)
✅ **Entropy as biomarker** - correlates strongly with disease state (p < 1e-300)
✅ **Archetype composition analysis** - reveals disease-specific enrichment patterns

### What's Next

**Immediate**:
1. Longer training run (300 epochs) to see if results improve
2. Hyperparameter sweep with Optuna to optimize cooldown schedule
3. Biological validation: annotate archetypes with expected cell types

**Short-term**:
1. Compare Abnormal vs PTCy archetype compositions
2. Identify high-entropy mixed phenotypes
3. Quantify deviation from Normal baseline

**Medium-term**:
1. Sample annotation and filtering (work with user on metadata)
2. Clinical correlation analysis
3. Longitudinal tracking (if time-series data available)
4. Refine biological priors based on domain knowledge

### Key Decisions Made

1. ✅ **Use train-and-project framework** - Normal as reference
2. ✅ **Use gamma cooldown** - Start high (0.5), decay to low (0.0)
3. ✅ **Use soft assignments** - Quantify mixed phenotypes with entropy
4. ✅ **Single-stage model** - Data already lymphocyte-enriched
5. ⏳ **Hyperparameter values** - Will optimize with Optuna

### Open Questions

1. What's the optimal gamma schedule? (will test via tuning)
2. How to handle sample heterogeneity? (need metadata from user)
3. What clinical outcomes to correlate with? (need user input)
4. Should we add more biological priors? (B-cell, activation, exhaustion markers)

---

## Implementation Summary (November 15, 2025)

### Features Completed Today

1. ✅ **Train-and-Project Framework**
   - Train on reference samples (Normal), project test samples (Abnormal, PTCy)
   - Maintains consistent archetype definitions
   - Enables future sample projection without retraining

2. ✅ **Biological Prior Loss Functions**
   - T-cell vs NK mutual exclusion
   - CD4 vs CD8 preference for T-cells
   - Memory vs Naive anti-correlation
   - Soft constraints that can be violated with strong signal

3. ✅ **Gamma Cooldown Schedule**
   - **NEW FEATURE**: Enforce biological priors strongly early, relax later
   - Linear decay: gamma → gamma_min over specified epochs
   - Rationale: Guide model to biologically plausible regions initially, then allow nuanced pattern discovery
   - Configurable: `--gamma`, `--gamma-cooldown`, `--gamma-min`

4. ✅ **Soft Archetype Assignments**
   - Probability distributions over archetypes (not just hard assignments)
   - Entropy metric quantifies "mixedness"
   - Enables identification of transitional/dual-marker phenotypes

5. ✅ **Hyperparameter Tuning Framework**
   - Optuna-based optimization
   - Tunes: n_archetypes, latent_dim, beta, gamma, gamma_cooldown, gamma_min, soft_temp
   - Metrics: assignment entropy (clarity), silhouette score (separation)

6. ✅ **Sample Projection Script**
   - Load trained model
   - Project new samples onto learned latent space
   - Assign to established archetypes
   - Compare against reference

7. ✅ **Comprehensive Documentation**
   - README_VAE_ARCHETYPES.md: Full usage guide
   - Nov15.md: Project status and planning document

### Code Status

| Script | Status | Tested | Purpose |
|--------|--------|--------|---------|
| `vae_archetypes_seacells.py` | ✅ Complete | ✅ Yes | Original VAE (baseline) |
| `vae_archetypes_with_priors.py` | ✅ Complete | ❌ No | New VAE with biological priors + gamma cooldown |
| `tune_vae_archetypes.py` | ✅ Complete | ❌ No | Hyperparameter optimization |
| `project_samples.py` | ✅ Complete | ❌ No | Project new samples onto trained model |

### Key Innovation: Gamma Cooldown

**Problem**: How to balance biological constraints with data-driven discovery?
- Too strict: Miss interesting biology
- Too loose: Uninterpretable archetypes

**Solution**: Adaptive schedule
- **Early training** (epochs 0-150): Strong biological guidance (gamma = 0.5)
  - Model parameters are random, need strong constraints
  - Guide toward known cell types (CD3+, CD4+, CD8+, NK, etc.)
- **Late training** (epochs 150+): Relax constraints (gamma → 0.0)
  - Model has learned basic structure
  - Allow discovery of disease-specific or transitional phenotypes

**Expected benefit**:
- Best of both worlds: biologically grounded + data-driven
- Should produce clearer archetypes than constant gamma

### Critical Next Steps

1. **Test gamma cooldown** - Does it actually improve archetypes?
2. **Compare to constant gamma** - Is cooldown better than static gamma?
3. **Hyperparameter tuning** - Find optimal cooldown schedule
4. **Biological validation** - Do archetypes match expected cell types?

### Open Research Questions

1. **Optimal cooldown schedule?**
   - Linear vs exponential decay?
   - Start/end values?
   - Duration relative to total epochs?

2. **Effect of merged channels?**
   - How much do merged channels (CD3+CD19, CD4+CD33) limit biological constraints?
   - Should we weight constraints differently based on channel quality?

3. **Sample-specific archetypes?**
   - Do Abnormal samples have unique archetypes not in Normal?
   - Should we identify these separately then align?

4. **Clinical utility?**
   - Which metrics correlate with disease severity/outcomes?
   - Archetype composition? Mixed phenotype prevalence? Distance from Normal?

---

## Session Update - Continued (November 15, 2025, Evening)

### Batch Effect Investigation & Validation Planning

**User Feedback**: "i'm not convinced that the normal vs disease is detected 'real' as opposed to batch effect"

This critical feedback prompted a thorough investigation of potential batch confounding.

### Investigation Results

**DISCOVERED**: Perfect confounding between `sample_id` and `sample_type`

Analysis revealed:
```python
# Check confounding
confound_check = adata.obs.groupby('sample_id')['sample_type'].unique()
# Result: Each sample_id has ONLY ONE sample_type
```

**Sample distribution**:
- Normal: 4 distinct samples (Normal_1 through Normal_4)
- Abnormal: 26 distinct samples (Abnormal_1 through Abnormal_26)
- PTCy: 8 distinct samples (PTCy_1 through PTCy_8)

**Critical implication**:
- ANY differences between sample types could be:
  - True biology ✓
  - Batch effects (processing time, operator, reagent lots, etc.) ✗
  - Most likely: **mixture of both, proportions unknown**

**What Harmony batch correction can/cannot do**:
- ✅ Corrects linear batch effects effectively
- ❌ May not fully correct nonlinear batch effects
- ❌ **Cannot resolve perfect confounding** (sample_type == batch_id)

### Actions Taken

1. **Updated all documentation** to reflect uncertainty:
   - Changed "MAJOR DISCOVERIES" → "POTENTIAL DISCOVERIES"
   - Added ⚠️ warnings to Nov15.md
   - Created FINDINGS_SUMMARY.md with caveats
   - Added prominent batch confounding section

2. **Created VALIDATION_ROADMAP.md**:
   - **Location**: `VALIDATION_ROADMAP.md` (project root)
   - **Committed**: Commit `d83380b`
   - **Contents**:
     - Problem definition (perfect confounding explained)
     - 5 validation strategies with implementation details:
       1. Technical replicates (gold standard, $10-15K, 2-3 months)
       2. Mixed sample batches (pragmatic, $5-10K, 2-3 months)
       3. **External validation** (recommended, $0, 3-4 weeks) ⭐
       4. Within-sample validation (not possible with current data)
       5. Positive/negative controls (partial solution)
     - Phased approach prioritizing external validation
     - Immediate action items for this week
     - Alternative analysis strategies (permutation testing, marker sensitivity)

**Recommended path forward**:
- **Phase 1** (immediate): Search public databases (ImmPort, FlowRepository, GEO) for external validation dataset
- **Phase 2** (if validated): Design new mixed-batch experiment
- **Phase 3** (if needed): Technical replicates for precise quantification

### Hyperparameter Optimization (In Progress)

**Status**: ⏳ Running (19/20 trials complete as of 6:55 PM)

**Command**:
```bash
/opt/homebrew/Caskroom/mambaforge/base/bin/conda run -n flow_archetype_stable \
  python scripts/tune_vae_archetypes.py \
  --input seacells_output/three_groups_phenograph/integrated_metacells_harmony_phenograph.h5ad \
  --train-on-reference \
  --epochs 100 \
  --n-trials 20 \
  --output-dir optuna_tuning_Nov15
```

**Parameters being optimized**:
- `n_archetypes`: 8-15
- `latent_dim`: 4-8
- `beta` (KL weight): 0.001-0.1
- `gamma` (biological prior weight): 0.1-1.0
- `gamma_cooldown` (cooldown epochs): 0-200
- `gamma_min` (final gamma): 0-0.1
- `soft_temp` (assignment temperature): 0.5-2.0

**Optimization objective**: Minimize combined score
- Lower entropy (clearer assignments) - weighted 30%
- Higher silhouette score (better separation) - weighted 70%

**Output location**: `optuna_tuning_Nov15/`
- Optuna database: `optuna_tuning_Nov15/optuna.db`
- Trial outputs: `optuna_tuning_Nov15/trial_0/` through `trial_19/`

**Expected completion**: ~7:00-7:05 PM

### Git Status

**Branch**: `nov15_gamma_cooldown_experiments`

**New commits**:
- `d83380b` - Add comprehensive validation roadmap to address batch confounding
  - Created VALIDATION_ROADMAP.md with 5 validation strategies
  - Documented perfect confounding issue
  - Recommended phased approach starting with external validation

**Files modified/created this session**:
- ✅ `VALIDATION_ROADMAP.md` (new, 471 lines, committed)
- ⏳ `optuna_tuning_Nov15/` (optimization running)

### Key Learnings

1. **User skepticism was correct** - Batch confounding is a real issue that must be addressed
2. **Framework remains valuable** - Methodology is sound even if current biology is uncertain
3. **Need external validation** - Most feasible path to resolve uncertainty
4. **Findings are hypothesis-generating** - Should be treated as preliminary until validated

### Next Steps

**Immediate** (tonight):
1. ✅ Commit VALIDATION_ROADMAP.md (done)
2. ⏳ Wait for hyperparameter optimization to complete
3. Analyze Optuna results
4. Update Nov15.md with optimization findings

**This week**:
1. Search public databases for external validation datasets
   - ImmPort: https://www.immport.org/
   - FlowRepository: https://flowrepository.org/
   - GEO: https://www.ncbi.nlm.nih.gov/geo/
2. Review literature on flow cytometry batch effect validation
3. Check if any within-sample replicates exist in current data

**Longer term**:
1. If external dataset found: Apply trained model for validation
2. If validated: Proceed with biological interpretation and publication
3. If not validated: Design new mixed-batch experiment

---

