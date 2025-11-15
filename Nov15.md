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

**Next**: Test new framework, run hyperparameter tuning, validate archetypes

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
