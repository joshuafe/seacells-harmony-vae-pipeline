# Validation Roadmap: Resolving Batch Confounding

**Date**: November 15, 2025
**Issue**: Perfect confounding between `sample_type` and `sample_id`
**Status**: Cannot definitively distinguish biological differences from batch effects

---

## The Problem

### Current Data Structure

| Sample Type | Number of Samples | Sample IDs |
|-------------|------------------|------------|
| Normal | 4 | Normal_1, Normal_2, Normal_3, Normal_4 |
| Abnormal | 26 | Abnormal_1 through Abnormal_26 |
| PTCy | 8 | PTCy_1 through PTCy_8 |

**Key Issue**: Each `sample_id` belongs to ONLY ONE `sample_type`

### What This Means

When we observe differences between Normal and Abnormal:
- ✅ Could be **true biology**: Disease changes lymphocyte phenotypes
- ❌ Could be **batch effects**: Different processing, operators, dates, reagent lots
- 🤷 **Most likely**: Mixture of both, but proportions unknown

Even with Harmony batch correction:
- Linear batch effects → corrected
- Nonlinear batch effects → may persist
- Sample-level confounding → **cannot be resolved computationally**

---

## Validation Strategies

### Strategy 1: Technical Replicates (Gold Standard)

**Concept**: Process the **same biological sample** multiple times in different batches

**Implementation**:
```
Sample: Patient X (Abnormal)
├── Replicate A: Processed on Monday, Batch 1
├── Replicate B: Processed on Tuesday, Batch 2
└── Replicate C: Processed on Wednesday, Batch 3
```

**Analysis**:
1. Compare replicates (A vs B vs C)
   - Variation within replicates = **technical/batch variation**
2. Compare across sample types
   - If variation **within** replicates < variation **between** sample types:
     - ✅ True biological differences dominate
   - If variation **within** replicates ≥ variation **between** sample types:
     - ❌ Batch effects dominate

**Required**:
- Take 3-5 samples from EACH group (Normal, Abnormal, PTCy)
- Process each sample 3 times in different batches
- Total: ~45 runs (15 samples × 3 replicates)

**Advantages**:
- ✅ Gold standard for batch effect assessment
- ✅ Quantifies technical variation precisely
- ✅ Can identify which markers are batch-sensitive

**Limitations**:
- Expensive (45 runs)
- Time-consuming
- Requires frozen aliquots of original samples

---

### Strategy 2: Mixed Sample Batches (Pragmatic)

**Concept**: Include samples from **all groups** in **each batch**

**Implementation**:
```
Batch 1:
├── Normal_1A
├── Abnormal_1A
└── PTCy_1A

Batch 2:
├── Normal_1B
├── Abnormal_1B
└── PTCy_1B

... (repeat for N batches)
```

**Analysis**:
1. Model: `outcome ~ sample_type + batch`
2. If `sample_type` significant after accounting for `batch`:
   - ✅ True biological differences
3. If `sample_type` not significant:
   - ❌ Differences were batch effects

**Required**:
- Redesign experiment from scratch
- Include all 3 sample types in every batch
- Minimum: 10 batches × 3 sample types = 30 samples

**Advantages**:
- ✅ Statistically robust
- ✅ Can separate batch from biology using regression
- ✅ Standard experimental design

**Limitations**:
- Requires NEW experiment
- Cannot be applied to existing data

---

### Strategy 3: External Validation Dataset (Feasible)

**Concept**: Test findings on **independent dataset** from different lab/center

**Implementation**:
1. Find published flow cytometry dataset with:
   - Normal and disease samples
   - Similar markers (CD3, CD4, CD8, CD45RA, CD45RO, etc.)
   - Processed by different lab/protocol
2. Apply our trained VAE model to this dataset
3. Check if same patterns emerge:
   - CD8+ accumulation in disease?
   - High entropy in disease?
   - Similar archetype enrichment?

**Analysis**:
- If **same patterns** in external data:
  - ✅ Likely real biology (replicated across batches/labs)
- If **different patterns**:
  - ❌ Original findings were batch-specific

**Required**:
- Identify suitable public dataset (ImmPort, FlowRepository, GEO)
- Align markers (may need imputation if panels differ)
- Apply projection script: `project_samples.py`

**Advantages**:
- ✅ Can use existing public data
- ✅ No new experiments needed
- ✅ Most cost-effective

**Limitations**:
- Marker panels may not match exactly
- Disease types may differ (limits direct comparison)
- Preprocessing differences

**Recommended Databases**:
- **ImmPort**: https://www.immport.org/
- **FlowRepository**: https://flowrepository.org/
- **GEO**: https://www.ncbi.nlm.nih.gov/geo/

---

### Strategy 4: Within-Sample Validation (Not Possible Here)

**Concept**: Compare technical replicates within the **same sample**

**Why not possible**:
- Our data: 1 metacell set per sample_id
- No within-sample replicates available
- Each sample already aggregated via SEACells

**If data were available**:
- Take same blood draw
- Split into technical replicates
- Process separately
- Compare variance within vs between samples

---

### Strategy 5: Positive/Negative Controls (Partial Solution)

**Concept**: Include samples with **known** phenotypes as controls

**Implementation**:
```
Batch 1:
├── Sample A (unknown)
├── Sample B (unknown)
├── Control: Healthy donor (expected: diverse phenotypes)
└── Control: Activated T-cells (expected: CD8+ high)
```

**Analysis**:
- Check if controls show expected phenotypes
- If controls behave as expected across batches:
  - ✅ Batch effects are minimal
- If controls vary by batch:
  - ❌ Significant batch effects present

**Required**:
- Obtain well-characterized control samples
- Include in every batch
- Process identically to unknowns

**Advantages**:
- ✅ Provides batch quality assessment
- ✅ Can detect gross batch effects

**Limitations**:
- Doesn't prove unknowns are batch-free
- Requires standardized controls
- Controls may not match disease phenotypes

---

## Recommended Approach (Prioritized)

### Phase 1: External Validation (Immediate, Low Cost)

**Action Items**:
1. **Search for external datasets** (1-2 weeks)
   - Keywords: "flow cytometry", "lymphocyte", "T-cell", "disease", "Normal"
   - Check ImmPort, FlowRepository, published papers
   - Look for datasets with ≥ 5 markers overlapping with ours

2. **Apply our model to external data** (1 week)
   - Use `project_samples.py` with our trained model
   - Compare archetype distributions
   - Check if CD8+ accumulation replicates

3. **Document results** (1 week)
   - If validated: ✅ Strong evidence for biology
   - If not validated: ❌ Likely batch effects
   - If mixed: 🤷 Need more data

**Estimated time**: 3-4 weeks
**Cost**: Computational only (~$0)

---

### Phase 2: Mixed-Batch Experiment (If Phase 1 Promising)

**Action Items**:
1. **Design new experiment** (2 weeks)
   - Power analysis: How many samples needed?
   - Randomization: How to assign samples to batches?
   - Blinding: Can operators be blinded to sample type?

2. **Execute experiment** (4-6 weeks)
   - Process 30-50 samples across 10 batches
   - Each batch contains all 3 sample types
   - Record all batch metadata (date, operator, reagent lot)

3. **Analysis** (2 weeks)
   - Apply same VAE framework
   - Statistical model: `outcome ~ sample_type + batch`
   - Quantify effect sizes: biology vs batch

**Estimated time**: 2-3 months
**Cost**: ~$5,000-10,000 (depends on sample acquisition, processing)

---

### Phase 3: Technical Replicates (If Definitive Answer Needed)

**Only if**:
- Phase 1 shows promising signal
- Phase 2 shows batch effects persist
- Need to quantify exact technical variation

**Action Items**:
1. Select representative samples (n=15)
2. Process each 3 times in different batches
3. Variance component analysis

**Estimated time**: 2-3 months
**Cost**: ~$10,000-15,000

---

## Immediate Actions (This Week)

### 1. Literature Search

**Goal**: Find papers using similar approaches to identify batch effects

**Queries**:
- "flow cytometry batch correction validation"
- "lymphocyte phenotype batch effects"
- "Harmony batch correction validation"
- "flow cytometry technical replicates"

**Expected**: 10-20 relevant papers with validation strategies

---

### 2. Public Data Search

**Goal**: Identify at least 1 external dataset for validation

**Databases to check**:

1. **ImmPort** (https://www.immport.org/)
   - Search: "lymphocyte" or "T cell"
   - Filter: Flow cytometry
   - Look for: Normal + disease samples

2. **FlowRepository** (https://flowrepository.org/)
   - Browse recent uploads
   - Filter: Human, lymphocyte
   - Download metadata first

3. **GEO** (https://www.ncbi.nlm.nih.gov/geo/)
   - Search: "flow cytometry" AND "lymphocyte" AND "disease"
   - Filter: Homo sapiens
   - Check marker overlap

**Expected**: 2-5 candidate datasets

---

### 3. Within-Sample Variability Check

**Goal**: Check if current data has ANY within-sample replicates we missed

**Action**:
```python
import scanpy as sc
adata = sc.read_h5ad('seacells_output/.../integrated_metacells_harmony_phenograph.h5ad')

# Check metacell IDs for patterns
print(adata.obs.index[:50])  # Look for duplicate patterns
print(adata.obs['sample_id'].value_counts())  # Check metacells per sample

# Check if any samples were run multiple times
sample_summary = adata.obs.groupby('sample_id').size()
print(sample_summary)
```

**If found**: Can use for within-sample validation!
**If not found**: Proceed with external validation

---

## Alternative Analysis Strategies

While waiting for validation data, we can strengthen our current analysis:

### 1. Marker-Specific Sensitivity Analysis

**Goal**: Identify which markers drive the differences

**Method**:
```python
# For each marker, train VAE excluding that marker
# Check if Normal vs disease separation persists

for excluded_marker in markers:
    train_vae(data_without_marker)
    check_separation()
    # If separation disappears: that marker was critical
```

**Interpretation**:
- If separation depends on 1-2 markers: **more likely batch effect**
- If separation consistent across many markers: **more likely biology**

---

### 2. Permutation Testing

**Goal**: Test if separation is stronger than expected by chance

**Method**:
```python
# Shuffle sample_type labels
for i in range(1000):
    shuffle_labels()
    train_vae()
    compute_separation()

# Compare real separation to null distribution
p_value = (n_null >= real_separation) / 1000
```

**Interpretation**:
- If p < 0.001: Separation is real (but could still be batch!)
- If p > 0.05: Separation not significant

---

### 3. Marker Correlation Analysis

**Goal**: Check if markers correlate as expected biologically

**Method**:
```python
# Known biological correlations:
# - CD4 and CD8 should anti-correlate
# - CD45RA and CD45RO should anti-correlate
# - CD3 and CD56 should anti-correlate (T vs NK)

# If these hold: biology likely preserved
# If broken: batch effects may have distorted relationships
```

---

## Summary & Recommendations

### Immediate Priority: **External Validation**

1. **This week**: Search public databases for validation dataset
2. **Next week**: Apply our model to external data
3. **Week 3**: Document results and decide next steps

### If External Validation Succeeds:
- ✅ Findings are likely real biology
- Proceed with biological interpretation
- Prepare manuscript

### If External Validation Fails:
- ❌ Findings are likely batch effects
- Design new mixed-batch experiment
- Use current work as methodology pilot

### If External Data Not Found:
- Plan mixed-batch experiment (Phase 2)
- Consider technical replicate subset (5 samples × 3 reps)
- Continue with permutation/sensitivity analyses

---

## Resources Needed

### Computational
- [x] Trained VAE model (✅ have it: `test_output/gamma_cooldown/vae_model.pt`)
- [x] Projection script (✅ have it: `scripts/project_samples.py`)
- [ ] External dataset (searching)

### Experimental (if needed)
- [ ] Budget: $5,000-15,000 for validation experiments
- [ ] Samples: Fresh or frozen aliquots
- [ ] Time: 2-3 months

### Personnel
- [ ] Bioinformatician: Data download and preprocessing
- [ ] Domain expert: Biological interpretation
- [ ] Lab technician: If new experiments needed

---

## Conclusion

**Bottom line**: The framework and methods are sound, but **biology is uncertain due to batch confounding**.

**Path forward**:
1. **Validate externally** (easiest, cheapest, fastest)
2. **If validated**: Strong evidence for biology → publish
3. **If not validated**: Likely batch effects → redesign experiment

**Timeline**:
- External validation: 3-4 weeks
- New experiment (if needed): 2-3 months
- Technical replicates (if needed): 2-3 months

**The framework remains valuable** regardless of validation outcome:
- Methodology is novel and rigorous
- Can be applied to better-designed data
- Gamma cooldown innovation is generalizable

---

**Next Steps**: See `NEXT_STEPS_VALIDATION.md` for detailed action items and timeline.
