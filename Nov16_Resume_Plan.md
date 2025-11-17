# Project Status & Resumption Plan (Nov 16)

This document summarizes the critical findings and the current state of the project. It provides a clear path to resume the work from where we left off.

## 1. Summary of Key Findings

Our recent work has uncovered and addressed a fundamental issue with the input data, leading to a significant overhaul of the analysis pipeline.

### 1.1. The Marker Panel Discrepancy
- **Problem**: We discovered that the `integrated_metacells_harmony_phenograph.h5ad` file, which has been the input for all VAE modeling so far, was generated using an **incorrect and incomplete 13-marker panel**.
- **Diagnosis**: The original SEACells generation script (`seacells_three_groups_phenograph.py`) contained a hard-coded, outdated marker list. This caused it to select only 13 of the available markers from the raw data, silently dropping several important channels.

### 1.2. Establishing Ground Truth
- **Action**: I created and ran an inspection script (`inspect_fcs_v2.py`) on the raw `.fcs` files located in `data/normalsgated/`.
- **Finding**: The raw data contains a **definitive 17-marker panel**. Your hypothesis about a "missing" 18th marker was also correct; a 'Viability' channel exists and was correctly excluded from the final list.
- **The Ground Truth 17-Marker Panel is**:
  ```python
  ['CD62L', 'CD152 R718', 'CD45', 'CD45RO', 'CD279+CD24', 'CD95', 
   'CD16+TIGIT', 'CD8+CD14', 'CD56', 'CD45RA', 'CD366', 'CD69', 
   'CD25 BB515', 'CD34+CD223', 'CD197', 'CD3+CD19', 'CD4+CD33']
  ```

### 1.3. The Disambiguation VAE Model
- **Action**: Based on your suggestion, I developed and successfully tested a new VAE architecture (`VAEDisambiguation`) that explicitly handles the merged marker channels.
- **Capability**: This model learns to de-convolve the 6 merged channels into their 12 constituent markers, resulting in a much more interpretable internal representation of 23 unique markers. This is a major methodological improvement.

---

## 2. Current State of the Codebase

All the following new scripts have been created and are saved on the current git branch (`disambiguation-model-dev`).

### 2.1. New Corrected Pipeline Scripts
A full, end-to-end pipeline has been created to process the data with the correct marker panel and run the new model.

- **`scripts/pipeline_step1_seacells.py`**:
  - **Purpose**: The new data generation script. It replaces the old, incorrect script.
  - **Function**: Reads raw `.fcs` files, runs SEACells using the **correct 17 markers**, combines metacells, and performs Harmony integration.
  - **Output**: A new, corrected `integrated_metacells_17_markers.h5ad` file.

- **`scripts/pipeline_step2_train_model.py`**:
  - **Purpose**: Trains the new `VAEDisambiguation` model.
  - **Function**: Takes the 17-marker `.h5ad` file from Step 1 as input and trains the advanced VAE model, which is architected to handle the 17 -> 23 marker disambiguation.

- **`scripts/pipeline_step3_tune_model.py`**:
  - **Purpose**: Performs a broad hyperparameter search for the new model.
  - **Function**: Uses Optuna to find the best model configurations based on archetype uniqueness and separation, in line with your latest guidance.

### 2.2. New Documentation
- **`PIPELINE_README.md`**: A new, detailed document that explains the entire 3-step workflow and provides the exact commands needed to run each script in the corrected pipeline.

---

## 3. Immediate Next Step: Regenerate the Data

**The project is currently paused at a critical juncture.** We have corrected the methodology, but we have not yet applied it to the data.

The immediate and necessary next step is to **execute Step 1 of the new pipeline**. This will regenerate the base `h5ad` file using the full 17-marker panel.

**This is a long-running, computationally intensive process.**

**Command to Resume:**
To proceed, you must run the following command. I was about to run this before you paused the work.

```bash
/opt/homebrew/Caskroom/mambaforge/base/bin/conda run -n flow_archetype_stable python scripts/pipeline_step1_seacells.py \
  --output-dir seacells_pipeline_output_17_markers \
  --n-parallel 6 \
  --compression-ratio 40.0
```

---

## 4. Path Forward After Data Regeneration

Once the command above completes successfully, the project can proceed as follows:

1.  **Train and Verify the New Model**: Run `scripts/pipeline_step2_train_model.py` on the newly generated `h5ad` file to train a single instance of the `VAEDisambiguation` model and analyze its disambiguated profiles.

2.  **Launch the Broad Hyperparameter Sweep**: Run `scripts/pipeline_step3_tune_model.py` to find the optimal model configurations. This will provide the foundation for the deep-dive analysis into mixed phenotypes that you requested.

3.  **Analyze Results**: Analyze the tuning results to select the most promising models for downstream biological interpretation.

The project is now on a solid foundation. The next step of regenerating the data is the key to unlocking the full potential of the new disambiguation model.
