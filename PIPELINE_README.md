# Corrected Analysis Pipeline (Nov 16)

This document outlines the corrected, end-to-end pipeline for generating metacells, running Harmony integration, and training the VAE Disambiguation model. This pipeline uses the ground-truth 17-marker panel extracted from the raw `.fcs` files.

## Overview

The workflow consists of three main steps, executed by three corresponding scripts:

1.  **SEACell Generation & Integration**: Reads raw `.fcs` files, generates metacells using the 17-marker panel, and runs Harmony batch correction.
2.  **VAE Model Training**: Trains the `VAEDisambiguation` model on the output from Step 1.
3.  **Hyperparameter Tuning**: Runs a broad Optuna search to find optimal hyperparameters for the `VAEDisambiguation` model.

---

## Step 1: Generate Analysis-Ready Metacells

This step uses the corrected 17-marker panel to generate a new `h5ad` file from the raw data. This script is **resumable** and **robust**. It will skip already processed samples and will strictly enforce that all 17 markers are present in a sample before processing, skipping invalid files.

**Script:** `scripts/pipeline_step1_seacells_resumable_v5.py`

**What it does:**
- Strictly validates that all 17 markers are present in each `.fcs` file.
- Runs SEACells on each valid sample individually.
- Combines the resulting metacells from all samples.
- Applies Harmony batch correction on the `sample_id` column.
- Performs Phenograph clustering on the integrated data.
- Saves the final, analysis-ready object.

**Command to run:**
```bash
/opt/homebrew/Caskroom/mambaforge/base/bin/conda run -n flow_archetype_stable python scripts/pipeline_step1_seacells_resumable_v3.py \
  --output-dir seacells_pipeline_output_17_markers \
  --n-parallel 6 \
  --compression-ratio 40.0 \
  --log-resources
```

**Key Outputs:**
- `seacells_pipeline_output_17_markers/integrated_metacells_17_markers.h5ad`: The main data file for the next step.

---

## Step 2: Train the VAE Disambiguation Model

This step trains the advanced VAE model on the newly generated 17-marker data.

**Script:** `scripts/pipeline_step2_train_model.py`

**What it does:**
- Loads the `integrated_metacells_17_markers.h5ad` file.
- Uses a VAE model with a Disambiguation Layer that maps the 17 input markers to 23 un-merged markers internally.
- Applies biological priors on the internal, un-merged representation.
- Saves the trained model and the annotated `.h5ad` file with archetype assignments.

**Command to run (for a test run):**
```bash
/opt/homebrew/Caskroom/mambaforge/base/bin/conda run -n flow_archetype_stable python scripts/pipeline_step2_train_model.py \
  --input seacells_pipeline_output_17_markers/integrated_metacells_17_markers.h5ad \
  --output-dir vae_pipeline_output_17_markers
```

**Key Outputs:**
- `vae_pipeline_output_17_markers/metacells_with_archetypes.h5ad`: Data with archetype assignments.
- `vae_pipeline_output_17_markers/archetype_profiles_disambiguated.csv`: The interpretable, un-merged marker profiles for each archetype.
- `vae_pipeline_output_17_markers/vae_model.pt`: The saved, trained model.

---

## Step 3: Tune Hyperparameters (Optional)

This step runs a broad hyperparameter search to find the best-performing models according to the new objective (uniqueness and separation).

**Script:** `scripts/pipeline_step3_tune_model.py`

**What it does:**
- Uses Optuna to systematically test a wide range of hyperparameters.
- For each trial, it calls the Step 2 script to train a model.
- It scores each model based on a combination of archetype uniqueness (low correlation) and cluster separation (silhouette score).

**Command to run:**
```bash
/opt/homebrew/Caskroom/mambaforge/base/bin/conda run -n flow_archetype_stable python scripts/pipeline_step3_tune_model.py \
  --input seacells_pipeline_output_17_markers/integrated_metacells_17_markers.h5ad \
  --output-dir vae_tuning_17_markers \
  --n-trials 100
```
**Note:** This is a long-running process and should be executed when significant computational time is available.