# Comprehensive Metacell Analysis Pipeline Guide

This document serves as a comprehensive guide to the entire metacell analysis pipeline, from initial sample processing and metacell generation to VAE model training, archetype analysis, clustering, and various diagnostics and visualizations. It aims to provide all necessary information to understand, run, and update the pipeline without additional external knowledge.

## 0. Overview

The pipeline is designed to analyze flow cytometry data, specifically focusing on bone marrow aspirates, to identify cellular archetypes and clusters, and to compare different biological cohorts (e.g., Normal, PTCy, Rux). It leverages SEACells for metacell generation, a custom Variational Autoencoder (VAE) with a disambiguation layer for archetype discovery, and various analytical scripts for in-depth interpretation.

The workflow is broadly divided into:
1.  **Data Preparation & Metacell Genesis (Step 1):** Processing raw FCS files into integrated metacells.
2.  **VAE Model Training & Archetype Generation (Step 2):** Training the VAE and deriving cellular archetypes.
3.  **Analysis & Visualization:** A suite of scripts for interpreting VAE archetypes, performing clustering analysis, and generating publication-ready plots.

All scripts are located in the `scripts/` directory, and output files are typically directed to specific analysis directories (e.g., `d28_vae_analysis/`, `d28_clustering_analysis/`).

## 1. Data Preparation and Metacell Genesis (Step 1)

This initial step processes raw `.fcs` files, generates metacells, and integrates them using Harmony batch correction.

*   **Script:** `scripts/pipeline_step1_seacells_resumable_v6.py`
*   **Purpose:** To convert raw flow cytometry data into a harmonized AnnData object of metacells, suitable for downstream VAE training and clustering. It ensures data quality by validating marker presence and handles batch effects.

### 1.1. Input Files

*   **Raw FCS files:** Located in `data/abnormalfcs_gated/` and `data/normalsgated/`. These are the primary input.
*   `sample_annotations.csv`: Located in the project root. This CSV contains metadata for each sample, crucial for cohort definition and batch correction. It includes columns like `sample_id`, `sample_type`, `gvh_prophylaxis`, `day_relative_to_transplant`, etc.

### 1.2. Key Parameters (in `pipeline_step1_seacells_resumable_v6.py`)

*   `--output-dir`: Directory to save all outputs from this step.
    *   **Default:** `seacells_pipeline_output_17_markers`
    *   **Impact:** Changing this will create a new output directory for the processed data.
*   `--n-parallel`: Number of parallel processes to use for SEACells generation.
    *   **Default:** `6`
    *   **Impact:** Affects computation time. Adjust based on available CPU cores.
*   `--compression-ratio`: Target compression ratio for SEACells.
    *   **Default:** `40.0`
    *   **Impact:** Higher values mean fewer, larger metacells; lower values mean more, smaller metacells. Affects granularity of analysis.
*   `--log-resources`: Flag to enable logging of resource usage.
    *   **Default:** `False`
    *   **Impact:** Useful for monitoring performance.

### 1.3. Outputs

*   `seacells_pipeline_output_17_markers/integrated_metacells_17_markers.h5ad`: The primary output. An AnnData object containing:
    *   `adata.X`: Scaled expression data for 17 markers.
    *   `adata.obs`: Observation metadata, including `sample_id`, `cohort`, and Harmony batch correction information.
    *   `adata.obsm['X_harmony']`: Harmony-corrected embeddings.
    *   `adata.uns`: Unstructured annotation, including `marker_names` (the 17-marker panel).
*   `seacells_pipeline_output_17_markers/processed_fcs_files.txt`: A log of successfully processed FCS files.
*   `seacells_pipeline_output_17_markers/skipped_fcs_files.txt`: A log of FCS files skipped due to missing markers or other issues.

### 1.4. Diagnostics/Visualizations

*   This script primarily focuses on data generation. Downstream scripts will generate visualizations based on this output.

## 2. VAE Model Training and Archetype Generation (Step 2)

This step trains the `VAEDisambiguation` model on the integrated metacells to discover cellular archetypes and their disambiguated marker profiles.

*   **Script:** `scripts/pipeline_step2_train_model.py`
*   **Purpose:** To learn a low-dimensional representation of metacells and a set of biologically interpretable archetypes, while disambiguating merged markers using biological priors.

### 2.1. Input Files

*   `seacells_pipeline_output_17_markers/integrated_metacells_17_markers.h5ad`: Output from Step 1.

### 2.2. Key Parameters (in `pipeline_step2_train_model.py`)

*   `--input`: Path to the input `.h5ad` file.
*   `--output-dir`: Directory to save all outputs from this step.
    *   **Default:** `vae_output` (but typically overridden to `d28_vae_analysis` or similar).
    *   **Impact:** Creates a new output directory.
*   `--n-archetypes`: Number of archetypes the VAE should learn.
    *   **Default:** `12`
    *   **Impact:** Crucial parameter. Affects the granularity of biological states identified. Changing this requires re-running downstream analyses.
*   `--epochs`: Number of training epochs for the VAE.
    *   **Default:** `100`
    *   **Impact:** More epochs can lead to better convergence but also overfitting.
*   `--latent-dim`: Dimensionality of the VAE's latent space.
    *   **Default:** `10`
    *   **Impact:** Affects the complexity of the learned representation.
*   `--beta`: Weight for the KL divergence term in the VAE loss.
    *   **Default:** `0.01`
    *   **Impact:** Controls the balance between reconstruction accuracy and latent space regularization.
*   `--gamma`: Weight for the biological prior loss.
    *   **Default:** `0.5`
    *   **Impact:** Controls how strongly biological priors (e.g., mutual exclusivity of markers) are enforced during disambiguation. **This is a critical parameter for improving biological interpretability.**
*   `--beta-warmup`: Number of epochs over which `beta` is linearly increased.
    *   **Default:** `100`
    *   **Impact:** Helps stabilize training.
*   `--gamma-cooldown`: Number of epochs over which `gamma` is linearly decreased.
    *   **Default:** `150`
    *   **Impact:** Allows the model to first learn general features, then refine with stronger biological constraints.
*   `--gamma-min`: Minimum value for `gamma` during cooldown.
    *   **Default:** `0.0`
    *   **Impact:** Ensures `gamma` doesn't drop to zero if a minimum influence is desired.
*   `--soft-assignment-temp`: Temperature for soft archetype assignments.
    *   **Default:** `1.0`
    *   **Impact:** Affects the "crispness" of archetype assignments. Higher values lead to softer assignments.

### 2.3. Outputs

*   `d28_vae_analysis/metacells_with_archetypes.h5ad`: An AnnData object containing:
    *   `adata.obsm['X_vae']`: VAE latent space embeddings.
    *   `adata.obs['archetype_X_prob']`: Soft assignment probabilities for each archetype (e.g., `archetype_0_prob`).
    *   `adata.obs['archetype_entropy']`: Entropy of archetype assignments (measure of "interstitialness").
    *   `adata.obs['archetype_id_hard']`: Hard assignment to the most probable archetype.
    *   `adata.obs['is_archetype']`: Boolean indicating if a metacell is an identified archetype center.
    *   `adata.obsm['X_umap']`: UMAP embeddings based on `X_vae`.
*   `d28_vae_analysis/archetype_profiles_disambiguated.csv`: A CSV file containing the disambiguated (un-merged) marker expression profiles for each learned archetype. This is key for biological interpretation.
*   `d28_vae_analysis/vae_model.pt`: The saved PyTorch state dictionary of the trained VAE model and the `StandardScaler` used for input data.

### 2.4. Diagnostics/Visualizations

*   **Loss Curves (Implicit):** The script prints training loss components (total, reconstruction, KL, biological prior loss) every 50 epochs, allowing monitoring of convergence.
*   **`archetype_profiles_disambiguated.csv`:** Direct diagnostic for biological interpretability of archetypes.

## 3. Archetype and Interstitial Analysis (Step 3 - Analysis Scripts)

This suite of scripts performs in-depth analysis of the VAE archetypes and the interstitial cell populations.

### 3.1. Input Files

*   `d28_vae_analysis/metacells_with_archetypes.h5ad`: Output from Step 2.
*   `d28_vae_analysis/archetype_usage_statistics.txt`: Statistical results from archetype usage comparison.
*   `d28_vae_analysis/interstitial_cells_statistics.txt`: Statistical results from interstitial cell comparison.

### 3.2. Scripts & Outputs

#### 3.2.1. `analyze_archetype_pairs.py`

*   **Purpose:** To analyze the relationships between pairs of archetypes, focusing on the 2D distribution of metacells that lie between them.
*   **Outputs:**
    *   `d28_vae_analysis/archetype_pair_distributions/archetype_pair_X_Y_2d_distribution.png` (for top 10 pairs): Scatter plots showing `proximity_ratio` (x-axis) vs. `archetype_entropy` (y-axis), colored by cohort. These plots are crucial for understanding the continuum of cell states between archetypes and how different cohorts populate these interstitial spaces.
*   **Key Parameters:** `--input`, `--output-dir`

#### 3.2.2. `compare_archetype_usage.py`

*   **Purpose:** Performs statistical comparisons of archetype usage across defined cohorts.
*   **Outputs:**
    *   `d28_vae_analysis/archetype_usage_statistics.txt`: A text file detailing Kruskal-Wallis test p-values and post-hoc Dunn's test p-values for each archetype across cohorts.
*   **Key Parameters:** `--input`, `--output-dir`

#### 3.2.3. `compare_interstitial_cells.py`

*   **Purpose:** Compares the proportion of interstitial cells (high entropy) across cohorts.
*   **Outputs:**
    *   `d28_vae_analysis/interstitial_cells_statistics.txt`: A text file detailing the proportion of interstitial cells per cohort and the Chi-squared test p-value.
*   **Key Parameters:** `--input`, `--output-dir`

#### 3.2.4. `plot_archetype_usage_by_cohort.py`

*   **Purpose:** Generates publication-ready bar plots showing the mean probability of each archetype per cohort, with error bars and statistical significance markers.
*   **Outputs:**
    *   `d28_vae_analysis/archetype_usage_plots/archetype_X_usage.png` (for each archetype): Bar plots.
*   **Key Parameters:** `--input`, `--stats-file`, `--output-dir`

#### 3.2.5. `plot_interstitial_proportions.py`

*   **Purpose:** Generates a publication-ready bar plot showing the proportion of interstitial cells per cohort, with the Chi-squared p-value.
*   **Outputs:**
    *   `d28_vae_analysis/interstitial_cells_proportions.png`: Bar plot.
*   **Key Parameters:** `--stats-file`, `--output-plot-file`

## 4. Clustering Analysis

This section details the scripts used for traditional unsupervised clustering analysis.

### 4.1. Input Files

*   `seacells_pipeline_output_17_markers/integrated_metacells_17_markers.h5ad`: Output from Step 1.

### 4.2. Scripts & Outputs

#### 4.2.1. `analyze_cluster_composition.py`

*   **Purpose:** Analyzes the cohort composition of each cluster (e.g., Phenograph or Leiden clusters).
*   **Outputs:**
    *   `d28_clustering_analysis/phenograph_composition.csv`: Cohort composition for Phenograph clusters.
    *   `d28_clustering_analysis/leiden_optimal_composition.csv`: Cohort composition for Leiden clusters.
*   **Key Parameters:** `--input`, `--output-dir`, `--cluster-method` (e.g., `phenograph`, `leiden`)

#### 4.2.2. `analyze_cluster_phenotypes.py`

*   **Purpose:** Characterizes the marker expression profiles of each cluster.
*   **Outputs:**
    *   `d28_clustering_analysis/phenograph_characterization_profiles.csv`: Mean marker expression for Phenograph clusters.
    *   `d28_clustering_analysis/leiden_characterization_profiles.csv`: Mean marker expression for Leiden clusters.
*   **Key Parameters:** `--input`, `--output-dir`, `--cluster-method`

## 5. Disambiguation Diagnostics

This script evaluates the quality of the VAE's marker disambiguation.

### 5.1. Input Files

*   `d28_vae_analysis/metacells_with_archetypes.h5ad`: Output from Step 2.
*   `d28_vae_analysis/archetype_profiles_disambiguated.csv`: Disambiguated archetype profiles from Step 2.

### 5.2. Script & Outputs

#### 5.2.1. `analyze_disambiguation_metrics.py`

*   **Purpose:** Assesses the biological plausibility of the disambiguated markers by checking for expected mutual exclusivity or co-expression patterns.
*   **Outputs:**
    *   `d28_vae_analysis/disambiguation_report.txt`: A text report summarizing correlations between disambiguated markers (e.g., CD3 vs CD19, CD4 vs CD8).
    *   `d28_vae_analysis/disambiguation_CDX_vs_CDY.png` (for selected pairs): Scatter plots of disambiguated marker expression.
*   **Key Parameters:** `--input-h5ad`, `--input-profiles`, `--output-dir`

## 6. Comprehensive Visualizations

These scripts generate various high-level visualizations to summarize the analysis.

### 6.1. Input Files

*   `d28_vae_analysis/metacells_with_archetypes.h5ad`: Output from Step 2.
*   `d28_vae_analysis/archetype_profiles_disambiguated.csv`: Disambiguated archetype profiles from Step 2.
*   `d28_clustering_analysis/metacells_with_clusters.h5ad`: AnnData object with clustering results (from `seacells_clustering_analysis.py`).

### 6.2. Scripts & Outputs

#### 6.2.1. `visualize_cluster_archetype_flow.py`

*   **Purpose:** Generates a Sankey diagram visualizing the flow of cells from clusters to archetypes, showing the relationship between discrete clusters and continuous archetypes.
*   **Outputs:**
    *   `d28_vae_analysis/cluster_archetype_sankey.html`: Interactive Sankey diagram.
*   **Key Parameters:** `--input-h5ad-vae`, `--input-h5ad-clusters`, `--output-dir`

#### 6.2.2. `visualize_disambiguated_markers.py`

*   **Purpose:** Projects the expression of disambiguated markers onto UMAP embeddings.
*   **Outputs:**
    *   `d28_vae_analysis/disambiguated_marker_umaps/CDX_umap.png` (for each un-merged marker): UMAP plots.
*   **Key Parameters:** `--input-h5ad`, `--input-profiles`, `--output-dir`

#### 6.2.3. `visualize_sample_archetypes.py`

*   **Purpose:** Generates a heatmap showing the archetype usage across individual samples.
*   **Outputs:**
    *   `d28_vae_analysis/sample_archetype_usage_heatmap.png`: Heatmap.
*   **Key Parameters:** `--input`, `--output-dir`

## 7. Running the Pipeline (General Instructions)

### 7.1. Environment Setup

Ensure you have the `flow_archetype_stable` conda environment set up. If not, create it using `environment.yml`.

```bash
conda env create -f environment.yml
conda activate flow_archetype_stable
```

### 7.2. Execution Order

The pipeline steps should generally be executed in the following order:

1.  **Step 1: Metacell Genesis**
    ```bash
    /opt/homebrew/Caskroom/mambaforge/base/bin/conda run -n flow_archetype_stable python scripts/pipeline_step1_seacells_resumable_v6.py \
      --output-dir seacells_pipeline_output_17_markers \
      --n-parallel 6 \
      --compression-ratio 40.0 \
      --log-resources
    ```
    *   **Note:** This script also performs initial Phenograph clustering and saves `seacells_pipeline_output_17_markers/metacells_with_clusters.h5ad` which is used by some visualization scripts.

2.  **Step 2: VAE Model Training**
    ```bash
    /opt/homebrew/Caskroom/mambaforge/base/bin/conda run -n flow_archetype_stable python scripts/pipeline_step2_train_model.py \
      --input seacells_pipeline_output_17_markers/integrated_metacells_17_markers.h5ad \
      --output-dir d28_vae_analysis \
      --epochs 200 \
      --gamma 1.0 # Consider increasing gamma for stronger prior enforcement
    ```
    *   **Note:** Adjust `--epochs` and `--gamma` as discussed in Section 2.2.

3.  **Analysis & Visualization Scripts (can be run in any order after Step 2)**

    *   **Archetype Usage Statistics:**
        ```bash
        /opt/homebrew/Caskroom/mambaforge/base/bin/conda run -n flow_archetype_stable python compare_archetype_usage.py \
          --input d28_vae_analysis/metacells_with_archetypes.h5ad \
          --output-dir d28_vae_analysis
        ```
    *   **Interstitial Cell Statistics:**
        ```bash
        /opt/homebrew/Caskroom/mambaforge/base/bin/conda run -n flow_archetype_stable python compare_interstitial_cells.py \
          --input d28_vae_analysis/metacells_with_archetypes.h5ad \
          --output-dir d28_vae_analysis
        ```
    *   **Plot Archetype Usage:**
        ```bash
        /opt/homebrew/Caskroom/mambaforge/base/bin/conda run -n flow_archetype_stable python plot_archetype_usage_by_cohort.py \
          --input d28_vae_analysis/metacells_with_archetypes.h5ad \
          --stats-file d28_vae_analysis/archetype_usage_statistics.txt \
          --output-dir d28_vae_analysis
        ```
    *   **Plot Interstitial Proportions:**
        ```bash
        /opt/homebrew/Caskroom/mambaforge/base/bin/conda run -n flow_archetype_stable python plot_interstitial_proportions.py \
          --stats-file d28_vae_analysis/interstitial_cells_statistics.txt \
          --output-plot-file d28_vae_analysis/interstitial_cells_proportions.png
        ```
    *   **Analyze Archetype Pairs (2D Distributions):**
        ```bash
        /opt/homebrew/Caskroom/mambaforge/base/bin/conda run -n flow_archetype_stable python analyze_archetype_pairs.py \
          --input d28_vae_analysis/metacells_with_archetypes.h5ad \
          --output-dir d28_vae_analysis
        ```
    *   **Disambiguation Metrics:**
        ```bash
        /opt/homebrew/Caskroom/mambaforge/base/bin/conda run -n flow_archetype_stable python analyze_disambiguation_metrics.py \
          --input-h5ad d28_vae_analysis/metacells_with_archetypes.h5ad \
          --input-profiles d28_vae_analysis/archetype_profiles_disambiguated.csv \
          --output-dir d28_vae_analysis
        ```
    *   **Visualize Cluster-Archetype Flow (Sankey):**
        ```bash
        /opt/homebrew/Caskroom/mambaforge/base/bin/conda run -n flow_archetype_stable python visualize_cluster_archetype_flow.py \
          --input-h5ad-vae d28_vae_analysis/metacells_with_archetypes.h5ad \
          --input-h5ad-clusters seacells_pipeline_output_17_markers/metacells_with_clusters.h5ad \
          --output-dir d28_vae_analysis
        ```
    *   **Visualize Disambiguated Markers on UMAP:**
        ```bash
        /opt/homebrew/Caskroom/mambaforge/base/bin/conda run -n flow_archetype_stable python visualize_disambiguated_markers.py \
          --input-h5ad d28_vae_analysis/metacells_with_archetypes.h5ad \
          --input-profiles d28_vae_analysis/archetype_profiles_disambiguated.csv \
          --output-dir d28_vae_analysis
        ```
    *   **Visualize Sample Archetype Usage:**
        ```bash
        /opt/homebrew/Caskroom/mambaforge/base/bin/conda run -n flow_archetype_stable python visualize_sample_archetypes.py \
          --input d28_vae_analysis/metacells_with_archetypes.h5ad \
          --output-dir d28_vae_analysis
        ```
    *   **Clustering Analysis (Phenograph/Leiden):**
        ```bash
        # For Phenograph
        /opt/homebrew/Caskroom/mambaforge/base/bin/conda run -n flow_archetype_stable python analyze_cluster_composition.py \
          --input seacells_pipeline_output_17_markers/metacells_with_clusters.h5ad \
          --output-dir d28_clustering_analysis \
          --cluster-method phenograph
        /opt/homebrew/Caskroom/mambaforge/base/bin/conda run -n flow_archetype_stable python analyze_cluster_phenotypes.py \
          --input seacells_pipeline_output_17_markers/metacells_with_clusters.h5ad \
          --output-dir d28_clustering_analysis \
          --cluster-method phenograph

        # For Leiden (if applicable, assuming Leiden clusters are also in metacells_with_clusters.h5ad)
        # Note: If Leiden clustering is not performed in pipeline_step1, you might need a separate script for it.
        # For this pipeline, Phenograph is generated in Step 1.
        ```

### 7.3. Updating the Pipeline

*   **Modifying Parameters:** Parameters for each script are exposed via `argparse`. Modify the command-line arguments when running the scripts.
*   **Changing Biological Priors:** The `DisambiguationBiologicalPriorLoss` class in `scripts/pipeline_step2_train_model.py` defines the biological priors. To add new priors or adjust their weights, modify this class directly.
    *   **Example:** To add mutual exclusivity between CD14 and CD3:
        ```python
        if 'CD14' in self.marker_idx and 'CD3' in self.marker_idx:
            cd14 = x_disambiguated[:, self.marker_idx['CD14']]
            cd3 = x_disambiguated[:, self.marker_idx['CD3']]
            total_loss += self.weights.get('cd14_cd3_exclusion', 1.0) * (torch.relu(cd14) * torch.relu(cd3)).mean()
        ```
        And add `'cd14_cd3_exclusion': 1.0` to the `weights` dictionary in the `__init__` method.
*   **Model Architecture:** The `VAEDisambiguation` class in `scripts/pipeline_step2_train_model.py` defines the VAE architecture. Changes to `latent_dim`, `hidden_dims`, or other layers require modifying this class.
*   **Adding New Analysis/Visualizations:** Create new Python scripts following the existing patterns, taking `.h5ad` files and other outputs as inputs, and saving new plots/reports to appropriate output directories.

### 7.4. File Locations Summary

*   **Raw Data:** `data/abnormalfcs_gated/`, `data/normalsgated/`
*   **Metadata:** `sample_annotations.csv`
*   **Scripts:** `scripts/` directory (e.g., `pipeline_step1_seacells_resumable_v6.py`, `pipeline_step2_train_model.py`, `analyze_archetype_pairs.py`, etc.)
*   **Step 1 Outputs:** `seacells_pipeline_output_17_markers/`
    *   `integrated_metacells_17_markers.h5ad`
    *   `metacells_with_clusters.h5ad` (contains Phenograph clusters from Step 1)
*   **Step 2 Outputs (VAE):** `d28_vae_analysis/` (or custom `--output-dir`)
    *   `metacells_with_archetypes.h5ad`
    *   `archetype_profiles_disambiguated.csv`
    *   `vae_model.pt`
    *   `archetype_usage_statistics.txt`
    *   `interstitial_cells_statistics.txt`
    *   `disambiguation_report.txt`
    *   `cluster_archetype_sankey.html`
    *   `sample_archetype_usage_heatmap.png`
    *   `interstitial_cells_proportions.png`
    *   `archetype_usage_plots/` (directory for individual archetype usage plots)
    *   `archetype_pair_distributions/` (directory for 2D distribution plots)
    *   `disambiguated_marker_umaps/` (directory for UMAPs of disambiguated markers)
*   **Clustering Outputs:** `d28_clustering_analysis/` (or custom `--output-dir`)
    *   `phenograph_composition.csv`
    *   `phenograph_characterization_profiles.csv`
    *   (Similar files for Leiden if generated)

This guide should provide a comprehensive understanding and operational instructions for the metacell analysis pipeline.
