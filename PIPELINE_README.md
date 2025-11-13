# SEACells + Harmony + Clustering + VAE Pipeline

Complete workflow for flow cytometry metacell analysis with batch correction and deep learning.

## Overview

This pipeline implements the following workflow:

1. **SEACells Metacell Generation**: Process individual samples with SEACells (n/75 cells per metacell)
2. **Harmony Batch Correction**: Integrate metacells across samples using Harmony
3. **Clustering Analysis**: Apply PhenoGraph and Leiden clustering to identify cell populations
4. **VAE Training**: Train a Variational Autoencoder on the integrated metacells with Optuna optimization

## Pipeline Architecture

```
FCS Files (Normal + Abnormal)
    ↓
SEACells (per-sample, n/75)
    ↓
Harmony Batch Correction
    ↓
Integrated Metacells
    ↓
├─→ Clustering (PhenoGraph + Leiden)
│       ↓
│   Cluster Assignments & Visualizations
│
└─→ VAE Training (with Optuna)
        ↓
    Learned Archetypes & Reconstructions
```

## Requirements

### Data Structure
```
data/
├── normalsgated/          # Normal sample FCS files
└── deidentify_abngated/   # Abnormal sample FCS files
```

### Dependencies

All dependencies are managed via conda environment (see `environment.yml`):
- Python 3.10
- SEACells 0.3.3
- Harmony 0.0.10
- PhenoGraph
- scanpy, anndata
- PyTorch
- Optuna
- And more (see `environment.yml` for complete list)

## Usage

### Option 1: Local Execution

#### Setup
```bash
# Create conda environment
conda env create -f environment.yml
conda activate metacell_vae_env

# Install SEACells from local directory
cd SEACells
pip install -e .
cd ..
```

#### Run Full Pipeline
```bash
# Using Python
python run_pipeline.py --n-trials 50

# Or using shell script
./run_pipeline.sh --n-trials 50
```

#### Run Specific Steps
```bash
# Skip SEACells if already generated
python run_pipeline.py --skip-seacells --n-trials 50

# Skip clustering analysis
python run_pipeline.py --skip-clustering --n-trials 50

# Run only clustering on existing data
python run_pipeline.py --skip-seacells --skip-clustering --n-trials 0
```

### Option 2: Docker Execution

#### Build Docker Image
```bash
docker build -t seacells-vae-pipeline .
```

#### Run Pipeline in Docker
```bash
# Full pipeline
docker run -v $(pwd)/data:/app/data \
           -v $(pwd)/seacells_output:/app/seacells_output \
           -v $(pwd)/metacell_vae_output:/app/metacell_vae_output \
           seacells-vae-pipeline --n-trials 50

# Skip SEACells step (use existing data)
docker run -v $(pwd)/seacells_output:/app/seacells_output \
           -v $(pwd)/metacell_vae_output:/app/metacell_vae_output \
           seacells-vae-pipeline --skip-seacells --n-trials 50
```

## Pipeline Details

### Step 1: SEACells + Harmony Batch Correction
**Script**: `scripts/seacells_batch_correction.py`

- Processes each FCS file individually with SEACells
- Uses default metacell size of n/75 (75 cells per metacell)
- Applies arcsinh transformation (scale factor: 5.0)
- Z-score normalization per marker
- PCA dimensionality reduction
- Harmony integration across samples
- Generates UMAP visualization

**Outputs**:
- `seacells_output/batch_corrected/[sample]_metacells.h5ad` - Individual metacell files
- `seacells_output/batch_corrected/integrated_metacells_harmony.h5ad` - Integrated data
- `seacells_output/batch_corrected/harmony_integration_overview.png` - Visualization

### Step 2: Clustering Analysis
**Script**: `scripts/seacells_clustering_analysis.py`

- **PhenoGraph**: Community detection with k=30 nearest neighbors
- **Leiden**: Multiple resolutions tested (0.5, 1.0, 1.5, 2.0)
- Clustering quality metrics: Silhouette, Calinski-Harabasz, Davies-Bouldin
- Cluster characterization by marker expression
- Normal/Abnormal composition analysis

**Outputs**:
- `seacells_output/clustering_analysis/umap_phenograph.png`
- `seacells_output/clustering_analysis/umap_leiden_optimal.png`
- `seacells_output/clustering_analysis/leiden_metrics.csv`
- `seacells_output/clustering_analysis/clustering_summary.txt`
- `seacells_output/clustering_analysis/leiden_characterization.png`
- `seacells_output/clustering_analysis/metacells_with_clusters.h5ad`

### Step 3: VAE Training
**Script**: `scripts/metacell_vae.py`

- Optuna-based hyperparameter optimization
- Heat diffusion preprocessing option
- Multiple architecture choices (standard, hierarchical, etc.)
- Reconstruction loss tracking
- Latent space visualization

**Outputs**:
- `metacell_vae_output/best_model.pt` - Trained model
- `metacell_vae_output/optuna_study.db` - Optimization history
- `metacell_vae_output/training_curves.png`
- `metacell_vae_output/latent_space_umap.png`

## Configuration

### SEACells Parameters
Modify in `scripts/seacells_batch_correction.py`:
- `n_cells_per_sample`: Cells to subsample per sample (default: 7500)
- `target_total_metacells`: Total metacells to generate (default: 1000)
- `n_parallel`: Parallel processes (default: 4)

### Clustering Parameters
```bash
python scripts/seacells_clustering_analysis.py \
    --input seacells_output/batch_corrected/integrated_metacells_harmony.h5ad \
    --output-dir seacells_output/clustering_analysis \
    --phenograph-k 30 \
    --leiden-resolutions 0.5 1.0 1.5 2.0
```

### VAE Parameters
```bash
python scripts/metacell_vae.py \
    --n-trials 50 \
    --input-file seacells_output/batch_corrected/integrated_metacells_harmony.h5ad
```

## Troubleshooting

### Issue: Out of Memory
- Reduce `n_cells_per_sample` in `seacells_batch_correction.py`
- Reduce `target_total_metacells`
- Use fewer parallel processes (`n_parallel`)

### Issue: PhenoGraph Not Found
```bash
conda activate metacell_vae_env
pip install phenograph
```

### Issue: SEACells Import Error
```bash
cd SEACells
pip install -e .
```

### Issue: Docker Build Fails
- Ensure sufficient disk space (>10GB)
- Check Docker memory allocation (>8GB recommended)
- Try building with `--no-cache` flag

## Output Summary

After successful completion, you'll have:

```
seacells_output/
├── batch_corrected/
│   ├── integrated_metacells_harmony.h5ad    # Main integrated data
│   ├── combined_metacells_pre_harmony.h5ad
│   └── [sample]_metacells.h5ad              # Individual samples
└── clustering_analysis/
    ├── metacells_with_clusters.h5ad         # Data with cluster labels
    ├── clustering_summary.txt               # Summary statistics
    ├── leiden_metrics.csv                   # Quality metrics
    └── [visualizations].png                 # Multiple plots

metacell_vae_output/
├── best_model.pt                            # Trained VAE
├── optuna_study.db                          # Optimization database
└── [results and visualizations]
```

## References

- SEACells: Persad et al., Nature Biotechnology 2023
- Harmony: Korsunsky et al., Nature Methods 2019
- PhenoGraph: Levine et al., Cell 2015
- Leiden: Traag et al., Scientific Reports 2019

## Support

For issues or questions:
1. Check the pipeline logs for detailed error messages
2. Verify input data structure and formats
3. Ensure all dependencies are installed correctly
4. Review individual script documentation in `scripts/`

---
*Generated by Claude Code*
*Last Updated: 2024-11-12*
