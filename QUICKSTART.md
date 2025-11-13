# Quick Start Guide

This guide will get you running the SEACells + Harmony + Clustering + VAE pipeline in minutes.

## Prerequisites

- **Data**: FCS files in `data/normalsgated/` and `data/deidentify_abngated/`
- **System**: 16GB+ RAM recommended, 8+ CPU cores optimal
- **Choice**: Docker (easier) or Conda (more control)

## Quick Start: Docker (Recommended)

### 1. Build the Docker image (one-time setup)
```bash
docker build -t seacells-vae-pipeline .
```
*This takes 10-20 minutes on first build*

### 2. Run the full pipeline
```bash
docker run \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/seacells_output:/app/seacells_output \
  -v $(pwd)/metacell_vae_output:/app/metacell_vae_output \
  seacells-vae-pipeline --n-trials 50
```

### 3. Check outputs
```bash
# Integrated metacells
ls -lh seacells_output/batch_corrected/integrated_metacells_harmony.h5ad

# Clustering results
ls -lh seacells_output/clustering_analysis/

# VAE results
ls -lh metacell_vae_output/
```

## Quick Start: Local Conda

### 1. Create environment (one-time setup)
```bash
conda env create -f environment.yml
conda activate metacell_vae_env
cd SEACells && pip install -e . && cd ..
```
*This takes 5-10 minutes*

### 2. Run the full pipeline
```bash
python run_pipeline.py --n-trials 50
```

**OR** use the shell script:
```bash
./run_pipeline.sh --n-trials 50
```

### 3. Check outputs
Same as Docker instructions above.

## Common Usage Patterns

### Pattern 1: Full Pipeline (First Time)
```bash
# Docker
docker run -v $(pwd)/data:/app/data \
           -v $(pwd)/seacells_output:/app/seacells_output \
           -v $(pwd)/metacell_vae_output:/app/metacell_vae_output \
           seacells-vae-pipeline --n-trials 50

# Local
python run_pipeline.py --n-trials 50
```

### Pattern 2: Re-run Only VAE (SEACells Already Done)
```bash
# Docker
docker run -v $(pwd)/seacells_output:/app/seacells_output \
           -v $(pwd)/metacell_vae_output:/app/metacell_vae_output \
           seacells-vae-pipeline --skip-seacells --skip-clustering --n-trials 100

# Local
python run_pipeline.py --skip-seacells --skip-clustering --n-trials 100
```

### Pattern 3: Re-run Clustering + VAE
```bash
# Docker
docker run -v $(pwd)/seacells_output:/app/seacells_output \
           -v $(pwd)/metacell_vae_output:/app/metacell_vae_output \
           seacells-vae-pipeline --skip-seacells --n-trials 50

# Local
python run_pipeline.py --skip-seacells --n-trials 50
```

### Pattern 4: Only Clustering (No VAE)
```bash
# Local only (Docker not needed for quick analysis)
python scripts/seacells_clustering_analysis.py \
  --input seacells_output/batch_corrected/integrated_metacells_harmony.h5ad \
  --output-dir seacells_output/clustering_analysis
```

## Expected Runtime

| Step | Time | Notes |
|------|------|-------|
| SEACells + Harmony | 15-30 min | Depends on # samples, uses 4 parallel processes |
| Clustering | 5-10 min | PhenoGraph + Leiden on ~1000 metacells |
| VAE Training (50 trials) | 1-2 hours | With Optuna optimization |

**Total for full pipeline**: ~2-3 hours

## Monitoring Progress

### Docker Logs
```bash
# Follow logs in real-time
docker logs -f <container_id>
```

### Local Progress
The pipeline prints detailed progress messages:
- SEACells: Per-sample progress with timestamps
- Clustering: Quality metrics for each resolution
- VAE: Trial-by-trial Optuna updates

## Outputs at Each Stage

### After Step 1 (SEACells + Harmony)
```
seacells_output/batch_corrected/
├── integrated_metacells_harmony.h5ad          ← Main output
├── harmony_integration_overview.png           ← Visualization
└── [Individual sample files...]
```

### After Step 2 (Clustering)
```
seacells_output/clustering_analysis/
├── metacells_with_clusters.h5ad               ← Clustered data
├── clustering_summary.txt                     ← Summary report
├── leiden_metrics.csv                         ← Quality metrics
├── umap_leiden_optimal.png                    ← Best clustering
├── leiden_characterization.png                ← Marker heatmap
└── leiden_composition.png                     ← Normal/Abnormal composition
```

### After Step 3 (VAE)
```
metacell_vae_output/
├── best_model.pt                              ← Trained model
├── optuna_study.db                            ← Optimization history
├── training_curves.png                        ← Loss plots
└── [Additional outputs...]
```

## Troubleshooting

### "Out of memory" error
```bash
# Reduce cells per sample
# Edit scripts/seacells_batch_correction.py:
# Change: n_cells_per_sample = 7500  →  n_cells_per_sample = 5000
# Change: target_total_metacells = 1000  →  target_total_metacells = 750
```

### "PhenoGraph not found"
```bash
conda activate metacell_vae_env
pip install phenograph
```

### Docker permission errors
```bash
# Linux only - fix output directory permissions
sudo chown -R $USER:$USER seacells_output metacell_vae_output
```

## Next Steps

After the pipeline completes:

1. **Review Clustering Results**
   - Open `seacells_output/clustering_analysis/clustering_summary.txt`
   - Check UMAP visualizations
   - Examine marker characterization heatmap

2. **Analyze VAE Results**
   - Review Optuna optimization history
   - Inspect learned latent representations
   - Evaluate reconstruction quality

3. **Further Analysis**
   - Load `metacells_with_clusters.h5ad` in Python/scanpy
   - Perform differential expression analysis
   - Identify marker signatures for each cluster

## Help

```bash
# Python script help
python run_pipeline.py --help

# Shell script help
./run_pipeline.sh --help

# Individual scripts
python scripts/seacells_batch_correction.py --help
python scripts/seacells_clustering_analysis.py --help
python scripts/metacell_vae.py --help
```

For detailed documentation, see: `PIPELINE_README.md`

---
*Ready to analyze flow cytometry data at scale!*
