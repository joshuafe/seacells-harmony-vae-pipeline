# SEACells + Harmony + Clustering + VAE Pipeline

Complete workflow for flow cytometry metacell analysis with batch correction and deep learning.

## Overview

This pipeline processes flow cytometry data through metacell generation, batch correction, clustering, and variational autoencoder training to identify cellular archetypes and patterns.

### Pipeline Architecture

```
FCS Files (Normal + Abnormal)
    ↓
SEACells Metacell Generation (per-sample, n/75)
    ↓
Harmony Batch Correction
    ↓
Integrated Metacells (~2500 across all samples)
    ↓
├─→ Clustering (PhenoGraph + Leiden)
│       ↓
│   Cell Population Identification
│
└─→ VAE Training (Optuna Optimization)
        ↓
    Learned Archetypes & Latent Space
```

## Features

- **Metacell Generation**: SEACells algorithm with per-sample processing
- **Batch Correction**: Harmony integration across multiple samples
- **Clustering**: PhenoGraph + Leiden with quality metrics
- **Deep Learning**: VAE with Optuna hyperparameter optimization
- **Docker Support**: Containerized for reproducibility
- **Comprehensive Logging**: Full pipeline monitoring

## Quick Start

### Prerequisites

```bash
# Data structure (not included in repo)
data/
├── normalsgated/          # Normal sample FCS files
└── deidentify_abngated/   # Abnormal sample FCS files
```

### Option 1: Docker (Recommended)

```bash
# Build image
./docker-run.sh build

# Run full pipeline
./docker-run.sh run --n-trials 50 2>&1 | tee logs/pipeline.log
```

### Option 2: Local Conda

```bash
# Create environment
conda env create -f environment.yml
conda activate metacell_vae_env

# Install SEACells
cd SEACells && pip install -e . && cd ..

# Run pipeline
python run_pipeline.py --n-trials 50 2>&1 | tee logs/pipeline.log
```

## Documentation

- **[QUICKSTART.md](QUICKSTART.md)** - Get started in 5 minutes
- **[PIPELINE_README.md](PIPELINE_README.md)** - Complete technical documentation

## Pipeline Steps

### Step 1: SEACells + Harmony
```bash
python scripts/seacells_batch_correction.py
```
- Processes each sample independently with SEACells
- Generates ~83 metacells per sample (2500 total)
- Applies Harmony batch correction
- Outputs: `integrated_metacells_harmony.h5ad`

### Step 2: Clustering Analysis
```bash
python scripts/seacells_clustering_analysis.py \
  --input seacells_output/batch_corrected/integrated_metacells_harmony.h5ad
```
- PhenoGraph community detection
- Leiden clustering at multiple resolutions
- Quality metrics (Silhouette, Calinski-Harabasz, Davies-Bouldin)
- Cluster characterization and visualization

### Step 3: VAE Training
```bash
python scripts/metacell_vae.py --n-trials 50
```
- Optuna-based hyperparameter optimization
- Heat diffusion preprocessing
- Latent space learning
- Archetype extraction

## Project Structure

```
.
├── scripts/                          # Pipeline scripts
│   ├── seacells_batch_correction.py # Step 1: SEACells + Harmony
│   ├── seacells_clustering_analysis.py # Step 2: Clustering
│   └── metacell_vae.py              # Step 3: VAE training
├── SEACells/                         # Local SEACells package
├── run_pipeline.py                   # Main pipeline orchestrator
├── run_pipeline.sh                   # Shell wrapper
├── docker-run.sh                     # Docker helper script
├── Dockerfile                        # Container definition
├── environment.yml                   # Conda environment
├── QUICKSTART.md                     # Quick start guide
├── PIPELINE_README.md                # Detailed documentation
└── README.md                         # This file
```

## Configuration

### Pipeline Parameters

Edit `scripts/seacells_batch_correction.py`:
```python
n_cells_per_sample = 7500      # Cells to subsample per sample
target_total_metacells = 2500  # Total metacells across all samples
n_parallel = 4                 # Parallel processes
```

### Command-Line Options

```bash
python run_pipeline.py --help

Options:
  --skip-seacells       Skip SEACells and Harmony step
  --skip-clustering     Skip clustering analysis
  --n-trials N          Number of Optuna trials (default: 50)
```

## Output Structure

```
seacells_output/
├── batch_corrected/
│   ├── integrated_metacells_harmony.h5ad    # Main integrated data
│   ├── combined_metacells_pre_harmony.h5ad
│   ├── harmony_integration_overview.png
│   └── [Individual sample metacells...]
└── clustering_analysis/
    ├── metacells_with_clusters.h5ad
    ├── clustering_summary.txt
    ├── leiden_metrics.csv
    └── [Visualization plots...]

metacell_vae_output/
├── best_model.pt
├── optuna_study.db
└── [Training results...]
```

## Dependencies

Key packages (see `environment.yml` for complete list):
- Python 3.10
- SEACells 0.3.3
- Harmony 0.0.10
- PhenoGraph
- scanpy, anndata
- PyTorch
- Optuna

## Expected Runtime

| Step | Time | Notes |
|------|------|-------|
| SEACells + Harmony | 30-45 min | 30 samples, 4 parallel processes |
| Clustering | 5-10 min | ~2500 metacells |
| VAE Training | 2-3 hours | 50 Optuna trials |
| **Total** | **3-4 hours** | Full pipeline |

## Docker Commands

```bash
# Build
docker build -t seacells-vae-pipeline .

# Run full pipeline
docker run --rm \
  -v $(pwd)/data:/app/data:ro \
  -v $(pwd)/seacells_output:/app/seacells_output \
  -v $(pwd)/metacell_vae_output:/app/metacell_vae_output \
  seacells-vae-pipeline --n-trials 50

# Interactive shell
docker run -it --rm \
  -v $(pwd)/data:/app/data:ro \
  -v $(pwd)/seacells_output:/app/seacells_output \
  -v $(pwd)/metacell_vae_output:/app/metacell_vae_output \
  --entrypoint /bin/bash \
  seacells-vae-pipeline
```

## Version Control

This project uses Git for version control. Data files and outputs are excluded via `.gitignore`.

### Workflow
```bash
# Make changes
git add .
git commit -m "Description of changes"
git push origin main

# Create feature branch
git checkout -b feature/new-analysis
# ... make changes ...
git commit -m "Add new analysis"
git push origin feature/new-analysis
```

## Troubleshooting

### Out of Memory
- Reduce `n_cells_per_sample` to 5000
- Reduce `target_total_metacells` to 1500
- Use fewer parallel processes

### Missing Dependencies
```bash
conda activate metacell_vae_env
pip install phenograph pynndescent
```

### Docker Issues
- Ensure >8GB RAM allocated to Docker
- Check disk space (>10GB needed)
- Try `docker system prune` if build fails

## References

- **SEACells**: Persad et al., Nature Biotechnology 2023
- **Harmony**: Korsunsky et al., Nature Methods 2019
- **PhenoGraph**: Levine et al., Cell 2015
- **Leiden**: Traag et al., Scientific Reports 2019

## License

[Add your license here]

## Citation

If you use this pipeline, please cite:
```
[Add citation information]
```

## Contact

[Add contact information]

---

**Last Updated**: 2024-11-12
**Maintained by**: [Your name/team]
