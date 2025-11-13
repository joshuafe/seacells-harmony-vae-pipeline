#!/bin/bash
# SEACells + Harmony + Clustering + VAE Pipeline
# This script provides a shell wrapper for the Python pipeline

set -e  # Exit on any error

# Default values
SKIP_SEACELLS=false
SKIP_CLUSTERING=false
N_TRIALS=50
COMPRESSION_RATIO=50.0

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --skip-seacells)
      SKIP_SEACELLS=true
      shift
      ;;
    --skip-clustering)
      SKIP_CLUSTERING=true
      shift
      ;;
    --n-trials)
      N_TRIALS="$2"
      shift 2
      ;;
    --compression-ratio)
      COMPRESSION_RATIO="$2"
      shift 2
      ;;
    --help)
      echo "Usage: $0 [OPTIONS]"
      echo ""
      echo "Options:"
      echo "  --skip-seacells       Skip SEACells and Harmony step"
      echo "  --skip-clustering     Skip clustering analysis step"
      echo "  --n-trials N          Number of Optuna trials (default: 50)"
      echo "  --compression-ratio R Cells per metacell (default: 50.0)"
      echo "  --help                Show this help message"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      echo "Use --help for usage information"
      exit 1
      ;;
  esac
done

# Activate conda environment if not already active
if [[ -z "$CONDA_DEFAULT_ENV" ]] || [[ "$CONDA_DEFAULT_ENV" != "metacell_vae_env" ]]; then
    echo "Activating conda environment: metacell_vae_env"
    source activate metacell_vae_env
fi

echo "================================================================================"
echo "SEACELLS + HARMONY + CLUSTERING + VAE PIPELINE"
echo "================================================================================"
echo ""
echo "Workflow:"
echo "  1. SEACells metacell generation (per-sample) + Harmony batch correction"
echo "  2. PhenoGraph and Leiden clustering analysis"
echo "  3. Variational Autoencoder (VAE) training with Optuna"
echo "================================================================================"

# Build command arguments
CMD_ARGS=()
if [ "$SKIP_SEACELLS" = true ]; then
    CMD_ARGS+=(--skip-seacells)
fi
if [ "$SKIP_CLUSTERING" = true ]; then
    CMD_ARGS+=(--skip-clustering)
fi
CMD_ARGS+=(--n-trials "$N_TRIALS")
CMD_ARGS+=(--compression-ratio "$COMPRESSION_RATIO")

# Run the Python pipeline
python run_pipeline.py "${CMD_ARGS[@]}"

echo ""
echo "================================================================================"
echo "Pipeline execution completed!"
echo "================================================================================"
