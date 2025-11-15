#!/bin/bash

# Test script for gamma cooldown feature
# Tests the new VAE framework with biological priors

echo "========================================================================"
echo "Testing VAE Archetypes with Gamma Cooldown"
echo "========================================================================"

# Set conda environment
CONDA_RUN="/opt/homebrew/Caskroom/mambaforge/base/bin/conda run -n flow_archetype_stable"

# Input data
INPUT="seacells_output/three_groups_phenograph/integrated_metacells_harmony_phenograph.h5ad"

# Test 1: Quick test with gamma cooldown (DEFAULT)
echo ""
echo "TEST 1: Quick run with gamma cooldown"
echo "----------------------------------------------------------------------"
$CONDA_RUN python scripts/vae_archetypes_with_priors.py \
  --input $INPUT \
  --train-on-reference \
  --latent-dim 5 \
  --n-archetypes 10 \
  --beta 0.01 \
  --gamma 0.5 \
  --gamma-cooldown 100 \
  --gamma-min 0.0 \
  --epochs 150 \
  --output-dir test_output/gamma_cooldown

echo ""
echo "✓ Test 1 complete"
echo ""

# Test 2: Constant gamma (no cooldown) for comparison
echo "TEST 2: Quick run with constant gamma (no cooldown)"
echo "----------------------------------------------------------------------"
$CONDA_RUN python scripts/vae_archetypes_with_priors.py \
  --input $INPUT \
  --train-on-reference \
  --latent-dim 5 \
  --n-archetypes 10 \
  --beta 0.01 \
  --gamma 0.1 \
  --gamma-cooldown 0 \
  --epochs 150 \
  --output-dir test_output/constant_gamma

echo ""
echo "✓ Test 2 complete"
echo ""

echo "========================================================================"
echo "All tests complete!"
echo "========================================================================"
echo ""
echo "Outputs:"
echo "  - test_output/gamma_cooldown/       (with cooldown)"
echo "  - test_output/constant_gamma/       (without cooldown)"
echo ""
echo "Compare:"
echo "  - Training curves: check gamma schedule plot"
echo "  - Archetypes: which produces more interpretable results?"
echo "  - Entropy: which has lower average entropy (clearer assignments)?"
echo ""
