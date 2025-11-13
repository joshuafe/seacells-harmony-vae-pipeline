"""
Complete SEACells + Harmony + Clustering + VAE Pipeline

This pipeline follows the workflow:
1. SEACells metacell generation (per-sample) + Harmony batch correction
2. PhenoGraph and Leiden clustering analysis
3. Variational Autoencoder (VAE) training

Author: Claude Code
Date: 2024
"""

import subprocess
import sys
import argparse
from pathlib import Path

def run_command(command):
    """Runs a command and prints its output in real-time."""
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip())
    rc = process.poll()
    return rc

def main():
    parser = argparse.ArgumentParser(
        description="Complete SEACells + Harmony + Clustering + VAE Pipeline"
    )
    parser.add_argument(
        "--skip-seacells",
        action="store_true",
        help="Skip SEACells and Harmony step (use existing integrated_metacells_harmony.h5ad)"
    )
    parser.add_argument(
        "--skip-clustering",
        action="store_true",
        help="Skip clustering analysis step"
    )
    parser.add_argument(
        "--n-trials",
        type=int,
        default=50,
        help="Number of Optuna trials for VAE hyperparameter optimization"
    )
    parser.add_argument(
        "--compression-ratio",
        type=float,
        default=50.0,
        help="SEACells compression ratio (cells per metacell). Default: 50.0"
    )
    args = parser.parse_args()

    print("=" * 80)
    print("SEACELLS + HARMONY + CLUSTERING + VAE PIPELINE")
    print("=" * 80)
    print("\nWorkflow:")
    print("  1. SEACells metacell generation (per-sample) + Harmony batch correction")
    print("  2. PhenoGraph and Leiden clustering analysis")
    print("  3. Variational Autoencoder (VAE) training with Optuna")
    print("=" * 80)

    # Step 1: Run SEACells per sample + Harmony batch correction
    if not args.skip_seacells:
        print("\n" + "=" * 80)
        print("STEP 1: SEACells + Harmony Batch Correction")
        print(f"Compression Ratio: {args.compression_ratio}:1")
        print("=" * 80)
        cmd1 = ["python", "scripts/seacells_batch_correction.py",
                "--compression-ratio", str(args.compression_ratio)]
        if run_command(cmd1) != 0:
            print("\nERROR: SEACells + Harmony step failed.")
            sys.exit(1)
        print("\n✓ Step 1 completed successfully")
    else:
        print("\n⊳ Skipping Step 1 (using existing integrated metacells)")

    # Step 2: Run clustering analysis
    if not args.skip_clustering:
        print("\n" + "=" * 80)
        print("STEP 2: Clustering Analysis (PhenoGraph + Leiden)")
        print("=" * 80)
        cmd2 = ["python", "scripts/seacells_clustering_analysis.py",
                "--input", "seacells_output/batch_corrected/integrated_metacells_harmony.h5ad",
                "--output-dir", "seacells_output/clustering_analysis"]
        if run_command(cmd2) != 0:
            print("\nERROR: Clustering analysis failed.")
            sys.exit(1)
        print("\n✓ Step 2 completed successfully")
    else:
        print("\n⊳ Skipping Step 2 (clustering analysis)")

    # Step 3: Run Metacell VAE
    print("\n" + "=" * 80)
    print("STEP 3: Variational Autoencoder (VAE) Training")
    print("=" * 80)
    cmd3 = ["python", "scripts/metacell_vae.py",
            "--n-trials", str(args.n_trials)]
    if run_command(cmd3) != 0:
        print("\nERROR: Metacell VAE training failed.")
        sys.exit(1)
    print("\n✓ Step 3 completed successfully")

    print("\n" + "=" * 80)
    print("PIPELINE COMPLETED SUCCESSFULLY!")
    print("=" * 80)
    print("\nOutputs:")
    print("  • Metacells: seacells_output/batch_corrected/integrated_metacells_harmony.h5ad")
    print("  • Clustering: seacells_output/clustering_analysis/")
    print("  • VAE Results: metacell_vae_output/")
    print("=" * 80)

if __name__ == "__main__":
    main()
