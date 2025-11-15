"""
Hyperparameter tuning for VAE archetypes using Optuna

This script runs a grid/random search over key hyperparameters:
- Number of archetypes
- Latent dimensions
- Beta (KL weight)
- Gamma (biological prior weight)
- Soft assignment temperature
"""

import optuna
import numpy as np
import pandas as pd
import scanpy as sc
import torch
import subprocess
from pathlib import Path
import json
import warnings
warnings.filterwarnings('ignore')


def objective(trial, base_args):
    """
    Optuna objective function

    We want to optimize for:
    1. Good reconstruction (low reconstruction loss)
    2. Biologically meaningful archetypes (low entropy on reference samples)
    3. Clear separation (high silhouette score in latent space)
    """

    # Sample hyperparameters
    n_archetypes = trial.suggest_int('n_archetypes', 8, 15)
    latent_dim = trial.suggest_int('latent_dim', 4, 8)
    beta = trial.suggest_float('beta', 0.001, 0.1, log=True)
    gamma = trial.suggest_float('gamma', 0.1, 1.0)
    gamma_cooldown = trial.suggest_int('gamma_cooldown', 0, 200)
    gamma_min = trial.suggest_float('gamma_min', 0.0, 0.1)
    soft_temp = trial.suggest_float('soft_temp', 0.5, 2.0)

    # Create trial-specific output directory
    trial_output = Path(base_args['output_dir']) / f"trial_{trial.number}"
    trial_output.mkdir(parents=True, exist_ok=True)

    # Build command
    cmd = [
        'python', 'scripts/vae_archetypes_with_priors.py',
        '--input', base_args['input'],
        '--latent-dim', str(latent_dim),
        '--n-archetypes', str(n_archetypes),
        '--beta', str(beta),
        '--gamma', str(gamma),
        '--gamma-cooldown', str(gamma_cooldown),
        '--gamma-min', str(gamma_min),
        '--soft-assignment-temp', str(soft_temp),
        '--epochs', str(base_args['epochs']),
        '--output-dir', str(trial_output)
    ]

    if base_args['train_on_reference']:
        cmd.append('--train-on-reference')

    # Run training
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600  # 10 min timeout
        )

        if result.returncode != 0:
            print(f"Trial {trial.number} failed: {result.stderr}")
            raise optuna.TrialPruned()

    except subprocess.TimeoutExpired:
        print(f"Trial {trial.number} timed out")
        raise optuna.TrialPruned()

    # Load results and compute metrics
    try:
        adata = sc.read_h5ad(trial_output / "metacells_with_archetypes.h5ad")

        # Metric 1: Average entropy (lower = more clear assignments)
        # Focus on reference samples
        ref_mask = adata.obs.get('is_reference', True)
        avg_entropy = adata.obs.loc[ref_mask, 'archetype_entropy'].mean()

        # Metric 2: Silhouette score in latent space
        from sklearn.metrics import silhouette_score
        latent = adata[ref_mask].obsm['X_vae']
        labels = adata.obs.loc[ref_mask, 'archetype_id_hard']

        if len(np.unique(labels)) > 1:
            silhouette = silhouette_score(latent, labels, metric='euclidean')
        else:
            silhouette = 0.0

        # Metric 3: Reconstruction quality (load from training)
        # We'll use the final reconstruction loss
        # For now, use entropy and silhouette

        # Combined score (lower is better for entropy, higher for silhouette)
        # Normalize to [0, 1] range
        # Typical entropy range: 0-2.5, silhouette: -1 to 1
        normalized_entropy = avg_entropy / 2.5
        normalized_silhouette = (silhouette + 1) / 2  # map to [0, 1]

        # Objective: minimize entropy, maximize silhouette
        score = normalized_entropy - 0.5 * normalized_silhouette

        # Store metrics for reporting
        trial.set_user_attr('avg_entropy', avg_entropy)
        trial.set_user_attr('silhouette', silhouette)
        trial.set_user_attr('output_dir', str(trial_output))

        print(f"Trial {trial.number}: score={score:.4f}, entropy={avg_entropy:.4f}, silhouette={silhouette:.4f}")

        return score

    except Exception as e:
        print(f"Trial {trial.number} evaluation failed: {e}")
        raise optuna.TrialPruned()


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Hyperparameter tuning for VAE archetypes")
    parser.add_argument(
        "--input",
        type=str,
        default="seacells_output/three_groups_phenograph/integrated_metacells_harmony_phenograph.h5ad",
        help="Input h5ad file"
    )
    parser.add_argument(
        "--train-on-reference",
        action="store_true",
        help="Train only on reference samples"
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=200,
        help="Training epochs (reduced for tuning)"
    )
    parser.add_argument(
        "--n-trials",
        type=int,
        default=50,
        help="Number of hyperparameter trials"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="optuna_vae_tuning",
        help="Base output directory for tuning"
    )

    args = parser.parse_args()

    print("="*80)
    print("VAE ARCHETYPE HYPERPARAMETER TUNING")
    print("="*80)
    print(f"Input: {args.input}")
    print(f"Train on reference: {args.train_on_reference}")
    print(f"Number of trials: {args.n_trials}")
    print(f"Epochs per trial: {args.epochs}")
    print("="*80)

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Base arguments for trials
    base_args = {
        'input': args.input,
        'train_on_reference': args.train_on_reference,
        'epochs': args.epochs,
        'output_dir': str(output_dir)
    }

    # Create Optuna study
    study = optuna.create_study(
        direction='minimize',
        study_name='vae_archetypes',
        storage=f'sqlite:///{output_dir}/optuna.db',
        load_if_exists=True
    )

    # Run optimization
    study.optimize(
        lambda trial: objective(trial, base_args),
        n_trials=args.n_trials,
        timeout=None,
        show_progress_bar=True
    )

    # Report results
    print("\n" + "="*80)
    print("TUNING RESULTS")
    print("="*80)
    print(f"Best trial: {study.best_trial.number}")
    print(f"Best score: {study.best_value:.4f}")
    print(f"\nBest hyperparameters:")
    for key, value in study.best_trial.params.items():
        print(f"  {key}: {value}")

    print(f"\nBest metrics:")
    print(f"  Average entropy: {study.best_trial.user_attrs['avg_entropy']:.4f}")
    print(f"  Silhouette score: {study.best_trial.user_attrs['silhouette']:.4f}")
    print(f"  Output: {study.best_trial.user_attrs['output_dir']}")

    # Save results
    results_df = study.trials_dataframe()
    results_df.to_csv(output_dir / "tuning_results.csv", index=False)
    print(f"\n✓ Saved results to {output_dir / 'tuning_results.csv'}")

    # Save best config
    best_config = {
        'params': study.best_trial.params,
        'metrics': {
            'score': study.best_value,
            'avg_entropy': study.best_trial.user_attrs['avg_entropy'],
            'silhouette': study.best_trial.user_attrs['silhouette']
        }
    }

    with open(output_dir / "best_config.json", 'w') as f:
        json.dump(best_config, f, indent=2)

    print(f"✓ Saved best config to {output_dir / 'best_config.json'}")
    print("="*80)


if __name__ == "__main__":
    main()
