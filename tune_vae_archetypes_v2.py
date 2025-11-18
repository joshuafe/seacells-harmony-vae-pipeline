"""
Hyperparameter tuning for VAE archetypes using Optuna (v2)

This script runs a grid/random search over key hyperparameters.
The objective function is redesigned to prioritize:
- Archetype uniqueness (low correlation between profiles)
- Clear separation in latent space (high silhouette score)
It no longer penalizes high entropy.
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
    Optuna objective function (v2)

    We want to optimize for:
    1. Unique archetypes (low correlation between marker profiles)
    2. Clear separation (high silhouette score in latent space)
    """

    # Sample hyperparameters from a broader search space
    n_archetypes = trial.suggest_int('n_archetypes', 8, 20)
    latent_dim = trial.suggest_int('latent_dim', 4, 12)
    beta = trial.suggest_float('beta', 0.0005, 0.1, log=True)
    gamma = trial.suggest_float('gamma', 0.0, 1.0)
    gamma_cooldown = trial.suggest_int('gamma_cooldown', 50, 250)
    gamma_min = trial.suggest_float('gamma_min', 0.0, 0.1)
    soft_temp = trial.suggest_float('soft_temp', 0.5, 2.5)

    # Create trial-specific output directory
    trial_output = Path(base_args['output_dir']) / f"trial_{trial.number}"
    trial_output.mkdir(parents=True, exist_ok=True)

    # Build command
    cmd = [
        'python',
        'scripts/vae_archetypes_with_priors.py',
        '--input',
        base_args['input'],
        '--latent-dim',
        str(latent_dim),
        '--n-archetypes',
        str(n_archetypes),
        '--beta',
        str(beta),
        '--gamma',
        str(gamma),
        '--gamma-cooldown',
        str(gamma_cooldown),
        '--gamma-min',
        str(gamma_min),
        '--soft-assignment-temp',
        str(soft_temp),
        '--epochs',
        str(base_args['epochs']),
        '--output-dir',
        str(trial_output)
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
        profiles = pd.read_csv(trial_output / "archetype_profiles.csv", index_col=0)

        # --- Metric 1: Silhouette score in latent space (higher is better) ---
        ref_mask = adata.obs.get('is_reference', True)
        latent = adata[ref_mask].obsm['X_vae']
        labels = adata.obs.loc[ref_mask, 'archetype_id_hard']
        if len(np.unique(labels)) > 1:
            silhouette = silhouette_score(latent, labels, metric='euclidean')
        else:
            silhouette = 0.0

        # --- Metric 2: Archetype uniqueness (lower correlation is better) ---
        marker_cols = [c for c in profiles.columns if c not in ['gini', 'entropy', 'n_metacells']]
        corr_matrix = profiles[marker_cols].T.corr(method='spearman')
        upper_tri = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
        mean_abs_corr = upper_tri.abs().mean().mean()
        
        # Handle case where there's only 1 archetype
        if pd.isna(mean_abs_corr):
            mean_abs_corr = 1.0 # Worst possible score

        # --- Combined score (Optuna minimizes) ---
        # We want to MAXIMIZE silhouette and MINIMIZE mean_abs_corr
        # So we minimize (1 - silhouette) and mean_abs_corr
        # Map silhouette from [-1, 1] to [0, 1] for the score
        normalized_silhouette_for_max = (silhouette + 1) / 2
        
        # Objective: 50% weight on separation, 50% on uniqueness
        score = 0.5 * (1 - normalized_silhouette_for_max) + 0.5 * mean_abs_corr

        # Store metrics for reporting
        trial.set_user_attr('silhouette', float(silhouette))
        trial.set_user_attr('mean_abs_corr', float(mean_abs_corr))
        trial.set_user_attr('avg_entropy', float(adata.obs.loc[ref_mask, 'archetype_entropy'].mean()))
        trial.set_user_attr('output_dir', str(trial_output))

        print(f"Trial {trial.number}: score={score:.4f}, silhouette={silhouette:.4f}, uniqueness={1-mean_abs_corr:.4f}")

        return score

    except Exception as e:
        print(f"Trial {trial.number} evaluation failed: {e}")
        raise optuna.TrialPruned()


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Hyperparameter tuning for VAE archetypes (v2)")
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
        default=150, # Reduced from 200 for faster tuning
        help="Training epochs per trial"
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
        default="optuna_tuning_v2",
        help="Base output directory for tuning"
    )

    args = parser.parse_args()

    print("="*80)
    print("VAE ARCHETYPE HYPERPARAMETER TUNING (v2)")
    print("="*80)
    print(f"Input: {args.input}")
    print(f"Train on reference: {args.train_on_reference}")
    print(f"Number of trials: {args.n_trials}")
    print(f"Epochs per trial: {args.epochs}")
    print(f"Output directory: {args.output_dir}")
    print("Objective: Maximize Silhouette Score & Archetype Uniqueness")
    print("="*80)

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    base_args = {
        'input': args.input,
        'train_on_reference': args.train_on_reference,
        'epochs': args.epochs,
        'output_dir': str(output_dir)
    }

    study = optuna.create_study(
        direction='minimize',
        study_name='vae_archetypes_v2',
        storage=f'sqlite:///{output_dir}/optuna_v2.db',
        load_if_exists=True
    )

    study.optimize(
        lambda trial: objective(trial, base_args),
        n_trials=args.n_trials,
        timeout=None,
        show_progress_bar=True
    )

    print("\n" + "="*80)
    print("TUNING RESULTS")
    print("="*80)
    print(f"Best trial: {study.best_trial.number}")
    print(f"Best score: {study.best_value:.4f}")
    print(f"\nBest hyperparameters:")
    for key, value in study.best_trial.params.items():
        print(f"  {key}: {value}")

    print(f"\nBest metrics:")
    print(f"  Silhouette score: {study.best_trial.user_attrs['silhouette']:.4f}")
    print(f"  Archetype Uniqueness (1-corr): {1 - study.best_trial.user_attrs['mean_abs_corr']:.4f}")
    print(f"  (for reference) Avg entropy: {study.best_trial.user_attrs['avg_entropy']:.4f}")
    print(f"  Output: {study.best_trial.user_attrs['output_dir']}")

    results_df = study.trials_dataframe()
    results_df.to_csv(output_dir / "tuning_results_v2.csv", index=False)
    print(f"\n✓ Saved results to {output_dir / 'tuning_results_v2.csv'}")

    best_config = {
        'params': study.best_trial.params,
        'metrics': {
            'score': study.best_value,
            'silhouette': study.best_trial.user_attrs['silhouette'],
            'mean_abs_corr': study.best_trial.user_attrs['mean_abs_corr'],
            'avg_entropy': study.best_trial.user_attrs['avg_entropy'],
        }
    }

    with open(output_dir / "best_config_v2.json", 'w') as f:
        json.dump(best_config, f, indent=2)

    print(f"✓ Saved best config to {output_dir / 'best_config_v2.json'}")
    print("="*80)


if __name__ == "__main__":
    main()
