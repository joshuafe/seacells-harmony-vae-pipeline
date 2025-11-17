"""
Metacell Pipeline v3: Hyperparameter Tuning for Disambiguation VAE

This script constitutes the third step of the corrected analysis pipeline. It:
1. Uses the Optuna framework to run a broad hyperparameter search.
2. Calls the Step 2 script (`pipeline_step2_train_model.py`) for each trial.
3. Optimizes for archetype uniqueness and cluster separation (silhouette score).
"""

import optuna
import numpy as np
import pandas as pd
import scanpy as sc
import subprocess
from pathlib import Path
import json
import warnings
from sklearn.metrics import silhouette_score
import argparse

warnings.filterwarnings('ignore')

def objective(trial, base_args):
    """Optuna objective function for the disambiguation model."""
    
    # Broader search space
    n_archetypes = trial.suggest_int('n_archetypes', 10, 25)
    latent_dim = trial.suggest_int('latent_dim', 8, 20)
    beta = trial.suggest_float('beta', 0.0001, 0.1, log=True)
    gamma = trial.suggest_float('gamma', 0.1, 2.5)
    gamma_cooldown = trial.suggest_int('gamma_cooldown', 50, 250)
    soft_temp = trial.suggest_float('soft_temp', 0.5, 2.5)

    trial_output = Path(base_args['output_dir']) / f"trial_{trial.number}"
    trial_output.mkdir(parents=True, exist_ok=True)

    cmd = [
        'python', 'scripts/pipeline_step2_train_model.py',
        '--input', base_args['input'],
        '--latent-dim', str(latent_dim),
        '--n-archetypes', str(n_archetypes),
        '--beta', str(beta),
        '--gamma', str(gamma),
        '--gamma-cooldown', str(gamma_cooldown),
        '--soft-assignment-temp', str(soft_temp),
        '--epochs', str(base_args['epochs']),
        '--output-dir', str(trial_output)
    ]
    if base_args['train_on_reference']:
        cmd.append('--train-on-reference')

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
        if result.returncode != 0:
            print(f"Trial {trial.number} failed: {result.stderr}")
            raise optuna.TrialPruned()
    except subprocess.TimeoutExpired:
        print(f"Trial {trial.number} timed out")
        raise optuna.TrialPruned()

    try:
        adata = sc.read_h5ad(trial_output / "metacells_with_archetypes.h5ad")
        profiles = pd.read_csv(trial_output / "archetype_profiles_disambiguated.csv", index_col=0)

        ref_mask = adata.obs['is_reference'] == 'True'
        latent = adata[ref_mask].obsm['X_vae']
        labels = adata.obs.loc[ref_mask, 'archetype_id_hard']
        
        silhouette = silhouette_score(latent, labels) if len(np.unique(labels)) > 1 else -1.0
        
        corr_matrix = profiles.T.corr(method='spearman')
        mean_abs_corr = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool)).abs().mean().mean()
        mean_abs_corr = 1.0 if pd.isna(mean_abs_corr) else mean_abs_corr

        score = 0.5 * (1 - (silhouette + 1) / 2) + 0.5 * mean_abs_corr

        trial.set_user_attr('silhouette', float(silhouette))
        trial.set_user_attr('uniqueness', 1.0 - float(mean_abs_corr))
        trial.set_user_attr('avg_entropy', float(adata.obs.loc[ref_mask, 'archetype_entropy'].mean()))
        
        print(f"Trial {trial.number}: score={score:.4f}, silhouette={silhouette:.4f}, uniqueness={1-mean_abs_corr:.4f}")
        return score

    except Exception as e:
        print(f"Trial {trial.number} evaluation failed: {e}")
        raise optuna.TrialPruned()

def main():
    parser = argparse.ArgumentParser(description="Metacell Pipeline v3: Hyperparameter Tuning")
    parser.add_argument("--input", type=str, default="seacells_pipeline_output/integrated_metacells_17_markers.h5ad")
    parser.add_argument("--train-on-reference", action="store_true", default=True)
    parser.add_argument("--epochs", type=int, default=150)
    parser.add_argument("--n-trials", type=int, default=100)
    parser.add_argument("--output-dir", type=str, default="vae_tuning_pipeline")
    args = parser.parse_args()

    print("="*80 + "\nMetacell Pipeline v3: Hyperparameter Tuning\n" + "="*80)
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    base_args = {k: v for k, v in vars(args).items() if k not in ['n_trials', 'output_dir']}
    base_args['output_dir'] = str(output_dir)

    study = optuna.create_study(direction='minimize', study_name='disambiguation_vae_pipeline',
                                storage=f'sqlite:///{output_dir}/tuning_pipeline.db', load_if_exists=True)
    
    study.optimize(lambda trial: objective(trial, base_args), n_trials=args.n_trials)

    # ... (Reporting logic as before) ...
    print("\n" + "="*80 + "\nTUNING COMPLETE\n" + "="*80)
    if study.best_trial:
        print(f"Best trial: {study.best_trial.number} with score {study.best_value:.4f}")
        # ... print best params and metrics ...
    
    results_df = study.trials_dataframe()
    results_df.to_csv(output_dir / "tuning_results.csv", index=False)
    print(f"\n✓ Full tuning results saved to {output_dir / 'tuning_results.csv'}")

if __name__ == "__main__":
    main()
