"""
Per-Sample SEACells with Harmony Batch Correction
Following the recommended workflow from Persad et al. Nature Biotechnology 2023
"""

import sys
import numpy as np
import pandas as pd
import scanpy as sc
import SEACells
from pathlib import Path
from fcsparser import parse
import warnings
warnings.filterwarnings('ignore')
from multiprocessing import Pool, Manager
import time
from datetime import datetime

# Set random seed for reproducibility
np.random.seed(42)

# Markers to use (same 18 as before)
MARKERS = [
    'CD45', 'CD3+CD19', 'CD4+CD33', 'CD8+CD14', 'CD16+TIGIT',
    'CD34+CD223', 'CD45RA', 'CD45RO', 'CD56', 'CD69',
    'CD95', 'CD152', 'CD197', 'CD278', 'CD366', 'HLA-DR', 'Ki-67', 'PD-1'
]

def log_progress(message, shared_list=None):
    """Thread-safe logging with timestamp."""
    timestamp = datetime.now().strftime("%H:%M:%S")
    msg = f"[{timestamp}] {message}"
    print(msg, flush=True)
    if shared_list is not None:
        shared_list.append(msg)

def process_single_sample(args):
    """
    Process a single sample with SEACells.
    Returns path to saved metacell h5ad file.
    """
    fcs_path, sample_id, sample_type, max_cells, compression_ratio, output_dir, shared_log = args

    try:
        log_progress(f"Starting {sample_id} ({sample_type})", shared_log)

        # Parse FCS file
        log_progress(f"  {sample_id}: Loading FCS file...", shared_log)
        meta, data = parse(str(fcs_path), reformat_meta=False)

        # Subset markers
        available_markers = [m for m in MARKERS if m in data.columns]
        data_subset = data[available_markers].copy()

        log_progress(f"  {sample_id}: Loaded {len(data_subset)} cells, {len(available_markers)} markers", shared_log)

        # Subsample if needed (cap at max_cells)
        if len(data_subset) > max_cells:
            sample_indices = np.random.choice(len(data_subset), size=max_cells, replace=False)
            data_subset = data_subset.iloc[sample_indices]
            log_progress(f"  {sample_id}: Subsampled to {max_cells} cells", shared_log)

        # Calculate n_metacells based on actual cell count and compression ratio
        n_cells_actual = len(data_subset)
        n_metacells = max(10, int(n_cells_actual / compression_ratio))  # Minimum 10 metacells
        log_progress(f"  {sample_id}: Will generate {n_metacells} metacells ({compression_ratio:.1f}:1 compression)", shared_log)

        # Create AnnData
        adata = sc.AnnData(X=data_subset.values,
                          var=pd.DataFrame(index=available_markers),
                          obs=pd.DataFrame(index=[f"{sample_id}_cell_{i}" for i in range(len(data_subset))]))
        adata.var_names = available_markers

        # Add metadata
        adata.obs['sample_id'] = sample_id
        adata.obs['sample_type'] = sample_type

        log_progress(f"  {sample_id}: Preprocessing...", shared_log)
        # Arcsinh transformation
        adata.X = np.arcsinh(adata.X / 5.0)

        # Z-score scaling
        sc.pp.scale(adata)

        # PCA
        sc.tl.pca(adata, svd_solver='arpack', n_comps=min(17, len(available_markers)-1))

        # Compute neighbors (required by SEACells)
        sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca')

        log_progress(f"  {sample_id}: Building SEACells model (n={n_metacells})...", shared_log)
        # Build SEACells model
        model = SEACells.core.SEACells(
            adata,
            build_kernel_on='X_pca',
            n_SEACells=n_metacells,
            n_waypoint_eigs=10,
            convergence_epsilon=1e-5,
            verbose=False
        )

        # Construct kernel and fit
        log_progress(f"  {sample_id}: Constructing kernel matrix...", shared_log)
        model.construct_kernel_matrix()

        log_progress(f"  {sample_id}: Fitting SEACells...", shared_log)
        model.fit()

        log_progress(f"  {sample_id}: SEACells converged", shared_log)

        # Get metacell assignments
        adata.obs['SEACell'] = model.get_hard_assignments()['SEACell']

        # Create metacell-level data
        log_progress(f"  {sample_id}: Aggregating metacells...", shared_log)
        metacell_ad = SEACells.core.summarize_by_SEACell(
            adata,
            SEACells_label='SEACell',
            summarize_layer='X',
            celltype_label=None
        )

        # Add sample metadata
        metacell_ad.obs['sample_id'] = sample_id
        metacell_ad.obs['sample_type'] = sample_type
        metacell_ad.obs['n_cells_in_metacell'] = adata.obs.groupby('SEACell').size().values

        # Save individual metacell h5ad
        sample_output = output_dir / f"{sample_id}_metacells.h5ad"
        metacell_ad.write_h5ad(sample_output)

        log_progress(f"✓ {sample_id}: Complete! Generated {metacell_ad.n_obs} metacells → {sample_output.name}", shared_log)

        return str(sample_output)

    except Exception as e:
        error_msg = f"✗ {sample_id}: ERROR - {str(e)}"
        log_progress(error_msg, shared_log)
        return None

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Per-sample SEACells with Harmony batch correction"
    )
    parser.add_argument(
        "--compression-ratio",
        type=float,
        default=50.0,
        help="Compression ratio (cells per metacell). Lower = more metacells. Default: 50.0"
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=10000,
        help="Maximum cells to use per sample. Default: 10000"
    )
    parser.add_argument(
        "--n-parallel",
        type=int,
        default=4,
        help="Number of samples to process in parallel. Default: 4"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="seacells_output/batch_corrected",
        help="Output directory. Default: seacells_output/batch_corrected"
    )

    args = parser.parse_args()

    print("="*80)
    print("PER-SAMPLE SEACELLS + HARMONY BATCH CORRECTION")
    print("="*80)
    print("Following Persad et al. (2023) workflow:")
    print("  1. Run SEACells independently on each sample")
    print("  2. Combine all metacells")
    print("  3. Apply Harmony batch correction")
    print("="*80)

    # Configuration from arguments
    max_cells_per_sample = args.max_cells
    compression_ratio = args.compression_ratio
    n_parallel = args.n_parallel

    # Output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Sample selection
    normal_dir = Path("data/normalsgated")
    abnormal_dir = Path("data/deidentify_abngated")

    # Get ALL normal samples
    normal_files = sorted(normal_dir.glob("*.fcs"))

    # Get ALL abnormal samples
    abnormal_files = sorted(abnormal_dir.glob("*.fcs"))

    print(f"\nConfiguration:")
    print(f"  Normal samples: {len(normal_files)}")
    print(f"  Abnormal samples: {len(abnormal_files)}")
    print(f"  Max cells per sample: {max_cells_per_sample:,}")
    print(f"  Compression ratio: {compression_ratio:.1f}:1 (constant)")
    print(f"  Expected metacells per sample: ~{int(max_cells_per_sample / compression_ratio)}")
    print(f"  Parallel processes: {n_parallel}")
    print(f"  Output: {output_dir}")
    print()
    print("NOTE: Each sample will generate metacells based on its actual cell count")
    print("      to maintain constant compression ratio across all samples.")
    print()

    # Prepare arguments for parallel processing
    manager = Manager()
    shared_log = manager.list()

    tasks = []

    # Add normal samples
    for i, fcs_path in enumerate(normal_files):
        sample_id = f"Normal_{i+1}"
        tasks.append((fcs_path, sample_id, "Normal", max_cells_per_sample,
                     compression_ratio, output_dir, shared_log))

    # Add abnormal samples
    for i, fcs_path in enumerate(abnormal_files):
        sample_id = f"Abnormal_{i+1}"
        tasks.append((fcs_path, sample_id, "Abnormal", max_cells_per_sample,
                     compression_ratio, output_dir, shared_log))

    print("="*80)
    print(f"PHASE 1: Running SEACells on {len(tasks)} samples in parallel")
    print("="*80)
    start_time = time.time()

    # Run parallel processing
    with Pool(processes=n_parallel) as pool:
        results = pool.map(process_single_sample, tasks)

    # Filter successful results
    successful_files = [r for r in results if r is not None]

    elapsed = time.time() - start_time
    print()
    print("="*80)
    print(f"PHASE 1 COMPLETE: {len(successful_files)}/{len(tasks)} samples processed")
    print(f"Total time: {elapsed/60:.1f} minutes")
    print("="*80)

    if len(successful_files) == 0:
        print("ERROR: No samples processed successfully!")
        return

    # Phase 2: Combine metacells
    print("\n" + "="*80)
    print("PHASE 2: Combining all metacells")
    print("="*80)

    all_metacells = []
    for h5ad_path in successful_files:
        adata = sc.read_h5ad(h5ad_path)
        all_metacells.append(adata)
        print(f"  Loaded {adata.n_obs} metacells from {Path(h5ad_path).name}")

    # Concatenate
    combined = sc.concat(all_metacells, join='outer', merge='same')
    combined.obs_names_make_unique()

    print(f"\nCombined: {combined.n_obs} metacells, {combined.n_vars} markers")
    print(f"  Normal metacells: {(combined.obs['sample_type'] == 'Normal').sum()}")
    print(f"  Abnormal metacells: {(combined.obs['sample_type'] == 'Abnormal').sum()}")

    # Save pre-correction
    combined.write_h5ad(output_dir / "combined_metacells_pre_harmony.h5ad")
    print(f"\nSaved: combined_metacells_pre_harmony.h5ad")

    # Phase 3: Harmony batch correction
    print("\n" + "="*80)
    print("PHASE 3: Applying Harmony batch correction")
    print("="*80)

    import harmonypy as hm

    # PCA on combined metacells
    print("  Running PCA on combined metacells...")
    sc.pp.scale(combined)
    n_pca_comps = min(12, combined.n_vars - 1, combined.n_obs - 1)
    print(f"  Using {n_pca_comps} PCA components (max for {combined.n_vars} markers)")
    sc.tl.pca(combined, n_comps=n_pca_comps)

    # Run Harmony
    print("  Running Harmony integration...")
    harmony_out = hm.run_harmony(
        combined.obsm['X_pca'],
        combined.obs,
        'sample_id',
        max_iter_harmony=20,
        verbose=False
    )

    combined.obsm['X_harmony'] = harmony_out.Z_corr.T
    print("  ✓ Harmony correction complete")

    # Compute neighbors and UMAP on Harmony-corrected space
    print("  Computing neighbors on Harmony-corrected space...")
    sc.pp.neighbors(combined, n_neighbors=15, use_rep='X_harmony')

    print("  Computing UMAP...")
    sc.tl.umap(combined)

    # Save final integrated data
    combined.write_h5ad(output_dir / "integrated_metacells_harmony.h5ad")
    print(f"\n✓ Saved: integrated_metacells_harmony.h5ad")

    # Phase 4: Generate visualizations
    print("\n" + "="*80)
    print("PHASE 4: Generating visualizations")
    print("="*80)

    import matplotlib.pyplot as plt
    import seaborn as sns

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    # 1. UMAP by sample type
    ax = axes[0, 0]
    for sample_type in ['Normal', 'Abnormal']:
        mask = combined.obs['sample_type'] == sample_type
        color = '#1f77b4' if sample_type == 'Normal' else '#ff7f0e'
        ax.scatter(combined.obsm['X_umap'][mask, 0],
                  combined.obsm['X_umap'][mask, 1],
                  c=color, label=sample_type, s=20, alpha=0.6, edgecolors='none')
    ax.set_xlabel('UMAP 1', fontweight='bold')
    ax.set_ylabel('UMAP 2', fontweight='bold')
    ax.set_title('Metacells by Sample Type\n(Harmony-Corrected)', fontweight='bold', fontsize=12)
    ax.legend(frameon=True, loc='best')
    ax.grid(True, alpha=0.3)

    # 2. UMAP by sample ID
    ax = axes[0, 1]
    sample_ids = combined.obs['sample_id'].unique()
    colors = plt.cm.tab20(np.linspace(0, 1, len(sample_ids)))
    for i, sample_id in enumerate(sample_ids):
        mask = combined.obs['sample_id'] == sample_id
        ax.scatter(combined.obsm['X_umap'][mask, 0],
                  combined.obsm['X_umap'][mask, 1],
                  c=[colors[i]], label=sample_id, s=15, alpha=0.5, edgecolors='none')
    ax.set_xlabel('UMAP 1', fontweight='bold')
    ax.set_ylabel('UMAP 2', fontweight='bold')
    ax.set_title('Metacells by Sample ID\n(Showing Batch Correction)', fontweight='bold', fontsize=12)
    ax.legend(frameon=True, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=7, ncol=1)
    ax.grid(True, alpha=0.3)

    # 3. Sample composition
    ax = axes[1, 0]
    composition = combined.obs.groupby(['sample_type', 'sample_id']).size().unstack(fill_value=0)
    composition.plot(kind='bar', stacked=False, ax=ax, width=0.7)
    ax.set_ylabel('Number of Metacells', fontweight='bold')
    ax.set_xlabel('Sample Type', fontweight='bold')
    ax.set_title('Metacells per Sample', fontweight='bold', fontsize=12)
    ax.legend(title='Sample ID', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=7)
    ax.grid(axis='y', alpha=0.3)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)

    # 4. Cells per metacell distribution
    ax = axes[1, 1]
    combined.obs.boxplot(column='n_cells_in_metacell', by='sample_type', ax=ax)
    ax.set_ylabel('Cells per Metacell', fontweight='bold')
    ax.set_xlabel('Sample Type', fontweight='bold')
    ax.set_title('Distribution of Metacell Sizes', fontweight='bold', fontsize=12)
    ax.grid(axis='y', alpha=0.3)
    plt.suptitle('')  # Remove automatic title

    plt.tight_layout()
    plt.savefig(output_dir / "harmony_integration_overview.png", dpi=300, bbox_inches='tight')
    print("  ✓ Saved: harmony_integration_overview.png")
    plt.close()

    # Summary statistics
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Total metacells: {combined.n_obs}")
    print(f"  Normal: {(combined.obs['sample_type'] == 'Normal').sum()}")
    print(f"  Abnormal: {(combined.obs['sample_type'] == 'Abnormal').sum()}")
    print(f"Markers: {combined.n_vars}")
    print(f"Samples: {len(combined.obs['sample_id'].unique())}")
    print(f"\nMean cells per metacell: {combined.obs['n_cells_in_metacell'].mean():.1f}")
    print(f"Median cells per metacell: {combined.obs['n_cells_in_metacell'].median():.1f}")
    print(f"\nTotal processing time: {(time.time() - start_time)/60:.1f} minutes")
    print("="*80)

    print("\nNext steps:")
    print("  1. Run clustering on integrated_metacells_harmony.h5ad")
    print("  2. Run Random Forest analysis")
    print("  3. Generate marker rank-lists")
    print()

if __name__ == "__main__":
    main()
