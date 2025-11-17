"""
Metacell Pipeline v1 (Resumable, Verbose, Robust) - v6 (Alias-Aware)

This v6 script introduces an alias-aware marker mapping function to handle
inconsistencies in column naming between different FCS file batches.
It is also resumable and will not re-process completed samples.
"""

import sys
import numpy as np
import pandas as pd
import scanpy as sc
import SEACells
from pathlib import Path
from fcsparser import parse
import warnings
from multiprocessing import Pool
import time
from datetime import datetime
import argparse
import harmonypy as hm
import phenograph
import psutil
import traceback

warnings.filterwarnings('ignore')
np.random.seed(42)

# --- Ground Truth Marker Panel & Known Aliases ---
CANONICAL_MARKERS = sorted([
    'CD62L', 'CD152 R718', 'CD45', 'CD45RO', 'CD279+CD24', 'CD95', 
    'CD16+TIGIT', 'CD8+CD14', 'CD56', 'CD45RA', 'CD366', 'CD69', 
    'CD25 BB515', 'CD34+CD223', 'CD197', 'CD3+CD19', 'CD4+CD33'
])

MARKER_ALIASES = {
    'CD279+CD24': ['CD24+CD279'],
    'CD152 R718': ['cCD152 R718']
}

def log_progress(message, log_resources=False):
    timestamp = datetime.now().strftime("%H:%M:%S")
    if log_resources:
        cpu = psutil.cpu_percent()
        mem = psutil.virtual_memory().percent
        msg = f"[{timestamp} | CPU: {cpu:5.1f}% | Mem: {mem:5.1f}%] {message}"
    else:
        msg = f"[{timestamp}] {message}"
    print(msg, flush=True)

def get_marker_mapping(fcs_columns):
    """
    Finds the corresponding column name for each canonical marker, checking aliases.
    """
    mapping = {}
    missing_markers = []
    
    for marker in CANONICAL_MARKERS:
        found = False
        # Create a list of names to check: the canonical name first, then any aliases
        names_to_check = [marker] + MARKER_ALIASES.get(marker, [])
        
        for name in names_to_check:
            for col in fcs_columns:
                if name in col:
                    mapping[marker] = col
                    found = True
                    break
            if found:
                break
        
        if not found:
            missing_markers.append(marker)
    
    if missing_markers:
        raise ValueError(f"Missing markers: {', '.join(missing_markers)}")
    return mapping

def process_single_sample(args_tuple):
    fcs_path, sample_id, sample_type, max_cells, compression_ratio, output_dir, log_res = args_tuple
    
    try:
        log_progress(f"Starting {sample_id} ({sample_type})", log_res)
        
        meta, data = parse(str(fcs_path), reformat_meta=False)
        
        try:
            marker_map = get_marker_mapping(data.columns)
        except ValueError as e:
            log_progress(f"✗ {sample_id}: SKIPPING - {e}", log_res)
            return None

        data_subset = data[marker_map.values()].copy()
        data_subset.columns = marker_map.keys()
        
        log_progress(f"  {sample_id}: Loaded {len(data_subset)} cells, {len(data_subset.columns)} markers", log_res)
        
        if len(data_subset) > max_cells:
            data_subset = data_subset.sample(n=max_cells, random_state=42)
            log_progress(f"  {sample_id}: Subsampled to {max_cells} cells", log_res)
            
        n_metacells = max(10, int(len(data_subset) / compression_ratio))
        
        adata = sc.AnnData(data_subset, dtype=np.float32)
        adata.obs['sample_id'] = sample_id
        adata.obs['sample_type'] = sample_type
        
        adata.X = np.arcsinh(adata.X / 5.0)
        sc.pp.scale(adata)
        sc.tl.pca(adata, svd_solver='arpack', n_comps=min(16, len(adata.var_names)-1))
        
        log_progress(f"  {sample_id}: Building SEACells model (n={n_metacells})...", log_res)
        model = SEACells.core.SEACells(adata, build_kernel_on='X_pca', n_SEACells=n_metacells, verbose=False)
        model.construct_kernel_matrix()
        model.fit()
        
        assignment_col_name = 'SEACell_assignment'
        adata.obs[assignment_col_name] = model.get_hard_assignments()['SEACell']
        metacell_ad = SEACells.core.summarize_by_SEACell(adata, SEACells_label=assignment_col_name, summarize_layer='X')
        
        metacell_ad.obs['sample_id'] = sample_id
        metacell_ad.obs['sample_type'] = sample_type
        
        counts = adata.obs[assignment_col_name].value_counts()
        metacell_ad.obs['n_cells_in_metacell'] = metacell_ad.obs.index.map(counts)

        sample_output = output_dir / f"{sample_id}_metacells.h5ad"
        metacell_ad.write_h5ad(sample_output)
        log_progress(f"✓ {sample_id}: Complete!", log_res)
        return str(sample_output)
        
    except Exception as e:
        log_progress(f"✗ {sample_id}: ERROR during processing - {str(e)}", log_res)
        traceback.print_exc()
        return None

def main():
    parser = argparse.ArgumentParser(description="Metacell Pipeline v6 (Alias-Aware)")
    parser.add_argument("--compression-ratio", type=float, default=40.0)
    parser.add_argument("--max-cells", type=int, default=20000)
    parser.add_argument("--n-parallel", type=int, default=4)
    parser.add_argument("--output-dir", type=str, default="seacells_pipeline_output_17_markers")
    parser.add_argument("--phenograph-k", type=int, default=30)
    parser.add_argument("--log-resources", action="store_true", help="Enable CPU and memory logging.")
    args = parser.parse_args()

    log_progress("Starting Metacell Pipeline v6 (Alias-Aware)", args.log_resources)
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    sample_dirs = {"Normal": Path("data/normalsgated"), "Abnormal": Path("data/deidentify_abngated"), "PTCy": Path("data/ptcy_fcs/PTCy Lymphocytes")}
    
    all_fcs_files = []
    for sample_type, directory in sample_dirs.items():
        if directory.exists():
            for i, fcs_path in enumerate(sorted(directory.glob("*.fcs"))):
                all_fcs_files.append({'type': sample_type, 'path': fcs_path, 'id': f"{sample_type}_{i+1}"})
    
    tasks = []
    for fcs_info in all_fcs_files:
        expected_output = output_dir / f"{fcs_info['id']}_metacells.h5ad"
        if expected_output.exists():
            log_progress(f"Skipping {fcs_info['id']}: Output file already exists.", args.log_resources)
        else:
            tasks.append((fcs_info['path'], fcs_info['id'], fcs_info['type'], args.max_cells, args.compression_ratio, output_dir, args.log_resources))

    if tasks:
        log_progress(f"Found {len(tasks)} new samples to process.", args.log_resources)
        start_time = time.time()
        with Pool(processes=args.n_parallel) as pool:
            results = pool.map(process_single_sample, tasks)
        log_progress(f"Sample processing phase complete in {time.time() - start_time:.1f}s", args.log_resources)
    else:
        log_progress("No new samples to process. Moving to integration.", args.log_resources)

    log_progress("Starting Integration Phase...", args.log_resources)
    successful_files = list(output_dir.glob("*_metacells.h5ad"))
    if not successful_files:
        log_progress("ERROR: No metacell files found to integrate. Pipeline halting.", args.log_resources)
        return

    log_progress(f"Found {len(successful_files)} successfully processed samples to integrate.", args.log_resources)
    all_metacells = [sc.read_h5ad(f) for f in successful_files]
    combined = sc.concat(all_metacells, join='inner', merge='same')
    combined.obs_names_make_unique()
    log_progress(f"Combined {combined.n_obs} metacells from {len(successful_files)} samples.", args.log_resources)
    
    sc.pp.scale(combined)
    sc.tl.pca(combined, n_comps=min(16, combined.n_vars - 1))
    
    log_progress("Running Harmony integration...", args.log_resources)
    harmony_out = hm.run_harmony(combined.obsm['X_pca'], combined.obs, 'sample_id', max_iter_harmony=20, verbose=False)
    combined.obsm['X_harmony'] = harmony_out.Z_corr.T
    
    sc.pp.neighbors(combined, n_neighbors=15, use_rep='X_harmony')
    sc.tl.umap(combined)
    
    log_progress(f"Running Phenograph with k={args.phenograph_k}...", args.log_resources)
    communities, _, _ = phenograph.cluster(combined.obsm['X_harmony'], k=args.phenograph_k, n_jobs=-1)
    combined.obs['phenograph_cluster'] = pd.Categorical(communities)
    
    final_output_path = output_dir / "integrated_metacells_17_markers.h5ad"
    combined.write_h5ad(final_output_path)
    log_progress(f"✓ Final analysis-ready file saved to: {final_output_path}", args.log_resources)

if __name__ == "__main__":
    main()
