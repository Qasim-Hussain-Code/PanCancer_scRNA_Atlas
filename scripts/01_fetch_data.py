# ==============================================================================
# PROJECT: Pan-Cancer scRNA-Seq Atlas
# SCRIPT: 01_fetch_data.py
# ==============================================================================

import scanpy as sc
import anndata as ad
import os
import warnings

warnings.filterwarnings('ignore')
sc.settings.verbosity = 3

# --- PATHS ---
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = os.path.join(BASE_DIR, "data", "raw")
sc.settings.writedir = RAW_DIR

print("--- [INITIATING HYBRID DATA CURATION] ---")

try:
    # -------------------------------------------------------------------------
    # DATASET 1: Healthy Immune Baseline (Scanpy API)
    # -------------------------------------------------------------------------
    print("\n[1/3] Fetching Cohort A: Healthy Immune Baseline...")
    adata_immune = sc.datasets.pbmc3k()
    adata_immune.obs['cohort'] = 'Immune_Baseline'
    adata_immune.obs['cancer_type'] = 'Healthy_Reference'
    
    immune_path = os.path.join(RAW_DIR, "cohort_A_healthy.h5ad")
    adata_immune.write(immune_path)
    print(f"      -> Curated: {adata_immune.n_obs} baseline cells.")

    # -------------------------------------------------------------------------
    # DATASET 2: Leukemia / Blood Cancer (Scanpy API)
    # -------------------------------------------------------------------------
    print("\n[2/3] Fetching Cohort B: Myeloid Leukemia...")
    adata_leukemia = sc.datasets.paul15()
    
    adata_leukemia.obs['cohort'] = 'Tumor_Site_1'
    adata_leukemia.obs['cancer_type'] = 'Leukemia'
    adata_leukemia.X = adata_leukemia.X.astype('float32')
    
    leuk_path = os.path.join(RAW_DIR, "cohort_B_leukemia.h5ad")
    adata_leukemia.write(leuk_path)
    print(f"      -> Curated: {adata_leukemia.n_obs} leukemia cells.")

    # -------------------------------------------------------------------------
    # DATASET 3: Tumor Immune Cell Atlas (Local Ingestion & Rebuild)
    # -------------------------------------------------------------------------
    print("\n[3/3] Processing Cohort C: Tumor Immune Cell Atlas (Local)...")
    
    tic_file = os.path.join(RAW_DIR, "TICAtlas_downsampled.h5ad")
    
    if not os.path.exists(tic_file):
        raise FileNotFoundError(f"\n[!] Missing {tic_file}.\nPlease ensure you downloaded the file from Zenodo and placed it in the data/raw/ folder.")
        
    adata_tumor = sc.read_h5ad(tic_file)
    
    # --- THE SANITIZER: Deep Memory Rebuild ---
    # 1. Strip the reserved '_index' name if it exists in standard columns
    if '_index' in adata_tumor.obs.columns:
        adata_tumor.obs = adata_tumor.obs.rename(columns={'_index': 'original_index'})
    if '_index' in adata_tumor.var.columns:
        adata_tumor.var = adata_tumor.var.rename(columns={'_index': 'original_index'})
        
    # 2. Re-instantiate the object from raw parts. 
    # This completely deletes the corrupted R-Seurat '.raw' dictionary.
    adata_tumor = ad.AnnData(
        X=adata_tumor.X.copy(), 
        obs=adata_tumor.obs.copy(), 
        var=adata_tumor.var.copy()
    )
    
    # 3. Clear Pandas index names
    adata_tumor.obs.index.name = None
    adata_tumor.var.index.name = None
    # ------------------------------------------
    
    # Subsample to keep the Pan-Cancer Atlas mathematically balanced
    if adata_tumor.n_obs > 3000:
        sc.pp.subsample(adata_tumor, n_obs=3000, random_state=42)
        
    adata_tumor.obs['cohort'] = 'Tumor_Site_2'
    adata_tumor.obs['cancer_type'] = 'Solid_Tumor_TME'
    
    tumor_path = os.path.join(RAW_DIR, "cohort_C_solid_tumor.h5ad")
    adata_tumor.write(tumor_path)
    print(f"      -> Curated: {adata_tumor.n_obs} solid tumor cells.")

    # -------------------------------------------------------------------------
    # SUMMARY
    # -------------------------------------------------------------------------
    total_cells = adata_immune.n_obs + adata_leukemia.n_obs + adata_tumor.n_obs
    print("\n========================================================")
    print(f"CURATION COMPLETE: {total_cells} Single Cells Processed.")
    print("Files successfully formatted and staged in: data/raw/")
    print("========================================================")

except Exception as e:
    print(f"\n[FATAL ERROR] Pipeline halted. Error: {e}")