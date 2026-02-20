# ==============================================================================
# PROJECT: Pan-Cancer scRNA-Seq Atlas
# SCRIPT: 02_build_end_level_atlas.py
# DESCRIPTION: QC, Normalization, scVI Deep Learning Integration, and UMAP.
# ==============================================================================

import scanpy as sc
import anndata as ad
import scvi
import os
import warnings

warnings.filterwarnings('ignore')
sc.settings.verbosity = 3

# --- PATHS ---
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = os.path.join(BASE_DIR, "data", "raw")
CURATED_DIR = os.path.join(BASE_DIR, "data", "curated")

print("--- [INITIATING AI INTEGRATION PIPELINE] ---")

# 1. LOAD & CONCATENATE
print("\n[1/6] Loading Curated Cohorts...")
adata_A = sc.read_h5ad(os.path.join(RAW_DIR, "cohort_A_healthy.h5ad"))
adata_B = sc.read_h5ad(os.path.join(RAW_DIR, "cohort_B_leukemia.h5ad"))
adata_C = sc.read_h5ad(os.path.join(RAW_DIR, "cohort_C_solid_tumor.h5ad"))

# Merge into one massive Atlas
adata = ad.concat([adata_A, adata_B, adata_C], label="batch")
adata.obs_names_make_unique()
print(f"      -> Pre-QC Atlas constructed: {adata.n_obs} cells.")

# Preserve raw counts for Deep Learning and later R/Seurat import
adata.layers["counts"] = adata.X.copy()

# 2. QUALITY CONTROL
print("\n[2/6] Executing Biological Quality Control...")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print(f"      -> Post-QC Atlas: {adata.n_obs} cells x {adata.n_vars} genes.")

# 3. NORMALIZATION & FEATURE SELECTION
print("\n[3/6] Normalizing and Isolating Highly Variable Genes...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True, layer="counts", flavor="seurat_v3", batch_key="batch")
print("      -> Locked in top 2000 highly variable features.")

# 4. DEEP LEARNING INTEGRATION (scVI)
print("\n[4/6] Booting Variational Autoencoder (scVI)...")
print("      (Training the neural network. This will take a few minutes...)")

scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
# Initialize the PyTorch model
model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="zinb")

# Train the AI (Using 50 epochs for local hardware efficiency)
model.train(max_epochs=50, early_stopping=True)

# Extract the mathematically perfectly aligned cells
adata.obsm["X_scVI"] = model.get_latent_representation()
print("      -> AI Integration Complete. Batch effects surgically removed.")

# 5. ECOLOGICAL MANIFOLD (UMAP & Clustering)
print("\n[5/6] Computing 3D UMAP and Leiden Clusters...")
# We use the AI's math, NOT standard PCA, to draw the map
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15)
sc.tl.umap(adata, min_dist=0.3)
sc.tl.leiden(adata, resolution=0.8, flavor="igraph")
print("      -> Spatial manifold coordinates locked.")

# 6. EXPORT THE ATLAS
print("\n[6/6] Saving the Final Grandmaster Atlas...")
FINAL_PATH = os.path.join(CURATED_DIR, "PanCancer_scAtlas_EndLevel.h5ad")
adata.write(FINAL_PATH)

print("\n========================================================")
print(f"SUCCESS: The End-Level Dataset is forged.")
print(f"Location: {FINAL_PATH}")
print("========================================================")