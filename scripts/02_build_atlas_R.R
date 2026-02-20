# ==============================================================================
# PROJECT: Pan-Cancer scRNA-Seq Atlas
# SCRIPT: 02_build_atlas_R.R
# DESCRIPTION: Seurat v5 Integration and Visualization
# ==============================================================================

library(Seurat)
library(harmony) 
library(ggplot2)
library(anndata) 

# --- PATHS ---
# Adjust this path if you are not running it from the project root
raw_dir <- "data/raw/"
out_dir <- "results/"

print("--- [INITIATING R-BASED INTEGRATION] ---")

# 2. LOAD DATA
load_data <- function(file, name) {
  obj_data <- read_h5ad(file)
  obj <- CreateSeuratObject(counts = t(obj_data$X), project = name)
  obj$cancer_type <- name
  return(obj)
}

immune <- load_data(paste0(raw_dir, "cohort_A_healthy.h5ad"), "Healthy")
leukemia <- load_data(paste0(raw_dir, "cohort_B_leukemia.h5ad"), "Leukemia")
tumor <- load_data(paste0(raw_dir, "cohort_C_solid_tumor.h5ad"), "Solid_Tumor")

# 3. MERGE & PROCESS
print("[1/3] Merging Cohorts...")
atlas <- merge(immune, y = c(leukemia, tumor), add.cell.ids = c("H", "L", "S"))

# --- THE SEURAT V5 FIX ---
# Join the distinct layers back into a single matrix for unified processing
atlas <- JoinLayers(atlas)
# -------------------------

# Standard Seurat Workflow
atlas <- NormalizeData(atlas)
atlas <- FindVariableFeatures(atlas, selection.method = "vst", nfeatures = 2000)
atlas <- ScaleData(atlas)
atlas <- RunPCA(atlas, npcs = 30)

# 4. HARMONY INTEGRATION
print("[2/3] Running Harmony Batch Correction...")
atlas <- RunHarmony(atlas, group.by.vars = "cancer_type")

# 5. UMAP & VISUALIZATION
print("[3/3] Generating Pan-Cancer Manifold...")
atlas <- RunUMAP(atlas, reduction = "harmony", dims = 1:30)
atlas <- FindNeighbors(atlas, reduction = "harmony", dims = 1:30)
atlas <- FindClusters(atlas, resolution = 0.5)

# Plotting
p1 <- DimPlot(atlas, reduction = "umap", group.by = "cancer_type", cols = c("#44AA99", "#CC6677", "#DDCC77")) + 
      ggtitle("Pan-Cancer scRNA-Seq Atlas (Integrated)") + theme_minimal()

# 6. SAVE
ggsave(paste0(out_dir, "PanCancer_Atlas_UMAP.png"), plot = p1, width = 8, height = 6)
saveRDS(atlas, file = "data/curated/Final_Atlas_Seurat.rds")

print("========================================================")
print("SUCCESS: Atlas forged in R. View results in results/ folder.")
print("========================================================")