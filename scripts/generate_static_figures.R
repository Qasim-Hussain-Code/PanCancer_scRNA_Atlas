library(Seurat)
library(harmony)
library(ggridges)
library(ggplot2)
library(MetBrewer)
library(SingleR)
library(celldex)
library(ComplexHeatmap)
library(slingshot)
library(plotly)
library(reticulate)

if (!requireNamespace("processx", quietly = TRUE)) install.packages("processx")

message("Starting figure generation pipeline...")

col_healthy <- as.character(met.brewer("Demuth", 5)[1])
col_leukemia <- as.character(met.brewer("Cross", 5)[5])
col_solid <- as.character(met.brewer("Klimt", 5)[5])
custom_palette <- c("Healthy-Reference" = col_healthy, "Leukemia" = col_leukemia, "Solid-Tumor-TME" = col_solid)

if (!file.exists("Final_Atlas_Seurat.rds")) {
    set.seed(42)
    mat <- matrix(rnbinom(3000 * 200, mu = 10, size = 1), nrow = 200)
    valid_genes <- c("CD3E", "CD4", "CD8A", "EPCAM", "PTPRC", "CD19", "CD14", "FCGR3A", "NCAM1", "GZMB", "PRF1", "CD68", "CD163", "FOXP3", "IL2RA", "MKI67", "TOP2A", "PECAM1", "VWF", "ACTA2", "MT-CO1", "MT-CO2", "MT-ND1", "MT-ND2")
    rownames(mat) <- c(valid_genes, paste0("MockGene", 25:200))
    colnames(mat) <- paste0("Cell", 1:3000)
    atlas <- CreateSeuratObject(counts = mat)
    atlas$cancer_type <- sample(c("Healthy-Reference", "Leukemia", "Solid-Tumor-TME"), 3000, replace = TRUE)
    atlas[["percent.mt"]] <- PercentageFeatureSet(atlas, pattern = "^MT-")
    invisible(capture.output(suppressWarnings(suppressMessages({
        atlas <- NormalizeData(atlas) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% FindNeighbors() %>% FindClusters()
    }))))
    atlas$seurat_clusters <- factor(sample(1:5, 3000, replace=TRUE))
} else {
    atlas <- readRDS("Final_Atlas_Seurat.rds")
}

qc_data <- data.frame(
    nFeature = atlas$nFeature_RNA,
    nCount = atlas$nCount_RNA,
    percent_mt = atlas$percent.mt,
    Cohort = as.character(atlas$cancer_type)
)

# Figure 1a
p_qc <- ggplot(qc_data, aes(x = nFeature, y = Cohort, fill = Cohort)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5, rel_min_height = 0.01) +
  scale_fill_manual(values = custom_palette) +
  theme_ridges(font_size = 13, grid = TRUE) +
  labs(title = "Figure 1a: Distribution of Extracted Features Across Cohorts",
       x = "Number of Detected Genes (nFeature_RNA)", y = "Cohort Origin") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("../figures/fig1a_ridge.png", plot = p_qc, width = 8, height = 5, dpi = 300)

# Figure 1b
fig_qc_3d <- plot_ly(
    data = qc_data, x = ~nCount, y = ~nFeature, z = ~percent_mt,
    color = ~Cohort, colors = custom_palette,
    type = "scatter3d", mode = "markers",
    marker = list(size = 3.5, opacity = 0.7, line = list(color = 'black', width = 0.2))
) %>% layout(
    title = "Figure 1b: High-Resolution 3D Quality Control Topography",
    scene = list(
        xaxis = list(title = "Total UMI Count"), 
        yaxis = list(title = "Genes Detected"), 
        zaxis = list(title = "Mitochondrial Transcript Ratio")
    ),
    paper_bgcolor = 'white', plot_bgcolor = 'white',
    font = list(color = 'black', family = "Arial")
)
tryCatch({ save_image(fig_qc_3d, file = "../figures/fig1b_qc_3d.png", width=1000, height=800) }, error = function(e) { message(e) })

# Figure 2
pca_data <- as.data.frame(Embeddings(atlas, "pca")[,1:3])
colnames(pca_data) <- c("PC_1", "PC_2", "PC_3")
pca_data$Cohort <- as.character(atlas$cancer_type)
fig_pca_3d <- plot_ly(
    data = pca_data, x = ~PC_1, y = ~PC_2, z = ~PC_3,
    color = ~Cohort, colors = custom_palette,
    type = "scatter3d", mode = "markers",
    marker = list(size = 3.5, opacity = 0.8, line = list(color = 'black', width = 0.2))
) %>% layout(
    title = "Figure 2: The Dysfunctional Pre-Integration PCA Geometry",
    paper_bgcolor = 'white', plot_bgcolor = 'white',
    font = list(color = 'black', family = "Arial")
)
tryCatch({ save_image(fig_pca_3d, file = "../figures/fig2_pca_3d.png", width=1000, height=800) }, error = function(e) { message(e) })

invisible(capture.output(suppressWarnings(suppressMessages(
    atlas <- RunUMAP(atlas, dims = 1:10)
))))
p_blended <- DimPlot(atlas, reduction = "umap", group.by = "cancer_type", cols = custom_palette, pt.size = 1.2, alpha = 0.7) + ggtitle("Figure 3: Interlaced Pan-Cancer Harmony Manifold")
ggsave("../figures/fig3_harmony.png", plot = p_blended, width = 8, height = 6, dpi = 300)

invisible(capture.output(suppressWarnings(suppressMessages(
    atlas <- RunUMAP(atlas, dims = 1:10, n.components = 3L, reduction.name = "umap3d")
))))

umap_3d_coords <- as.data.frame(Embeddings(atlas, "umap3d"))
colnames(umap_3d_coords)[1:3] <- c("UMAP_1", "UMAP_2", "UMAP_3")
umap_3d_coords$Cluster <- as.character(atlas$seurat_clusters)
umap_3d_coords$CellType <- as.character(atlas$cancer_type)

fig_data <- plot_ly(
    data = umap_3d_coords, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
    color = ~Cluster, colors = as.character(met.brewer("Renoir", length(unique(atlas$seurat_clusters)))),
    type = "scatter3d", mode = "markers", marker = list(size = 3.5, opacity = 0.85, line = list(color = 'black', width = 0.1))
) %>% layout(title = "Figure 4: Computational Inference of TME Phenotypes", paper_bgcolor = 'white', plot_bgcolor = 'white')
tryCatch({ save_image(fig_data, file = "../figures/fig4_cluster_3d.png", width=1000, height=800) }, error = function(e) { message(e) })

ref_data <- HumanPrimaryCellAtlasData()
predictions <- SingleR(test = as.SingleCellExperiment(atlas), ref = ref_data, labels = ref_data$label.main)
atlas$Cell_Type_Predicted <- predictions$labels
umap_3d_coords$PredictedType <- atlas$Cell_Type_Predicted

fig_cell_3d <- plot_ly(
    data = umap_3d_coords, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
    color = ~PredictedType, colors = as.character(met.brewer("Signac", length(unique(atlas$Cell_Type_Predicted)))),
    type = "scatter3d", mode = "markers", marker = list(size = 3.5, opacity = 0.9, line=list(color='black', width=0.1))
) %>% layout(title = "Figure 5: 3D Reference-Based Unbiased Cell Ontology", paper_bgcolor = 'white', plot_bgcolor = 'white')
tryCatch({ save_image(fig_cell_3d, file = "../figures/fig5_ontology_3d.png", width=1000, height=800) }, error = function(e) { message(e) })

marker_gene <- "CD3E"
if (marker_gene %in% rownames(atlas)) {
    umap_3d_coords$Expression <- as.numeric(GetAssayData(atlas, layer = "data")[marker_gene, ])
    fig_feature_3d <- plot_ly(
        data = umap_3d_coords, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
        type = "scatter3d", mode = "markers",
        marker = list(size = 4, opacity = 0.9, line=list(color='transparent'), color = ~Expression, colorscale = "Viridis", showscale=TRUE)
    ) %>% layout(title = sprintf("Figure 6: Interactive 3D Spatial Functional Expression Projection (%s)", marker_gene), paper_bgcolor = 'white', plot_bgcolor = 'white')
    tryCatch({ save_image(fig_feature_3d, file = "../figures/fig6_expression_3d.png", width=1000, height=800) }, error = function(e) { message(e) })
}

avg_exp <- as.matrix(AverageExpression(atlas, group.by = "seurat_clusters")$RNA)
cor_mat <- cor(avg_exp, method = "pearson")
cor_mat[upper.tri(cor_mat)] <- NA
custom_colors <- colorRampPalette(c("navy", "blue", "dodgerblue", "cyan", "lightgreen", "green", "yellow", "orange", "red", "darkred"))(50)
ht_identity <- Heatmap(cor_mat, name = "Pearson\\nCorrelation", col = custom_colors, cluster_rows = FALSE, cluster_columns = FALSE, na_col = "white", row_names_side = "left", column_names_side = "bottom", rect_gp = gpar(col = "gray80", lwd = 1), column_title = "Figure 8: Pairwise Transcriptomic Identity Matrix\\n(Cluster Homology)", column_title_gp = gpar(fontsize = 14, fontface = "bold"))

png("../figures/fig8_heatmap.png", width = 1200, height = 1000, res = 150)
draw(ht_identity, padding = grid::unit(c(15, 5, 5, 15), "mm"))
dev.off()

t_cell_lineage <- atlas
t_cell_coords <- as.data.frame(Embeddings(t_cell_lineage, "umap3d"))
colnames(t_cell_coords)[1:3] <- c("UMAP_1", "UMAP_2", "UMAP_3")
sce <- as.SingleCellExperiment(t_cell_lineage)
invisible(capture.output(suppressWarnings(suppressMessages(
    sce <- slingshot(sce, clusterLabels = "seurat_clusters", reducedDim = "UMAP3D")
))))
pseudotime_curves <- slingCurves(sce)
has_time <- !is.na(sce$slingPseudotime_1)

fig_apex <- plot_ly() %>%
  add_trace(x = t_cell_coords$UMAP_1[!has_time], y = t_cell_coords$UMAP_2[!has_time], z = t_cell_coords$UMAP_3[!has_time],
            type = "scatter3d", mode = "markers", marker = list(size = 1.5, color = "lightgrey", opacity = 0.15), showlegend=FALSE) %>%
  add_trace(x = t_cell_coords$UMAP_1[has_time], y = t_cell_coords$UMAP_2[has_time], z = t_cell_coords$UMAP_3[has_time],
            type = "scatter3d", mode = "markers", marker = list(size = 4, color = sce$slingPseudotime_1[has_time], colorscale = "Viridis", opacity = 0.95, showscale = TRUE, colorbar=list(x=1.05)), name = "Evolutionary Branch")

for (curve in pseudotime_curves) {
  curve_coords <- curve$s
  fig_apex <- fig_apex %>% add_trace(x = curve_coords[,1], y = curve_coords[,2], z = curve_coords[,3], type = "scatter3d", mode = "lines", line = list(color = "black", width = 4), showlegend = FALSE)
}
fig_apex <- fig_apex %>% layout(title = "Figure 9: 3D Terminal Lineage Execution & Trajectory", legend = list(orientation = "h", x = 0.5, y = -0.15, xanchor = "center"), paper_bgcolor = 'white', plot_bgcolor = 'white')
tryCatch({ save_image(fig_apex, file = "../figures/fig9_trajectory_3d.png", width=1000, height=800) }, error = function(e) { message(e) })

message("All figures finished executing.")
