# 3D Pan-Cancer scRNA-Seq Atlas: Unveiling the Tumor Microenvironment via Advanced Manifold Learning and Evolutionary Trajectories

## Abstract and Study Rationale
The tumor microenvironment (TME) represents a profoundly complex and dynamic ecological niche, distinctly characterized by profound immunosuppression, metabolic competition, and continuous cellular exhaustion. Traditional bulk RNA sequencing methodologies inherently mask this critical intratumoral heterogeneity by providing merely population-level expression averages, thereby obscuring the rare, mechanistic cell states driving therapeutic resistance. 

To comprehensively dissect these intricate multicellular interactions at an individual cell level, this computational pipeline delineates the construction, algorithmic integration, and downstream bioinformatic interpretation of a high-resolution Pan-Cancer Single-Cell RNA-Sequencing (scRNA-Seq) Atlas. 

By systematically integrating normal physiological immune baselines with distinct malignant lineages (specifically mapping circulating hematological profiles against solid tumor infiltrates), we mathematically isolate and analyze the transcriptomic signatures unique to tumor-infiltrating immune cells versus their healthy, circulating counterparts. 

## Architectural Overview
The core repository executes the analytical framework required to normalize, batch-correct, and algorithmically predict immune subsets across multiple cohorts. The pipeline deployed herein encompasses rigorous multi-parameter quality control, non-linear manifold learning, advanced structural batch-effect correction utilizing the Harmony algorithm, unbiased artificial intelligence-driven cell-type annotation via Spearman correlation modeling against primary reference atlases, comprehensive non-parametric differential expression profiling, and pseudotemporal evolutionary trajectory inference.

All scripts execute natively within an R-based ecosystem utilizing Seurat object architectures, mapped precisely against advanced interactive JavaScript visualization engines (Plotly) bounded to Kaggle Jupyter structures.

---

### Phase 1: Robust Quality Control & Data Topography

Prior to engaging in high-dimensional matrix factorization and non-linear dimensional reduction, granular quality control assessment is paramount to prevent mathematically derived artifacts downstream. We evaluate the captured read depth against detected genomic features to isolate dying cells (high mitochondrial leakage) or multiplexed doublets (abnormally high transcript volumes).

While standard violin plots over-smooth multimodal distributions, we utilize high-resolution probability density mappings and a dynamic 3D scatter topology. This reveals the precise capture depth geometry across our respective clinical cohorts. 

![Figure 1a: Quality Control Distribution](figures/fig1a_ridge.png)
*Figure 1a: Ridge density mappings of unique feature detection distributions across healthy baselines, leukemic captures, and solid TME infiltrates.*

![Figure 1b: 3D Topography](figures/fig1b_qc_3d.png)
*Figure 1b: 3D interactive quality control topology matrix correlating UMI depth, total features, and apoptotic mitochondrial leakage.*

### Phase 2: Unharmonized Baseline Geometry vs. Harmonized Manifold

Merging disparate single-cell cohorts inherently introduces profound technical batch effects resulting from disparate sequencing platforms, enzyme kinetics, and independent cell processing environments. In an uncorrected state, variance is dominated by platform origin rather than underlying biological truths.

Prior to integration, the intrinsic 3D Principal Component Analysis (PCA) space maps this extreme divergence. Following identification of mutual variance vectors, we systematically mitigate this utilizing the **Harmony** architectural algorithm. Harmony iteratively identifies mutually nearest neighbors across disparate datasets, computing a smooth, integrated vector embedding. 

![Figure 2: Dysfunctional PCA Space](figures/fig2_pca_3d.png)
*Figure 2: The pre-integration PCA landscape demonstrating cohort-driven artifact clustering.*

![Figure 3: Harmonized Contour Map](figures/fig3_harmony.png)
*Figure 3: Alpha-blended two-dimensional confirmation of biological continuity across integrated coordinates.*

### Phase 3: High-Fidelity 3D Topological Inference

Conveying high-dimensional biological continuity within a severely restricted two-dimensional Euclidean plane enforces extreme structural distortion; mathematically obligating distinct evolutionary lineages to overlap visually.

To mitigate this mathematical compression limit and observe accurate transcriptomic developmental pathways without profound spatial loss, we extract distinct latent vectors (`n.components = 3L`) via UMAP inference mechanisms, mapped directly into multidimensional rotatable visual arrays.

![Figure 4: Computational Immune Topography](figures/fig4_cluster_3d.png)
*Figure 4: Algorithmic identification of unsupervised immune phenotypes across integrated cohorts.*

### Phase 4: Probabilistic Reference-Based Cell Ontology & Molecular Geography

Relying on human visual confirmation of marker genes to annotate thousands of unsupervised clusters represents a significant bottleneck and source of systemic mathematical bias in scRNA-Seq pipelines. Utilizing the `SingleR` paradigm against human primary reference data eradicates this gating error entirely, correlating empirical transcriptomic values globally.

![Figure 5: Spatial AI Diagnostics](figures/fig5_ontology_3d.png)
*Figure 5: High-fidelity ontology projections bypassing human visual bias.*

![Figure 6: Spatial Transcriptomics](figures/fig6_expression_3d.png)
*Figure 6: Direct functional validation mapped across geometric clusters.*

### Phase 5: Deciphering the TME Immunosuppressive Identity & Global Homology

The central inquiry of modern tumor immunology revolves around how the highly hostile, hypoxic, and nutrient-deprived Tumor Microenvironment (TME) actively corrupts and systematically metabolizes otherwise functional immune cells into chronically exhausted or hyper-regulatory phenotypes.

We derive the foundational mathematical identities defining this transformation by computing a high-dimensional non-parametric distinct abundance matrix. Furthermore, to rigorously establish the global homology between these subsets, we compute correlation grids directly mirroring high-resolution viral or species alignments.

![Figure 8: Transcriptomic Similarity Matrix](figures/fig8_heatmap.png)
*Figure 8: Lower-triangular transcriptomic pairwise identity matrix highlighting functional and evolutionary homology.*

### Phase 6: Deep Trajectory Inference and Pseudotemporal Mapping

Cellular biology fundamentally rejects discrete clustering models. Differentiation, physiological exhaustion, and oncological transformation do not manifest instantaneously; rather, they unfold sequentially across smooth temporospatial transcriptomic gradients. By fitting simultaneous principal curves across multidimensional probability spaces, cellular distances are mathematically quantized along the graph structure, deriving explicit developmental "Pseudotime."

![Figure 9: Mathematical Evolution and Vector Trajectory](figures/fig9_trajectory_3d.png)
*Figure 9: Terminal Lineage Trajectories computing pseudotemporal cellular sliding.*

---

## 3. High-Fidelity Research Synthesis

Through the meticulous execution of this highly scaled Single-Cell computational architecture, we have comprehensively mapped the highly discordant transcriptomic geography defining human immune adaptation across major oncological landscapes.

**Core Scientific Achievements:**
1. **Algorithmic Convergence:** Successful execution of structural harmonizing matrices (`Harmony`), enabling deeply disparate primary sequence platforms to align biologically, devoid of compounding technical artifacts.
2. **Computational Ontology:** Implementing multi-dimensional AI inference (`SingleR`) unequivocally eradicated subjective gating interpretation, locking all immune clusters to stringent primary reference profiles via explicit correlation matrices.
3. **Tumor-Induced Reprogramming Maps:** The construction of the Interactive 3D Differential Landscape identified the explicit transcriptomic signatures separating circulating healthy subsets from solid-tumor specific functional collapse.
4. **Lineage Quantization:** Advanced Minimum-Spanning Tree deployment successfully measured and displayed terminal exhaustion vectors, empirically proving the continuous temporal slide triggered by sustained antigen presence within the solid tumor microenvironment.

## Pipeline Execution & Static Verification

To securely render interactive visualization components for local review or markdown conversion without triggering massive JSON payload lockups within web interfaces, a dedicated companion script `scripts/generate_static_figures.R` has been integrated directly into the repository. 

To systematically populate `.png` architectural components for GitHub landing pages or local visual confirmation:
1. Ensure `plotly`, `reticulate`, and the python module `kaleido` are active on your workstation.
2. Run the identical Seurat mathematical architecture natively to print statically, bypassing generic HTML wrapper constraints:
```bash
Rscript scripts/generate_static_figures.R
```

## Large File Governance
**Note:** Explicit `.gitignore` limitations have been applied to this repository. Target dataset volumes defining high-density matrices (`.rds`, `.h5ad`, `.csv`) and their compiled model architectures necessarily scale far beyond the strict 100 MB per-file threshold standard across free-tier GitHub repositories. Repository execution pathways must source physical data independently; placeholder architectural nodes (`/data/`, `/models/`) indicate explicit local directory dependencies containing internal `README.md` descriptors, left functionally open for local data loading prior to execution.
