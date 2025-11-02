# Single-Cell RNA-Seq Analysis of Ovarian Dataset (Seurat → Monocle → STRING)

## Introduction

This study investigates transcriptional heterogeneity in an ovarian single-cell RNA-seq dataset.  
The goal is to characterize **immune-to-proliferative transcriptional transitions** using Seurat (for clustering) and Monocle (for pseudotime analysis).

---

## Methodology

1. **Preprocessing & Quality Control (QC)**
   - Removed low-quality cells based on mitochondrial content and total detected features.
   - Ensured consistent nFeature_RNA and nCount_RNA distributions after filtering.

2. **Normalization, HVG Selection, PCA**
   - Normalized data and identified ~2,000 HVGs.
   - Performed PCA to capture major axes of variance.

3. **Elbow Plot & PCA Loadings**
   - Determined top 15 PCs for clustering.
   - Identified key variance-driving genes through PCA loadings.

4. **Clustering & UMAP Visualization**
   - Constructed nearest-neighbor graphs using first 15 PCs.
   - Used UMAP to visualize distinct clusters and identify biologically relevant groups.

5. **Marker Gene Identification**
   - Used `FindAllMarkers()` to find cluster-specific genes.
   - Visualized expression with dot plots and heatmaps.

6. **Pseudotime Trajectory Analysis (Monocle 3)**
   - Converted Seurat object into a Monocle `CellDataSet`.
   - Ordered cells along pseudotime to infer progression patterns.
   - Visualized 12 key transcription factors along pseudotime.

7. **Functional Enrichment & Network Analysis**
   - Conducted GO/KEGG enrichment of key gene sets.
   - Used STRING to identify hub interactions among TFs.

---

## Graph Explanations

- **QC Plots:** Visualize improvements in mitochondrial % and gene distribution after filtering.
- **PCA/Elbow:** Highlights main components driving cell variance and supports clustering depth.
- **UMAP:** Displays distinct immune, stromal, and epithelial clusters.
- **Pseudotime:** Orders cells along differentiation; early immune → late proliferative transitions observed.
- **TF Trends:** 12 TFs show phase-wise regulation consistent with pseudotime trajectory.
- **Heatmap & Marker Plots:** Validate the identity and distinct expression modules of each cluster.
- **STRING Network:** FOXA2–E2F8–MYBL2 hub reinforces proliferative cluster role.
- **Venn Diagram:** Shows overlap and exclusivity among gene sets.

---

## Results Summary

- **Immune Activation (Early pseudotime):** High *BATF*, *XCL1*, *CCL5* expression — cytokine and inflammation pathways.  
- **Proliferation (Late pseudotime):** High *E2F8*, *MYBL2*, *FOXA2* expression — cell-cycle and hypoxia pathways.  
- **Network Integration:** FOXA2–E2F8–MYBL2 hub connects transcriptional drivers of proliferation.  
- **GO/KEGG Enrichment:** Cytokine and immune terms early; cell-cycle, DNA replication, and hypoxia terms later.

---

## Final Hypothesis

> The ovarian single-cell transcriptome follows a **continuous immune-to-proliferative trajectory**.  
> Cells begin with cytokine and immune response activation, later transitioning into proliferative, hypoxia-driven states.  
> This progression mirrors tumor adaptation and potential therapy resistance, driven by transcription factors *FOXA2*, *E2F8*, and *MYBL2*.

---

## Future Work

- Integrate ligand–receptor communication networks (e.g., *CellChat*).  
- Cross-validate trajectory modules across multiple ovarian datasets.  
- Extend pseudotime modeling to drug-response predictions.
