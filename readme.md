# Single-Cell RNA-Seq Analysis of Ovarian Dataset (Seurat → Monocle → STRING)

**Author:** [@sruthi41210](https://github.com/sruthi41210)  
**Tools:** R, Seurat, Monocle 3, STRING, ggplot2  
**Focus:** Transcriptional programs, cell-type annotation, and pseudotime trajectory in ovarian single-cell RNA-seq data.

---

##  Project Overview

This project analyzes an ovarian single-cell RNA-seq dataset to explore **transcriptional transitions** from immune to proliferative cell states using Seurat and Monocle.

**Pipeline Summary**
1. **Quality Control (QC)** – Filtering low-quality cells by mitochondrial %, feature count, and outliers.  
2. **Normalization & PCA** – Feature scaling, HVG selection, PCA, and elbow analysis.  
3. **Clustering & UMAP** – Dimensional reduction and visualization of major cell populations.  
4. **Marker Identification** – Differential gene expression to identify top cluster markers.  
5. **Annotation & Integration** – Annotating cell types and aligning with known immune and proliferative markers.  
6. **Pseudotime Trajectory (Monocle 3)** – Ordering cells along differentiation trajectories.  
7. **Enrichment & PPI (STRING)** – GO/KEGG enrichment and network-level validation.

---

##  Key Results

- **Immune → Proliferative Transition:**  
  Early pseudotime regions enriched for immune-related TFs (e.g., *BATF*, *XCL1*, *CCL5*), later regions showed proliferative TFs (*E2F8*, *FOXA2*, *MYBL2*).  

- **Enrichment Patterns:**  
  - Early trajectory → cytokine signaling, immune activation pathways.  
  - Late trajectory → cell-cycle, hypoxia, and DNA repair pathways.  

- **STRING Network Insight:**  
  *FOXA2–E2F8–MYBL2* network acts as a central proliferative hub.

---

##  Hypothesis

> Ovarian single-cell transcriptomes show a **sequential transcriptional reprogramming** from immune-primed to proliferative states.  
> This shift likely represents tumor progression from immune activation toward sustained proliferation, consistent with the enriched cytokine-to-cell-cycle axis.

---

##  Visual Outputs

### QC & Filtering
![QC Before](assets/images/step1_before_filter.png)
![QC After](assets/images/step1-after_filter.png)
![Feature Scatter](assets/images/feature_scatter_before_filtering.png)

### PCA & Dimensionality
![PCA Scatter](assets/images/step%202%20PCA%20Scatter%20Plot.png)
![PCA Loadings](assets/images/step2_PCA_Loadings_Plot.png)
![Elbow Plot](assets/images/step%202%20elbow%20plot.png)

### Clustering & Annotation
![UMAP](assets/images/step%203%20Umap.png)
![Annotated UMAP](assets/images/annotated%20umap.png)

### Pseudotime
![Pseudotime Overlay](assets\images\combine%20with%20cluster%20labels%20pseodotime.png)
![Pseudotime Trends](assets/images/pseodotime%20plots.png)
![TF Trends (12)](assets/images/12%20tf%20factors%20pseodotime%20analysis.png)

### Marker/Heatmap/Networks
![Top 5 Markers](assets/images/top%205%20marker%20per%20clustor.png)
![Heatmap](assets/images/heatmap.png)
![PPI Network](assets/images/ppi%20network.png)
![Venn](assets/images/venn%20diagram.png)

> Full captions and interpretations are in [`CAPTIONS.md`](CAPTIONS.md).

---

##  Repository Structure
```
scRNAseq-Ovarian-Seurat-Monocle/
│
├─ assets/
│  └─ images/          # Analysis figures (QC, PCA, UMAP, pseudotime, heatmap, STRING)
│
├─ code/
│  ├─ step1_qc.R
│  ├─ step2_pca.R
│  ├─ step3_umap.R
│  ├─ step4_markers.R
│  ├─ heatmap.R
│  ├─ rename_and_annotate.R
│  ├─ analysis1_pseudotime.R
│  ├─ garnet_try.R
│  └─ scrna_workflow.R
│
├─ REPORT.md           # Full workflow + hypothesis
├─ CAPTIONS.md         # Captions & interpretations for each figure
├─ .gitignore
└─ README.md
```

---

##  Related Files

-  [`REPORT.md`](REPORT.md) – Extended write-up of the experiment and findings.  
-  [`CAPTIONS.md`](CAPTIONS.md) – Per-figure explanations and insights.  
-  R scripts under [`/code`](code/) – End-to-end pipeline implementation.

---

## one-line summary

> **Single-Cell RNA-Seq Analysis (Seurat → Monocle → STRING):**  
> Identified immune-to-proliferative transcriptional transitions in ovarian dataset; performed QC, clustering, pseudotime analysis, and enrichment validation.

---

##  Citation

If you reference this project in a portfolio or academic submission, please cite.

```
Sruthi S. (2025). Single-Cell RNA-Seq Analysis of Ovarian Dataset (Seurat → Monocle → STRING).  
GitHub: https://github.com/sruthi41210/scRNAseq-Ovarian-Seurat-Monocle
```
