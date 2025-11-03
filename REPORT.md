# Ovarian Cancer Single-Cell RNA-Seq Analysis — Zhang et al. (2022)

## Overview

This report presents a comprehensive single-cell RNA-sequencing (scRNA-seq) analysis of ovarian cancer using the publicly available dataset **GSE244574** from Zhang et al. (2022). The goal of this study was to investigate the transcriptional heterogeneity within the tumor microenvironment, focusing on transcription factor (TF)-driven transitions from immune-active to proliferative tumor states. The analysis integrates Seurat-based clustering, transcription factor identification through Venny and STRING-db, functional enrichment via ShinyGO, and pseudotime trajectory modeling with Monocle3.

---

## 1. Dataset and Tools

- **Dataset:** GSE244574 (Zhang et al., *Science Advances*, 2022)
- **Objective:** To identify transcriptional regulators controlling immune-to-tumor transitions in ovarian cancer.
- **Tools Used:** R (Seurat, Monocle3), STRING-db, ShinyGO, ggplot2, Venny
- **Data Type:** Single-cell RNA-seq (scRNA-seq) expression matrix from ovarian tumor tissue samples.

---

## 2. Clustering and Annotation (Seurat)

Dimensionality reduction was performed using PCA and UMAP, followed by graph-based clustering. The top marker genes per cluster were identified using `FindAllMarkers()` and used for biological annotation.

**Cell Types Identified:**
- **Immune:** T cells (CD3G, GIMAP7), NK cells (XCL1, GNLY), macrophages (C1QA, TREM2), dendritic cells (FLT3)
- **Tumor:** Tumor epithelial (PAX2, FOXA2), proliferative/hypoxic tumor (MKI67, AURKB, E2F8)
- **Stromal:** CAFs (MMP11, POSTN), fibroblasts (OSR1), mesothelial cells (CALML5)
- **Others:** B cells, plasma cells, mast cells, club/secretory epithelial

**Inference:** The UMAP projection demonstrated clear separation between immune and tumor compartments, confirming robust clustering and distinct subpopulation identity.

---

## 3. Top Marker Gene Extraction and Visualization

The top 10 differentially expressed genes (DEGs) per cluster were extracted based on average log2 fold-change values. Heatmaps and dot plots were used to visualize cluster-specific gene expression patterns.

**Key Observations:**
- Immune clusters expressed *BATF*, *XCL1*, and *CD3G*.
- Tumor clusters were characterized by *FOXA2*, *MKI67*, and *MYBL2*.
- Stromal clusters expressed *OSR1* and *POSTN*, indicating matrix remodeling roles.

These patterns confirmed the biological validity of the cluster assignments and aligned with the Zhang et al. (2022) annotations.

---

## 4. Transcription Factor Filtering and Network Analysis

To identify key regulatory factors, a curated human TF list was cross-referenced with the DEGs using **Venny**.

**12 transcription factors (TFs)** were identified:
```
BATF, XCL1, OSR1, PLA2G10, EME1, PAX2, FOXA2, LTF, HLA-DQB2, E2F8, MYBL2, HIST1H1B
```

**Functional Grouping:**
- **Immune-related TFs:** BATF, XCL1, HLA-DQB2
- **Proliferative TFs:** E2F8, MYBL2, FOXA2, HIST1H1B
- **Stromal/Other TFs:** OSR1, PAX2, LTF, PLA2G10, EME1

STRING-db network analysis revealed two major hubs:
1. **Immune Hub:** BATF and HLA-DQB2 connected to immune activation and T-cell signaling.
2. **Proliferation Hub:** E2F8 and MYBL2 central to cell-cycle regulation and replication stress.

An additional link was observed between **PLA2G10** (inflammatory signaling) and the proliferative module, suggesting crosstalk between inflammation and tumor proliferation.

---

## 5. Gene Ontology and Pathway Enrichment (ShinyGO)

Functional enrichment analysis was performed using **ShinyGO** on the 12 TFs and their co-expressed genes.

**Top GO Biological Processes:**
- Leukocyte migration and immune activation
- Cytokine-mediated signaling pathway
- Regulation of cell proliferation and apoptosis

**Top KEGG Pathways:**
- IL-17 signaling pathway
- Cytokine-cytokine receptor interaction
- PD-L1 expression and PD-1 checkpoint pathway
- Graft-versus-host disease (GvHD)

**Interpretation:** The enrichment of immune activation and checkpoint pathways suggests an early hyperactive immune state followed by an immune exhaustion phase. The GvHD-like signature reflects immune overstimulation commonly seen in chronic tumor inflammation.

---

## 6. Pseudotime Trajectory Analysis (Monocle3)

Monocle3 was used to order cells along a pseudotime trajectory beginning from T cell clusters. Expression patterns of the 12 transcription factors were plotted across pseudotime to model temporal transcriptional dynamics.

**Key Results:**
- **Early pseudotime TFs:** BATF, XCL1 — enriched in T and NK cells, promoting cytokine signaling and immune activation.
- **Intermediate TFs:** FOXA2 — associated with epithelial transition and potential epithelial-mesenchymal transition (EMT).
- **Late pseudotime TFs:** E2F8, MYBL2 — enriched in proliferative tumor clusters, driving cell-cycle progression and hypoxia adaptation.

**Inference:** The trajectory highlights a clear immune → proliferative shift. Early immune-regulatory TFs diminish as proliferative and survival TFs dominate, suggesting transcriptional reprogramming underlying tumor progression.

---

## 7. Comparison with Literature

The Zhang et al. (2022) study characterized immune and stromal diversity but did not explicitly quantify transcriptional dynamics over pseudotime. The present analysis extends their findings by uncovering **temporal regulation of transcription factors** that align with functional transitions from immune engagement to tumor proliferation.

**Key Literature Correlations:**
- *BATF* regulates CD8+ T cell activation (PMID: 31036946)
- *MYBL2* correlates with poor prognosis in ovarian cancer (PMID: 29211714)
- *FOXA2* controls epithelial plasticity and tumor progression (PMID: 25979562)

This study integrates and expands these findings into a coherent pseudotemporal model of transcriptional reprogramming.

---

## 8. Hypothesis and Conclusion

> Ovarian tumors exhibit a pseudotemporal transition from immune-active to proliferative and hypoxic states, regulated by a stage-specific transcription factor network. Early TFs (BATF, XCL1) mediate immune activation, while later TFs (FOXA2, E2F8, MYBL2) promote proliferation and immune suppression.

This model supports a sequential transcriptional shift underlying tumor evolution. The immune overactivation observed early in pseudotime may drive chronic inflammation, eventually leading to immune exhaustion and tumor immune escape.

**Therapeutic Implications:**
- Targeting early immune TFs could enhance anti-tumor responses.
- Inhibiting late proliferative TFs could prevent tumor expansion and resistance.

---

## 9. Future Directions

- Validate TF expression dynamics experimentally using CRISPR or RNA interference in ovarian tumor lines.
- Perform SCENIC or regulon-based TF activity analysis to refine transcriptional regulation networks.
- Integrate spatial transcriptomics or proteomic data to confirm localization of pseudotime states.
- Extend the analysis to other cancers (e.g., breast, endometrial) to evaluate the generality of the immune-to-proliferative transition model.

---

## Reference

Zhang, Q., et al. (2022). *Landscape and dynamics of single immune cells in ovarian cancer.* Science Advances, 8(9), eabm1831.  
DOI: [10.1126/sciadv.abm1831](https://www.science.org/doi/10.1126/sciadv.abm1831)
