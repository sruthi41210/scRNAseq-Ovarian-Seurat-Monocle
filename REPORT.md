# REPORT — Single-Cell RNA-Seq (Seurat → Monocle)

## Methods (Detailed)
- **Preprocessing (Seurat):** Normalize → Identify HVGs (~2k) → Scale → PCA (up to 50 PCs) → Elbow plot.
- **Embedding & Clustering:** UMAP from top PCs → FindNeighbors/FindClusters (res≈0.5, tune as needed).
- **Markers & Annotation:** FindAllMarkers per cluster; annotate using canonical markers.
- **Pseudotime (Monocle 3):** Convert to CDS → preprocess (dims=50) → UMAP → learn_graph → order_cells (root on immune-like cluster).
- **Enrichment:** GO/KEGG on cluster markers/pseudotime modules; STRING PPI for top gene sets.

## Results (Narrative)
- Early pseudotime enriched for immune/cytokine TFs (e.g., BATF, XCL1); later stages show proliferative/cell-cycle TFs (e.g., MYBL2, E2F8, FOXA2).
- Enrichment suggests a shift from immune-responsive states to proliferative/hypoxia programs along trajectory.

## Final Hypothesis
Immune-primed cell states dominate the early trajectory, transitioning toward proliferation-associated programs later, consistent with a progression from response to stimulus → growth and division. This aligns with GO/KEGG terms (cytokine signaling early; cell cycle/hypoxia later) observed in the figures.

> Replace/augment this text with your exact phrasing from the handwritten submission if needed.
