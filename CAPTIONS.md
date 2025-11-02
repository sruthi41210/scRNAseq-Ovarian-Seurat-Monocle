# Figure Captions — Single-Cell RNA-Seq Analysis

Each figure corresponds to a major step of the Seurat → Monocle → STRING pipeline.

---

## step1_before_filter.png
- Shows: QC metrics before filtering (nFeature_RNA, nCount_RNA, percent.mt)
- Takeaway: Identified outliers and high mitochondrial content cells.

## step1_after_filter.png
- Shows: QC metrics after filtering.
- Takeaway: Clean dataset with improved distribution.

## feature_scatter_before_filter.png
- Shows: Correlation between features and counts.
- Takeaway: Visualized quality variation across cells.

## step2_PCA_Scatter_Plot.png
- Shows: PCA plot after normalization and scaling.
- Takeaway: Separation of cell subpopulations begins to appear.

## step2_PCA_Loadings_Plot.png
- Shows: PCA loadings per component.
- Takeaway: Top genes contributing to variance.

## step2_elbow_plot.png
- Shows: Elbow plot for PC selection.
- Takeaway: Selected top ~15 PCs for clustering.

## step3_Umap.png
- Shows: UMAP projection of all clusters.
- Takeaway: Distinct spatial segregation of populations.

## annotated_umap.png
- Shows: Annotated UMAP with cell type labels.
- Takeaway: Assigned biological identity to each cluster.

## combine_with_cluster_labels_pseudotime.png
- Shows: Monocle pseudotime overlay on UMAP.
- Takeaway: Smooth trajectory progression from immune to proliferative cells.

## pseudotime_plots.png
- Shows: Gene expression trends along pseudotime.
- Takeaway: Early immune → late proliferative patterns.

## 12_tf_factors_pseudotime_analysis.png
- Shows: Top 12 TFs dynamic along pseudotime.
- Takeaway: Key regulators driving trajectory transitions.

## heatmap.png
- Shows: Expression heatmap of top cluster markers.
- Takeaway: Strong gene segregation confirming cluster validity.

## ppi_network.png
- Shows: STRING protein–protein interaction network.
- Takeaway: FOXA2–E2F8–MYBL2 form a central proliferative hub.

## top_5_marker_per_cluster.png
- Shows: Dot plot for top markers per cluster.
- Takeaway: Cluster-specific expression patterns.

## venn_diagram.png
- Shows: Overlap between marker sets.
- Takeaway: Shared vs unique features across clusters.

