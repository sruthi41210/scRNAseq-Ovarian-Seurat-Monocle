# Clustering using top 20 PCs
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Run UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# Plot UMAP with cluster labels
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5)
