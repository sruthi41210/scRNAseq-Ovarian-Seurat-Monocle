# Normalize, find variable features, scale
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

# Run PCA on those variable features
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = 50)


# Visualize PCA loadings
VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
DimPlot(seurat_obj, reduction = "pca")

# Elbow Plot to pick # of PCs
ElbowPlot(seurat_obj, ndims = 50)
