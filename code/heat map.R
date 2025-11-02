

# Get scaled gene names properly in Seurat v5
scaled_genes <- rownames(GetAssayData(seurat_obj, layer = "scale.data"))

# Filter genes that are valid
valid_genes <- intersect(top_markers$gene, scaled_genes)

# Plot heat map
DoHeatmap(seurat_obj, features = valid_genes) +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5))
  
