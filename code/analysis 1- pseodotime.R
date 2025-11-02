library(Seurat)
library(monocle3)
library(SeuratWrappers)

cds <- as.cell_data_set(seurat_obj)


cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, reduction_method = "UMAP")

cds <- cluster_cells(cds)

cds <- learn_graph(cds)
cds <- order_cells(cds)

plot_genes_in_pseudotime(cds[c("BATF", "E2F8", "MYBL2", "XCL1", "FOXA2"), ])


rowData(cds)$gene_short_name <- rownames(cds)


plot_genes_in_pseudotime(cds[c(
  "BATF", "XCL1", "OSR1", "PLA2G10", "EME1", "PAX2",
  "FOXA2", "LTF", "HLA-DQB2", "E2F8", "MYBL2", "HIST1H1B"
), ],
color_cells_by = "pseudotime")
