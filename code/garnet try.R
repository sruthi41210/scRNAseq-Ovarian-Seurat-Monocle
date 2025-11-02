install.packages("remotes")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Now install dependencies
BiocManager::install(c(
  "monocle", 
  "DelayedArray", 
  "DelayedMatrixStats", 
  "org.Hs.eg.db", 
  "org.Mm.eg.db", 
  "AnnotationDbi", 
  "Biobase"
))
classifier <- readRDS("hsLung_20191017.RDS")

counts <- GetAssayData(seurat_obj, layer = "counts")
pd <- new("AnnotatedDataFrame", data = seurat_obj@meta.data)

fd <- data.frame(gene_short_name = rownames(counts))
rownames(fd) <- rownames(counts)
fd <- new("AnnotatedDataFrame", data = fd)

library(monocle)
cds <- newCellDataSet(counts, phenoData = pd, featureData = fd)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)


cds <- classify_cells(cds, classifier, db = org.Hs.eg.db)


pd <- new("AnnotatedDataFrame", data = seurat_obj@meta.data)

fd <- data.frame(gene_short_name = rownames(counts))
rownames(fd) <- rownames(counts)
fd <- new("AnnotatedDataFrame", data = fd)


