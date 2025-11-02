# Load required libraries
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

# Set working directory
setwd("C:/Users/sruth/Data_Zhang2022_Ovarian")

# Load files
count_matrix <- readMM("Exp_data_UMIcounts.mtx")
genes <- read.table("Genes.txt", header = FALSE, stringsAsFactors = FALSE)
cells <- read.csv("Cells.csv", header = TRUE)

# Assign row and column names to count matrix
rownames(count_matrix) <- genes$V1
colnames(count_matrix) <- cells$cell_name

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = count_matrix, project = "Zhang_Ovarian", min.cells = 3)

# Add metadata to Seurat object
metadata_cols <- colnames(cells)[colnames(cells) != "cell_name"]
for (col in metadata_cols) {
  seurat_obj[[col]] <- cells[[col]]
}

# If percent.mt is already present in your metadata, this is not needed:
if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
}

# QC Plots (before filtering)
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Apply QC filtering
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# QC Plots (after filtering)
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# View summary
seurat_obj

