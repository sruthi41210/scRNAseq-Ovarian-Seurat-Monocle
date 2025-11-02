new_cluster_ids <- c(
  "T cell",                 
  "NK cell",                
  "Macrophage",            
  "Fibroblast / Stromal",   
  "B cell",                 
  "Stromal / CAF",          
  "Mesothelial",            
  "Tumor",                  
  "Secretory Tumor",        
  "Inflammatory-like",      
  "Epithelial",             
  "Dendritic Cell",         
  "Tumor / Hypoxic",        
  "Mast Cell",              
  "Club / Secretory Epi",   
  "Plasma Cell",            
  "Proliferating Tumor",    
  "pDC / Precursor",        
  "Dividing Tumor",         
  "Unknown_19"              # <-- placeholder for the 20th cluster
)

names(new_cluster_ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new_cluster_ids)

DimPlot(seurat_obj, label = TRUE, pt.size = 0.5) +
  ggplot2::ggtitle("UMAP with Annotated Cell Types")

