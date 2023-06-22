# Libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(STutility)
set.seed(42)

# I will load the merged final object
# This is the Seurat version of the object on Zenodo
load("spat.data_sct_log311022.RData")

# Let's get the metadata now
meta <- spat.data_sct_log@meta.data

# Subset with cluster 4
c4_obj <- subset(spat.data_sct_log, seurat_clusters == "4")
c4_meta <- c4_obj@meta.data

# I get all the subset files where folded areas
# were manually filtered out using a local shiny app
# Unfortunately this part is not available as it was manually done
files <- list.files("20230428_Analysis_mouse_brain/Data_analysis/Subset_sample/")
  
for (file in files) {
  load(file = paste0("20230428_Analysis_mouse_brain/Data_analysis/Subset_sample/",file))
}

# We will filter the spots that are in the folded regions.
# For this we need to match the rownames. So I fix the subset data

a2_idents <- paste0(
  gsub("1_virgin[0-9]","2",se_A2_OB_virgin20_select@meta.data$labels), "_",
  gsub("_1","",rownames(se_A2_OB_virgin20_select@meta.data))
)

b2_idents <- paste0(
  gsub("_mother_[0-9]","",se_B2_OB_mother_20_select@meta.data$labels), "_",
  gsub("_1","",rownames(se_B2_OB_mother_20_select@meta.data))
)

c2_idents <- paste0(
  gsub("1_virgin[0-9]","2",se_C2_OB_virgin_30_select@meta.data$labels), "_",
  gsub("_1","",rownames(se_C2_OB_virgin_30_select@meta.data))
)

d2_idents <- paste0(
  gsub("_mother[0-9]","",se_D2_OB_mother_30_select@meta.data$labels), "_",
  gsub("_1","",rownames(se_D2_OB_mother_30_select@meta.data))
)

# I will now subset the seurat data with these indices
all_idents <- c(a2_idents, b2_idents, c2_idents, d2_idents)

# I will make some comparisons between original data and
# folds removed data.

# First I am curious about the flat counts of cluster 4
c4_nonfold <- subset(c4_obj, cells = all_idents )

sprintf("Total spots in Cluster 4: %d",length(rownames(c4_obj@meta.data)))
sprintf("Unfolded total spots in Cluster 4: %d",length(rownames(c4_nonfold@meta.data)))
sprintf("Spots coming from folded areas in Cluster 4: %d",
        length(rownames(c4_obj@meta.data))-length(rownames(c4_nonfold@meta.data)))

c4_meta %>% group_by(orig.ident) %>% count()
c4_nonfold@meta.data %>% group_by(orig.ident) %>% count()

# Now checking with the whole data
all_unfold <- subset(spat.data_sct_log, cells = all_idents)

sprintf("Total spots: %d",length(rownames(spat.data_sct_log@meta.data)))
sprintf("Unfolded total spots: %d",length(rownames(all_unfold@meta.data)))
sprintf("Spots coming from folded areas: %d",
        length(rownames(spat.data_sct_log@meta.data))-length(rownames(all_unfold@meta.data)))

spat.data_sct_log@meta.data %>% group_by(orig.ident) %>% count()
all_unfold@meta.data %>% group_by(orig.ident) %>% count()

# The proportions
proportions_table_all <- as.data.frame.matrix(
  prop.table(table(meta$orig.ident, 
                   meta$seurat_clusters), 
             margin = 1)
)

proportions_table_unfold <- as.data.frame.matrix(
  prop.table(table(all_unfold@meta.data$orig.ident, 
                   all_unfold@meta.data$seurat_clusters), 
             margin = 1)
)

# Now I start the sensitivity analysis with the object
# that does not contain any spots in the folded regions

# We do a mild filtering
new_all <- RenameAssays(all_unfold, Spatial = "RNA")
new_all <- subset(new_all, subset = nFeature_RNA > 500 & nFeature_RNA < 8000)
# 4 spots in total are deleted this way
new_all <- SCTransform(new_all)

# Now I want to do clustering again to see if we get a similar cluster
new_all_cluster <- RunPCA(new_all, assay = "SCT", verbose = FALSE)
ElbowPlot(new_all_cluster) # 9 is the elbow of the plot
new_all_cluster <- FindNeighbors(new_all_cluster, reduction = "pca", dims = 1:9)
new_all_cluster <- FindClusters(new_all_cluster, verbose = FALSE, resolution = 0.5)
new_all_cluster <- RunUMAP(new_all_cluster, reduction = "pca", dims = 1:9)

# First checking the new UMAP
pdf(file="sensitivity_umap.pdf",
    width = 16, height = 9)
UMAPPlot(new_all_cluster, label = T, seed = 42)
dev.off()

# Plotting the clusters
pdf(file="mouse_sensitivity_all_clusters.pdf",
    width = 16, height = 9)
SpatialDimPlot(new_all_cluster, crop = F)
dev.off()

# I also plot every cluster separately, side by side
pdf(file="mouse_sensitivity_separate_clusters.pdf",
    width = 16, height = 9)
for (i in 0:(length(unique(new_all_cluster@meta.data$seurat_clusters))-1)){
  a <- SpatialDimPlot(crop = F, new_all_cluster, cells.highlight = CellsByIdentities(object = new_all_cluster, idents = i), facet.highlight = TRUE, ncol = 1, images = "slice1") + scale_fill_manual(values=c("red","NA"))
  b <- SpatialDimPlot(crop= F, new_all_cluster, cells.highlight = CellsByIdentities(object = new_all_cluster, idents = i), facet.highlight = TRUE, ncol = 1, images = "slice1_B2") + scale_fill_manual(values=c("red","NA"))
  c <- SpatialDimPlot(crop = F, new_all_cluster, cells.highlight = CellsByIdentities(object = new_all_cluster, idents = i), facet.highlight = TRUE, ncol = 1, images = "slice1_C2") + scale_fill_manual(values=c("red","NA"))
  d <- SpatialDimPlot(crop= F, new_all_cluster, cells.highlight = CellsByIdentities(object = new_all_cluster, idents = i), facet.highlight = TRUE, ncol = 1, images = "slice1_D2") + scale_fill_manual(values=c("red","NA"))
  print(cowplot::plot_grid(a,b,c,d,nrow = 1, ncol = 4))
}
# crop=F was added later since we wanted to see the real sizes
# of the samples.
dev.off()


# Finding the cluster markers
new_cluster_markers <- FindAllMarkers(new_all_cluster, logfc.threshold = 0.1, assay = "SCT")
new_cluster_markers <- new_cluster_markers[new_cluster_markers$p_val_adj <= 0.05,]
write.csv(new_cluster_markers, "sensitivity_cluster_markers.csv")

# Saving the Seurat object
saveRDS(new_all_cluster, "sensitivity_object.RDS")


# We also check the proportions
as.data.frame.matrix(
  prop.table(table(new_all_cluster@meta.data$orig.ident, 
                   new_all_cluster@meta.data$seurat_clusters), 
             margin = 1)
)

as.data.frame.matrix(
  prop.table(table(spat.data_sct_log@meta.data$orig.ident, 
                   spat.data_sct_log@meta.data$seurat_clusters), 
             margin = 1)
)

# I will now create the heatmap for the sensitivity analysis
# Get top 10 markers for each cluster
top10 <- new_cluster_markers %>% filter(avg_log2FC>0.4) %>% group_by(cluster) %>% top_n(n = 10, wt = p_val_adj)
gene_count <- table(top10$gene)

# Get the genes that only appear once
unique_genes <- names(gene_count)[gene_count == 1]

# Filter the top markers to only include unique markers
unique_top10 <- top10 %>% filter(gene %in% unique_genes)
genes.use <- unique_top10$gene

# Plot the heatmap
pdf(file="sensitivity_heatmap.pdf",
    width = 9, height = 12)
DoHeatmap(object = new_all_cluster, features = genes.use, )
dev.off()
