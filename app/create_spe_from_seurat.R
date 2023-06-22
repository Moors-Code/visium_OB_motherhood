library(SpatialExperiment)
library(Seurat)
library(dplyr)
library("BiocFileCache")
library("rtracklayer")
library("lobstr")
library("spatialLIBD")
library("ggplot2")

# load the data from original analysis
load("spat.data_sct_log311022.RData")

# Below function is taken and adapted from:
# https://github.com/drighelli/SpatialExperiment/issues/115

a20b <- subset(x = spat.data_sct_log, subset = orig.ident == "A2_OB_virgin20")
b20b <- subset(x = spat.data_sct_log, subset = orig.ident == "B2_OB_mother20")
c30b <- subset(x = spat.data_sct_log, subset = orig.ident == "C2_OB_virgin30")
d30b <- subset(x = spat.data_sct_log, subset = orig.ident == "D2_OB_mother30")

## Function
# The function is somewhat hard coded to fit our Seurat object.
# Please be aware of possible swaps in coordinates and differences in
# image IDs
seurat_to_spe <- function(seu, sample_id, img_id, assay="Spatial") {
  ## Convert to SCE
  sce <- Seurat::as.SingleCellExperiment(seu, assay=assay)

  ## Extract spatial coordinates
  spatialCoords <- as.matrix(
    seu@images[[img_id]]@coordinates[, c("imagecol", "imagerow")])
  # This colname change is due to hardcoded spatialLIBD things...
  colnames(spatialCoords) <- c("pxl_col_in_fullres", "pxl_row_in_fullres")

  ## Extract and process image data
  img <- SpatialExperiment::SpatialImage(
    x = as.raster(seu@images[[img_id]]@image))

  imgData <- DataFrame(
    sample_id = sample_id,
    image_id = "image",
    data = I(list(img)),
    scaleFactor = seu@images[[img_id]]@scale.factors$lowres)

  # Convert to SpatialExperiment
  spe <- SpatialExperiment(
    assays = assays(sce),
    rowData = rowData(sce),
    colData = colData(sce),
    metadata = metadata(sce),
    reducedDims = reducedDims(sce),
    altExps = altExps(sce),
    sample_id = sample_id,
    spatialCoords = spatialCoords,
    imgData = imgData
  )
  # indicate all spots are on the tissue
  spe$in_tissue <- 1
  spe$sample_id <- sample_id
  # Return Spatial Experiment object
  return(spe)
}


spe1 <- seurat_to_spe(a20b, "A2_OB_virgin20", "slice1", "Spatial")
spe2 <- seurat_to_spe(b20b, "B2_OB_mother20", "slice1_B2", "Spatial")
spe3 <- seurat_to_spe(c30b, "C2_OB_virgin30", "slice1_C2", "Spatial")
spe4 <- seurat_to_spe(d30b, "D2_OB_mother30", "slice1_D2", "Spatial")

spe_ls <- list("A2_OB_virgin20" = spe1,
              "B2_OB_mother20" = spe2,
              "C2_OB_virgin30" = spe3,
              "D2_OB_mother30" = spe4)

spe <- Reduce(cbind, spe_ls)
spe


# Next steps are here:
# http://bioconductor.org/packages/release/data/experiment/vignettes/spatialLIBD/inst/doc/TenX_data_download.html


# These are required by the spatialLIBD
spe$key <- paste0(colnames(spe), "_", spe$sample_id)
stopifnot(!any(duplicated(spe$key)))
spe$sum_umi <- colSums(counts(spe))
spe$sum_gene <- colSums(counts(spe) > 0)

# Adding the gene information
gtf <-
  rtracklayer::import(
    "genes.gtf"
  )
gtf <- gtf[gtf$type == "gene"]

## Set the names to be the gene IDs
names(gtf) <- gtf$gene_name

spe <- spe[-grep("FALSE",rownames(spe) %in% gtf$gene_name),]


## Match the genes
match_genes <- match(rownames(spe), gtf$gene_name)

stopifnot(all(!is.na(match_genes)))

## Keep only some columns from the gtf (you could keep all of them if you want)
mcols(gtf) <-
  mcols(gtf)[, c(
    "source",
    "type",
    "gene_id",
    "gene_version",
    "gene_name",
    "gene_type"
  )]

rowRanges(spe) <- gtf[match_genes]

## Inspect the gene annotation data we added
rowRanges(spe)

## Add information used by spatialLIBD
rowData(spe)$gene_search <- paste0(
  rowData(spe)$gene_name, "; ", rowData(spe)$gene_id
)

## Compute chrM expression and chrM expression ratio
is_mito <- which(seqnames(spe) == "chrM")
spe$expr_chrM <- colSums(counts(spe)[is_mito, , drop = FALSE])
spe$expr_chrM_ratio <- spe$expr_chrM / spe$sum_umi


## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) == 0)

## Number of genes with no counts
length(no_expr)
#> [1] 12136

## Compute the percent of genes with no counts
length(no_expr) / nrow(spe) * 100
#> [1] 37.63684
spe <- spe[-no_expr, , drop = FALSE]

## Remove spots without counts
summary(spe$sum_umi)


if (any(spe$sum_umi == 0)) {
  spots_no_counts <- which(spe$sum_umi == 0)
  ## Number of spots with no counts
  print(length(spots_no_counts))
  ## Percent of spots with no counts
  print(length(spots_no_counts) / ncol(spe) * 100)
  spe <- spe[, -spots_no_counts, drop = FALSE]
}

spe$ManualAnnotation <- "NA"
dim(spe)
lobstr::obj_size(spe) / 1024^2

check_spe(spe)
# Wants gene id as rownames
rownames(spe) <- rowData(spe)$gene_id
check_spe(spe) # It doesn't give any errors now!

# Update: The app still works even if check_spe fails

# Creating some random variable for testing if the library works

set.seed(42)
spe$random_cluster <- sample(1:7, ncol(spe), replace = TRUE)

spatialLIBD::vis_clus(
  spe = spe,
  sampleid = "A2_OB_virgin20",
  clustervar = "random_cluster",
  image_id = "slice1"
)
check_spe(spe)

# Running the default app here.
# From here onwards, cosmetics of the app can be edited by changing
# the www folder that app has.
if (interactive()) {
  run_app(
    spe,
    sce_layer = NULL,
    modeling_results = NULL,
    sig_genes = NULL,
    title = "spatialLIBD: visium OB motherhood",
    spe_discrete_vars = c("SCT_snn_res.0.5", "seurat_clusters",
                          "SCT_snn_res.0.8","ManualAnnotation"),
    spe_continuous_vars = c("percent.ribo", "nCount_SCT","sum_umi",
                            "sum_gene", "expr_chrM", "expr_chrM_ratio"),
    default_cluster = "seurat_clusters"
  )
}

# This object is uploaded to Zenodo
# https://zenodo.org/record/8002261/files/spe.rds?download=1
saveRDS(spe, "spe.RDS")




  ##################################################
 # This part is for the sensitivity analysis part #
##################################################

# Below part is basically copied from above.

sensitivity <- readRDS("sensitivity_object.RDS")
sensitivity <- RenameAssays(object = sensitivity, RNA = "Spatial")

a20b <- subset(x = sensitivity, subset = orig.ident == "A2_OB_virgin20")
b20b <- subset(x = sensitivity, subset = orig.ident == "B2_OB_mother20")
c30b <- subset(x = sensitivity, subset = orig.ident == "C2_OB_virgin30")
d30b <- subset(x = sensitivity, subset = orig.ident == "D2_OB_mother30")

## Function
# The function is somewhat hard coded to fit our Seurat object.
# Please be aware of possible swaps in coordinates and differences in
# image IDs
seurat_to_spe <- function(seu, sample_id, img_id, assay="Spatial") {
  ## Convert to SCE
  sce <- Seurat::as.SingleCellExperiment(seu, assay=assay)
  
  ## Extract spatial coordinates
  spatialCoords <- as.matrix(
    seu@images[[img_id]]@coordinates[, c("imagecol", "imagerow")])
  # This colname change is due to hardcoded spatialLIBD things...
  colnames(spatialCoords) <- c("pxl_col_in_fullres", "pxl_row_in_fullres")
  
  ## Extract and process image data
  img <- SpatialExperiment::SpatialImage(
    x = as.raster(seu@images[[img_id]]@image))
  
  imgData <- DataFrame(
    sample_id = sample_id,
    image_id = "image",
    data = I(list(img)),
    scaleFactor = seu@images[[img_id]]@scale.factors$lowres)
  
  # Convert to SpatialExperiment
  spe <- SpatialExperiment(
    assays = assays(sce),
    rowData = rowData(sce),
    colData = colData(sce),
    metadata = metadata(sce),
    reducedDims = reducedDims(sce),
    altExps = altExps(sce),
    sample_id = sample_id,
    spatialCoords = spatialCoords,
    imgData = imgData
  )
  # indicate all spots are on the tissue
  spe$in_tissue <- 1
  spe$sample_id <- sample_id
  # Return Spatial Experiment object
  return(spe)
}


spe1 <- seurat_to_spe(a20b, "A2_OB_virgin20", "slice1", "Spatial")
spe2 <- seurat_to_spe(b20b, "B2_OB_mother20", "slice1_B2", "Spatial")
spe3 <- seurat_to_spe(c30b, "C2_OB_virgin30", "slice1_C2", "Spatial")
spe4 <- seurat_to_spe(d30b, "D2_OB_mother30", "slice1_D2", "Spatial")

spe_ls <- list("A2_OB_virgin20" = spe1,
               "B2_OB_mother20" = spe2,
               "C2_OB_virgin30" = spe3,
               "D2_OB_mother30" = spe4)

spe <- Reduce(cbind, spe_ls)
spe


# Next steps are here:
# http://bioconductor.org/packages/release/data/experiment/vignettes/spatialLIBD/inst/doc/TenX_data_download.html


# These are required by the spatialLIBD
spe$key <- paste0(colnames(spe), "_", spe$sample_id)
stopifnot(!any(duplicated(spe$key)))
spe$sum_umi <- colSums(counts(spe))
spe$sum_gene <- colSums(counts(spe) > 0)

# Adding the gene information
gtf <-
  rtracklayer::import(
    "genes.gtf"
  )
gtf <- gtf[gtf$type == "gene"]

## Set the names to be the gene IDs
names(gtf) <- gtf$gene_name

spe <- spe[-grep("FALSE",rownames(spe) %in% gtf$gene_name),]


## Match the genes
match_genes <- match(rownames(spe), gtf$gene_name)

stopifnot(all(!is.na(match_genes)))

## Keep only some columns from the gtf (you could keep all of them if you want)
mcols(gtf) <-
  mcols(gtf)[, c(
    "source",
    "type",
    "gene_id",
    "gene_version",
    "gene_name",
    "gene_type"
  )]

rowRanges(spe) <- gtf[match_genes]

## Inspect the gene annotation data we added
rowRanges(spe)

## Add information used by spatialLIBD
rowData(spe)$gene_search <- paste0(
  rowData(spe)$gene_name, "; ", rowData(spe)$gene_id
)

## Compute chrM expression and chrM expression ratio
is_mito <- which(seqnames(spe) == "chrM")
spe$expr_chrM <- colSums(counts(spe)[is_mito, , drop = FALSE])
spe$expr_chrM_ratio <- spe$expr_chrM / spe$sum_umi


## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) == 0)

## Number of genes with no counts
length(no_expr)
#> [1] 12136

## Compute the percent of genes with no counts
length(no_expr) / nrow(spe) * 100
#> [1] 37.63684
spe <- spe[-no_expr, , drop = FALSE]

## Remove spots without counts
summary(spe$sum_umi)


if (any(spe$sum_umi == 0)) {
  spots_no_counts <- which(spe$sum_umi == 0)
  ## Number of spots with no counts
  print(length(spots_no_counts))
  ## Percent of spots with no counts
  print(length(spots_no_counts) / ncol(spe) * 100)
  spe <- spe[, -spots_no_counts, drop = FALSE]
}

spe$ManualAnnotation <- "NA"
dim(spe)
lobstr::obj_size(spe) / 1024^2

check_spe(spe)
# Wants gene id as rownames
rownames(spe) <- rowData(spe)$gene_id
check_spe(spe) # It doesn't give any errors now!


# Creating some random variable for testing if the library works

set.seed(42)
spe$random_cluster <- sample(1:7, ncol(spe), replace = TRUE)

spatialLIBD::vis_clus(
  spe = spe,
  sampleid = "A2_OB_virgin20",
  clustervar = "random_cluster",
  image_id = "slice1"
)
check_spe(spe)

# Running the default app here.
# From here onwards, cosmetics of the app can be edited by changing
# the www folder that app has.
if (interactive()) {
  run_app(
    spe,
    sce_layer = NULL,
    modeling_results = NULL,
    sig_genes = NULL,
    title = "spatialLIBD: visium OB motherhood - sensitivity",
    spe_discrete_vars = c("SCT_snn_res.0.5", "seurat_clusters",
                          "SCT_snn_res.0.8","ManualAnnotation"),
    spe_continuous_vars = c("percent.ribo", "nCount_SCT","sum_umi",
                            "sum_gene", "expr_chrM", "expr_chrM_ratio"),
    default_cluster = "seurat_clusters"
  )
}

# This is the object that is available on Zenodo
# https://zenodo.org/record/8002261/files/sensitivity_spe.RDS?download=1
saveRDS(spe, "sensitivity_spe.RDS")
