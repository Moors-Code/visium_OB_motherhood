#BiocManager version 3.16
#SpatialLIBD version 1.13.6
#R 4.2.2

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("spatialLIBD")

#set the working directory to the folder where "sensitivity_spe.RDS" is stored.
setwd("")

# We load the libraries
library("spatialLIBD")
library("markdown") ## due to a bug it needs to be called

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## Load the data (all paths are relative to this script's location)
## spe.rds is a SpatialExperiment object adjusted for spatialLIBD library.
# sensitivity_spe.rds is the sensitivity analysis part.
## Genes with no expression levels were filtered out to shrink the size
## of the object
# spe_wrapper <- readRDS("spe.rds")
spe_wrapper <- readRDS("sensitivity_spe.RDS")
vars <- colnames(colData(spe_wrapper))

## Run the app. This will start a separate window
# www folder for the docs_path can be either the default one from spatialLIBD
# or the one from GitHub can be downloaded and decompressed
run_app(
  spe_wrapper,
  sce_layer = NULL,
  auto_crop_default = FALSE,
  modeling_results = NULL,
  sig_genes = NULL,
  title = "spatialLIBD: Mouse brain results",
  spe_discrete_vars = c("SCT_snn_res.0.5", "seurat_clusters",
                        "SCT_snn_res.0.8","ManualAnnotation"),
  spe_continuous_vars = c("nCount_SCT", "percent.ribo","sum_umi",
                          "sum_gene", "expr_chrM", "expr_chrM_ratio"),
  default_cluster = "seurat_clusters",
  docs_path = "www"
)
