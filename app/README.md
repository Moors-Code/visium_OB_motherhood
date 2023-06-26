# Visualisation App



`create_spe_from_seurat.R`


## Code description

1. `create_spe_from_seurat.R`
* Uses [shiny](https://shiny.rstudio.com/) web application,
* Uses Bioconductor package at
    [bioconductor.org/packages/spatialLIBD](http://bioconductor.org/packages/spatialLIBD)
    (or from [here](http://research.libd.org/spatialLIBD/))
* Uses `spe.rds` form accessable via [Zenodo][1] as input.
 
[1]:https://zenodo.org/record/8002261

## SessionInfo

<pre>R version 4.1.0 (2021-05-18)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.2 LTS

Matrix products: default
BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=C              LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices datasets  utils     methods   base     

loaded via a namespace (and not attached):
  [1] rappdirs_0.3.3                spatialLIBD_1.11.1            rtracklayer_1.54.0            scattermore_0.8               SpatialExperiment_1.4.0      
  [6] R.methodsS3_1.8.2             SeuratObject_4.1.3            tidyr_1.2.1                   ggplot2_3.4.0                 bit64_4.0.5                  
 [11] knitr_1.41                    irlba_2.3.5.1                 DelayedArray_0.20.0           R.utils_2.12.2                data.table_1.14.6            
 [16] KEGGREST_1.34.0               RCurl_1.98-1.9                doParallel_1.0.17             generics_0.1.3                BiocGenerics_0.40.0          
 [21] ScaledMatrix_1.2.0            cowplot_1.1.1                 usethis_2.1.6                 RSQLite_2.2.19                RANN_2.6.1                   
 [26] future_1.29.0                 bit_4.0.5                     spatstat.data_3.0-0           xml2_1.3.3                    httpuv_1.6.6                 
 [31] SummarizedExperiment_1.24.0   assertthat_0.2.1              viridis_0.6.2                 xfun_0.35                     jquerylib_0.1.4              
 [36] promises_1.2.0.1              fansi_1.0.3                   restfulr_0.0.15               dbplyr_2.2.1                  igraph_1.3.5                 
 [41] DBI_1.1.3                     htmlwidgets_1.5.4             spatstat.geom_3.0-3           stats4_4.1.0                  paletteer_1.5.0              
 [46] purrr_0.3.5                   ellipsis_0.3.2                benchmarkmeData_1.0.4         dplyr_1.0.10                  deldir_1.0-6                 
 [51] sparseMatrixStats_1.6.0       MatrixGenerics_1.6.0          vctrs_0.5.1                   SingleCellExperiment_1.16.0   Biobase_2.54.0               
 [56] ROCR_1.0-11                   abind_1.4-5                   cachem_1.0.6                  progressr_0.11.0              sctransform_0.3.5            
 [61] GenomicAlignments_1.30.0      goftest_1.2-3                 cluster_2.1.2                 ExperimentHub_2.2.1           dotCall64_1.0-2              
 [66] lazyeval_0.2.2                crayon_1.5.2                  spatstat.explore_3.0-5        edgeR_3.36.0                  pkgconfig_2.0.3              
 [71] GenomeInfoDb_1.30.1           vipor_0.4.5                   nlme_3.1-152                  pkgload_1.3.2                 rlang_1.0.6                  
 [76] globals_0.16.2                lifecycle_1.0.3               miniUI_0.1.1.1                filelock_1.0.2                BiocFileCache_2.2.1          
 [81] rsvd_1.0.5                    AnnotationHub_3.2.2           rprojroot_2.0.3               polyclip_1.10-4               matrixStats_0.63.0           
 [86] lmtest_0.9-40                 Matrix_1.5-3                  Rhdf5lib_1.16.0               zoo_1.8-11                    beeswarm_0.4.0               
 [91] ggridges_0.5.4                png_0.1-7                     viridisLite_0.4.1             rjson_0.2.21                  bitops_1.0-7                 
 [96] R.oo_1.25.0                   KernSmooth_2.23-20            spam_2.9-1                    rhdf5filters_1.6.0            Biostrings_2.62.0            
[101] blob_1.2.3                    DelayedMatrixStats_1.16.0     stringr_1.4.1                 parallelly_1.32.1             spatstat.random_3.0-1        
[106] S4Vectors_0.32.4              beachmat_2.10.0               scales_1.2.1                  memoise_2.0.1                 magrittr_2.0.3               
[111] plyr_1.8.8                    ica_1.0-3                     zlibbioc_1.40.0               compiler_4.1.0                dqrng_0.3.0                  
[116] BiocIO_1.4.0                  RColorBrewer_1.1-3            fitdistrplus_1.1-8            Rsamtools_2.10.0              cli_3.4.1                    
[121] XVector_0.34.0                listenv_0.8.0                 patchwork_1.1.2               pbapply_1.6-0                 MASS_7.3-54                  
[126] tidyselect_1.2.0              stringi_1.7.8                 yaml_2.3.6                    BiocSingular_1.10.0           locfit_1.5-9.6               
[131] ggrepel_0.9.2                 grid_4.1.0                    sass_0.4.4                    lobstr_1.1.2                  tools_4.1.0                  
[136] future.apply_1.10.0           parallel_4.1.0                rstudioapi_0.14               foreach_1.5.2                 gridExtra_2.3                
[141] Rtsne_0.16                    DropletUtils_1.14.2           digest_0.6.30                 BiocManager_1.30.19           shiny_1.7.3                  
[146] Rcpp_1.0.9                    GenomicRanges_1.46.1          scuttle_1.4.0                 BiocVersion_3.14.0            later_1.3.0                  
[151] shinyWidgets_0.7.5            RcppAnnoy_0.0.20              httr_1.4.4                    AnnotationDbi_1.56.2          colorspace_2.0-3             
[156] XML_3.99-0.12                 fs_1.5.2                      tensor_1.5                    reticulate_1.26               IRanges_2.28.0               
[161] splines_4.1.0                 fields_14.1                   statmod_1.4.37                uwot_0.1.14                   rematch2_2.1.2               
[166] spatstat.utils_3.0-1          scater_1.22.0                 sp_1.5-1                      renv_0.16.0                   sessioninfo_1.2.2            
[171] plotly_4.10.1                 xtable_1.8-4                  jsonlite_1.8.3                benchmarkme_1.0.8             R6_2.5.1                     
[176] pillar_1.8.1                  htmltools_0.5.3               mime_0.12                     glue_1.6.2                    fastmap_1.1.0                
[181] DT_0.26                       BiocParallel_1.28.3           BiocNeighbors_1.12.0          interactiveDisplayBase_1.32.0 codetools_0.2-18             
[186] maps_3.4.1                    utf8_1.2.2                    bslib_0.4.1                   lattice_0.20-44               spatstat.sparse_3.0-0        
[191] tibble_3.1.8                  ggbeeswarm_0.6.0              curl_4.3.3                    leiden_0.4.3                  attempt_0.3.1                
[196] config_0.3.1                  magick_2.7.3                  golem_0.3.5                   survival_3.2-11               limma_3.50.3                 
[201] roxygen2_7.2.2                desc_1.4.2                    munsell_0.5.0                 rhdf5_2.38.1                  GenomeInfoDbData_1.2.7       
[206] iterators_1.0.14              HDF5Array_1.22.1              reshape2_1.4.4                gtable_0.3.1                  Seurat_4.3.0                 
</pre>





