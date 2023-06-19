# PREGNANCY-RESPONSIVE POOLS OF ADULT NEURAL STEM CELLS FOR TRANSIENT NEUROGENESIS IN MOTHERS


This repository contains the processing pipeline for 10x Visium from [xxxxxxx] [1] structured into 3 modules Rmarkdown scripts.

[1]:https to paper
<pre>`1_Normalisation_subset_merge.Rmd`
`2_DGE_Table_Paper.Rmd`
`120123_Figure_svg.Rmd`</pre>


## Data Availability

All spatial transcriptomics (10x Visium) used for this study is publicly available on GEO under the accession token  [GSE224360][2].

[2]:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224360

1. QC, Normalization and Subset: `120123_Normalisation_subset_merge_paper.Rmd`
    * Input file (assumes path: "")
    * Reads in the output of the spaceranger pipeline, and returns a Seurat object. 
    * Removal of connective tissue and joint dimensional reduction to cluster on the underlying RNA expression data.

2. Differential gene expression: `120123_DGE_Table_Paper.Rmd`
    * Takes merged file as input and filters for state/cluster to compere.
    * FindMarkers() function is used to generate .csv files of idents of interests.


3. Figures from the paper: `120123_Figure_svg.Rmd`
    * Takes merged file as input 
    * Uses .xlsx file as input for figuer xx

Figures: XXX Were created by shiny app

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

other attached packages:
 [1] cowplot_1.1.1      viridis_0.6.2      viridisLite_0.4.1  readxl_1.4.1       ggrepel_0.9.2      rmarkdown_2.19     dplyr_1.0.10      
 [8] patchwork_1.1.2    ggplot2_3.4.0      SeuratObject_4.0.1 Seurat_4.0.0      

loaded via a namespace (and not attached):
  [1] Rtsne_0.16           colorspace_2.0-3     deldir_1.0-6         ellipsis_0.3.2       ggridges_0.5.4       rstudioapi_0.14     
  [7] spatstat.data_3.0-0  farver_2.1.1         leiden_0.4.3         listenv_0.9.0        bit64_4.0.5          fansi_1.0.3         
 [13] codetools_0.2-18     splines_4.1.0        knitr_1.41           polyclip_1.10-4      jsonlite_1.8.3       ica_1.0-3           
 [19] cluster_2.1.2        png_0.1-8            uwot_0.1.14          shiny_1.7.4          sctransform_0.3.2    BiocManager_1.30.19 
 [25] compiler_4.1.0       httr_1.4.4           assertthat_0.2.1     Matrix_1.5-1         fastmap_1.1.0        lazyeval_0.2.2      
 [31] limma_3.48.3         cli_3.5.0            later_1.3.0          htmltools_0.5.4      tools_4.1.0          igraph_1.3.5        
 [37] gtable_0.3.1         glue_1.6.2           RANN_2.6.1           reshape2_1.4.4       Rcpp_1.0.9           spatstat_1.64-1     
 [43] scattermore_0.8      cellranger_1.1.0     vctrs_0.5.0          nlme_3.1-152         lmtest_0.9-40        xfun_0.34           
 [49] stringr_1.5.0        globals_0.16.2       mime_0.12            miniUI_0.1.1.1       lifecycle_1.0.3      irlba_2.3.5.1       
 [55] renv_0.16.0          goftest_1.2-3        future_1.30.0        MASS_7.3-54          zoo_1.8-11           scales_1.2.1        
 [61] promises_1.2.0.1     spatstat.utils_3.0-1 parallel_4.1.0       RColorBrewer_1.1-3   yaml_2.3.6           reticulate_1.26     
 [67] pbapply_1.6-0        gridExtra_2.3        rpart_4.1-15         stringi_1.7.8        rlang_1.0.6          pkgconfig_2.0.3     
 [73] matrixStats_0.62.0   evaluate_0.19        lattice_0.20-44      ROCR_1.0-11          purrr_0.3.5          tensor_1.5          
 [79] labeling_0.4.2       htmlwidgets_1.6.0    bit_4.0.5            tidyselect_1.2.0     parallelly_1.33.0    RcppAnnoy_0.0.20    
 [85] plyr_1.8.7           magrittr_2.0.3       R6_2.5.1             generics_0.1.3       DBI_1.1.3            pillar_1.8.1        
 [91] withr_2.5.0          mgcv_1.8-35          fitdistrplus_1.1-8   survival_3.2-11      abind_1.4-5          tibble_3.1.8        
 [97] future.apply_1.10.0  crayon_1.5.2         hdf5r_1.3.7          KernSmooth_2.23-20   utf8_1.2.2           plotly_4.10.1       
[103] grid_4.1.0           data.table_1.14.4    digest_0.6.30        xtable_1.8-4         tidyr_1.2.1          httpuv_1.6.6        
[109] munsell_0.5.0 </pre>

Usage of more recent version of certain R packages can lead to difference in results especially `Seurat` and `sctransform`.

## Citation 
[SPATIO-TEMPORAL_RECRUITMENT_OF_ADULT_NEURAL_STEM_CELLS_FOR_TRANSIENT_NEUROGENESIS_DURING_PREGNANCY] [1]








# To do:

*change names of flies

*Set file directery apropriet/ change save files 

*Comment better

*Insert link to paper 

*App
