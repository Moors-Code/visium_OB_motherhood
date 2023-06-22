# Sensitivity analysis 



<pre>`sensitivity_check.R`</pre>


## Code description

1. `sensitivity_check.R`
   

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
[1] SeuratObject_4.1.3 sp_1.6-0          

loaded via a namespace (and not attached):
  [1] nlme_3.1-152           matrixStats_0.63.0     spatstat.sparse_3.0-1  RcppAnnoy_0.0.20       RColorBrewer_1.1-3     httr_1.4.5            
  [7] sctransform_0.3.5      tools_4.1.0            utf8_1.2.3             R6_2.5.1               irlba_2.3.5.1          KernSmooth_2.23-20    
 [13] uwot_0.1.14            lazyeval_0.2.2         colorspace_2.1-0       tidyselect_1.2.0       gridExtra_2.3          compiler_4.1.0        
 [19] progressr_0.13.0       cli_3.6.1              spatstat.explore_3.1-0 plotly_4.10.1          Seurat_4.3.0           scales_1.2.1          
 [25] lmtest_0.9-40          spatstat.data_3.0-1    ggridges_0.5.4         pbapply_1.7-0          goftest_1.2-3          stringr_1.5.0         
 [31] digest_0.6.31          spatstat.utils_3.0-2   pkgconfig_2.0.3        htmltools_0.5.5        parallelly_1.35.0      fastmap_1.1.1         
 [37] htmlwidgets_1.6.2      rlang_1.1.0            rstudioapi_0.13        shiny_1.7.4            generics_0.1.3         zoo_1.8-12            
 [43] jsonlite_1.8.4         spatstat.random_3.1-4  ica_1.0-3              dplyr_1.1.2            magrittr_2.0.3         patchwork_1.1.2       
 [49] Matrix_1.5-4           Rcpp_1.0.10            munsell_0.5.0          fansi_1.0.4            abind_1.4-5            reticulate_1.28       
 [55] lifecycle_1.0.3        stringi_1.7.12         MASS_7.3-54            Rtsne_0.16             plyr_1.8.8             grid_4.1.0            
 [61] parallel_4.1.0         listenv_0.9.0          promises_1.2.0.1       ggrepel_0.9.3          deldir_1.0-6           miniUI_0.1.1.1        
 [67] lattice_0.20-44        cowplot_1.1.1          splines_4.1.0          tensor_1.5             pillar_1.9.0           igraph_1.4.2          
 [73] spatstat.geom_3.1-0    future.apply_1.10.0    reshape2_1.4.4         codetools_0.2-18       leiden_0.4.3           glue_1.6.2            
 [79] data.table_1.14.8      renv_0.16.0            png_0.1-8              vctrs_0.6.2            httpuv_1.6.9           polyclip_1.10-4       
 [85] gtable_0.3.3           RANN_2.6.1             purrr_1.0.1            tidyr_1.3.0            scattermore_0.8        future_1.32.0         
 [91] ggplot2_3.4.2          mime_0.12              xtable_1.8-4           later_1.3.0            survival_3.2-11        viridisLite_0.4.1     
 [97] tibble_3.2.1           cluster_2.1.2          globals_0.16.2         fitdistrplus_1.1-11    ellipsis_0.3.2         ROCR_1.0-11                       
</pre>





