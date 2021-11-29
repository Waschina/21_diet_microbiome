# 21_diet_microbiome
Metabolic modeling of microbial community metabolism (Diet-Microbiome-Study 2021)

## Prerequisites 

#### Data (model)

[gapseq](https://doi.org/10.1186/s13059-021-02295-1) reconstructions can be obtained from zenodo ([![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5735482.svg)](https://doi.org/10.5281/zenodo.5735482)). Extract the archive to your local PC and adjust the path in the script `R/00_1_prepare_MAG_info_models.R` (line 38 and 42) to point towards your 'models' directory.

#### R and R-packages

The simulations are performed in R (v 4.1.1 or higher) and require the following packages "data.table", "stringr", "sybil", "cplexAPI", "zoo", "lme4", "ggplot2", "parallel", and "ggtext".

For specific package version numbers that were for the analysis presented in the manuscript, please see below under "R session information".

In addition the in-house R-package "MicrobiomeGS2"  is needed. The source package archive can be downloaded [here](https://www.nutrinf.uni-kiel.de/MicrobiomeGS2_0.1.0.tar.gz).

#### Solver

The simulations require the CPLEX solver included in the *IBM ILOG CPLEX Optimization Studio*. 

## Performing metabolic simulations

Initializing input data, meta information, and models:

```R
# prepare diet info data
source("analysis/v1/scripts/00_1_prepare_Diet_info.R")

# prepare MAG metabolic models and their isolation phenotypes + meta info
source("analysis/v1/scripts/00_1_prepare_MAG_info_models.R")

# sample meta data
source("analysis/v1/scripts/00_1_prepare_metadata.R")
```

Now, we can start the time-consuming simulations of microbial community metabolism: (~45 minutes).

```R
# Community FBA
source("analysis/v1/scripts/03_communityFBA_fluxes.R")
```

And it's done.

## R session information

```
R version 4.1.1 (2021-08-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.3 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggtext_0.1.1        ggplot2_3.3.5       lme4_1.1-27.1      
 [4] zoo_1.8-9           cplexAPI_1.4.0      MicrobiomeGS2_0.1.1
 [7] sybil_2.2.0         lattice_0.20-45     Matrix_1.3-4       
[10] stringr_1.4.0       data.table_1.14.2  

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.7       compiler_4.1.1   pillar_1.6.4     nloptr_1.2.2.3  
 [5] tools_4.1.1      digest_0.6.28    boot_1.3-28      lifecycle_1.0.1 
 [9] tibble_3.1.6     gtable_0.3.0     nlme_3.1-153     pkgconfig_2.0.3 
[13] rlang_0.4.12     CHNOSZ_1.4.1     parallel_4.1.1   Rdisop_1.52.0   
[17] xml2_1.3.2       withr_2.4.2      dplyr_1.0.7      generics_0.1.1  
[21] vctrs_0.3.8      gridtext_0.1.4   tidyselect_1.1.1 grid_4.1.1      
[25] glue_1.5.0       R6_2.5.1         fansi_0.5.0      minqa_1.2.4     
[29] farver_2.1.0     purrr_0.3.4      magrittr_2.0.1   scales_1.1.1    
[33] ellipsis_0.3.2   MASS_7.3-54      splines_4.1.1    colorspace_2.0-2
[37] labeling_0.4.2   utf8_1.2.2       stringi_1.7.5    munsell_0.5.0   
[41] markdown_1.1     crayon_1.4.2 
```

