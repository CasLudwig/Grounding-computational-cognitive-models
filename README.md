# Grounding computational cognitive models

Repository accompanying [Ludwig, Stuchlý, & Malhotra](https://osf.io/preprints/psyarxiv/vur6t). 

- Code for generating Figure 4 (mixture of Gaussians example).
- Code for generating Figures 5 and C1 (simulated expanded judgement paradigm).

## Mixture of Gaussians

The [mixGaussiansPF](https://github.com/CasLudwig/Grounding-computational-cognitive-models/tree/main/mixGaussiansPF) folder contains all the code for generating Figure 4 in the paper. Assuming you have all the relevant libraries and Stan working together with R (see below for setup info), you should simply be able to run [timeVarNormalSimulation.R](https://github.com/CasLudwig/Grounding-computational-cognitive-models/blob/main/mixGaussiansPF/timeVarNormalSimulation.R). 

## Simulated expanded judgement paradigm

The [expandedJudgementSimulations](https://github.com/CasLudwig/Grounding-computational-cognitive-models/tree/main/expandedJudgementSimulations) folder contains code and data for generating Figures 5 and C1 (Appendix C) in the paper. Running this code is a bit more involved, but we have tried to make things easier by providing a [RNotebook](https://github.com/CasLudwig/Grounding-computational-cognitive-models/blob/main/expandedJudgementSimulations/simulateExpJudgementData.Rmd) with extensive commentary. You do not need to run the code in this notebook, but it provides insight in how the (simulated) data reported in the paper were generated. The notebook contains all the instructions for recreating the figures, but in brief:

- Figure 5: run [trajectoryExpJudgement.R](https://github.com/CasLudwig/Grounding-computational-cognitive-models/blob/main/expandedJudgementSimulations/trajectoryExpJudgement.R).
- Figure C1: run [rewardRatesExpJudgement.R](https://github.com/CasLudwig/Grounding-computational-cognitive-models/blob/main/expandedJudgementSimulations/rewardRatesExpJudgement.R).

## Disclaimer

If the code does not work because your setup is different from mine (see below), there is not much I can do about this although I'm happy to try and help. If the code does not work properly because there are errors, please let me know and I will of course fix as soon as I can.

## Environment

Code was written in the following environment. Not all packages below are loaded all the time, but they are all the ones that are loaded at some point or another.

R version 4.4.0 (2024-04-24)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.6.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/London
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1      purrr_1.0.2        tidyr_1.3.1        tibble_3.2.1       tidyverse_2.0.0   
 [8] tmvtnorm_1.6       gmm_1.8            sandwich_3.1-0     Matrix_1.7-0       mvtnorm_1.2-5      Hmisc_5.1-3        bayesplot_1.11.1  
[15] coda_0.19-4.1      rstan_2.32.6       StanHeaders_2.32.8 latex2exp_0.9.6    gridExtra_2.3      reshape_0.8.9      viridis_0.6.5     
[22] viridisLite_0.4.2  dplyr_1.1.4        readr_2.1.5        ggplot2_3.5.1     

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1   loo_2.7.0          fastmap_1.2.0      digest_0.6.35      rpart_4.1.23       timechange_0.3.0   lifecycle_1.0.4   
 [8] cluster_2.1.6      magrittr_2.0.3     compiler_4.4.0     rlang_1.1.3        tools_4.4.0        utf8_1.2.4         yaml_2.3.8        
[15] data.table_1.15.4  knitr_1.47         ggsignif_0.6.4     htmlwidgets_1.6.4  pkgbuild_1.4.4     plyr_1.8.9         abind_1.4-5       
[22] withr_3.0.0        foreign_0.8-86     nnet_7.3-19        grid_4.4.0         fansi_1.0.6        ggpubr_0.6.0       colorspace_2.1-0  
[29] inline_0.3.19      scales_1.3.0       cli_3.6.2          rmarkdown_2.27     generics_0.1.3     RcppParallel_5.1.7 rstudioapi_0.16.0 
[36] tzdb_0.4.0         parallel_4.4.0     matrixStats_1.3.0  base64enc_0.1-3    vctrs_0.6.5        jsonlite_1.8.8     carData_3.0-5     
[43] car_3.1-2          hms_1.1.3          rstatix_0.7.2      Formula_1.2-5      htmlTable_2.4.2    glue_1.7.0         codetools_0.2-20  
[50] stringi_1.8.4      gtable_0.3.5       QuickJSR_1.1.3     munsell_0.5.1      pillar_1.9.0       htmltools_0.5.8.1  R6_2.5.1          
[57] evaluate_0.23      lattice_0.22-6     backports_1.5.0    broom_1.0.6        Rcpp_1.0.12        checkmate_2.3.1    xfun_0.44         
[64] zoo_1.8-12         pkgconfig_2.0.3

Stan version 2.32.2 (May 2023)
