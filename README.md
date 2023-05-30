# TransCistor

TransCistor is an analytical framework which aims to classify cis- or trans-acting long non-coding RNAs(lncRNAs). It includes two modules; digital and analogue. There is only one file required as input: a processed whole-transcriptome list of target/not-target genes (“regulation file”).
https://www.biorxiv.org/content/10.1101/2022.09.18.508380v1
## Description

TransCistor-digital tests for statistical enrichment in proximal genes using the hypergeometric distribution. Proximity is defined based on Topologically Associating Domains (TADs). \
TransCistor-analogue, defines a distance statistic as the mean TSS-to-TSS distance of all same-chromosome targets of a given lncRNA. To estimate statistical significance, a null distribution is calculated by randomisation of target labels.

## Getting Started

### Dependencies

TransCistor is built as a webserver (https://transcistor.unibe.ch/) or standalone script.
* The webserver is a shiny app and has no dependencies.
* The standalone version is built in R and the main dependencies are R packages: tidyverse, dplyr, rlang, rlist, etc.


### Executing/Installing TransCistor

* For single-case studies we highly recommend to use the webserver version of TransCistor (https://transcistor.unibe.ch/). 
* Standalone version is recommended for analysing multiple datasets. Executing simply requires the start_TransCistor.R script. The tested dataset can be changed by manipulating the metadata and RegulationFiles directory

``` 
TransCistor <- function(input.file, id.type, species, cell = "H1-NPC",
                        lncRNA.name, lncRNA.chr, lncRNA.tss, lncRNA.strand = "+",
                        enricher.threshold = 0, simulations = 100000 )
```

You can use the plots.R file to get results for single lncRNAs.  You may perform the analysis as shown below. It should take a few seconds in a normal computer. \

``` 
### Input the Identifier for the lncRNA####
input='UMLILO'

####Read the metadata to get the regulation file path. Either read the given metadata or make a custom one based on the given format#### 

metadata_file <- 'metadata.txt'
metadata <- read.table(metadata_file, header = TRUE, sep = "\t",na.strings=c("","NA"))
file_name <- metadata[metadata$SYMBOL==input,]$File_Name

name <- metadata[i[num],]$SYMBOL
chr <- metadata[i[num],]$CHR
tss <- metadata[i[num],]$TSS
strand <- metadata[i[num],]$strand
file.type <- metadata[i[num],]$FILE_TYPE
species <- metadata[i[num],]$SPECIES
#TransCistor function call
results <- TransCistor(input.file = current_file, 
                         id.type = file.type, 
                         species = species, 
                         lncRNA.name = name,
                         lncRNA.chr = chr, 
                         lncRNA.tss = tss, 
                         lncRNA.strand = strand, 
                         enricher.threshold = enrichr_threshold
  )

```

...or visit the webserver to get publucation ready figures!

## Main contributors' names and contact info. 

(more people have contributed in this work, see https://www.biorxiv.org/content/10.1101/2022.09.18.508380v1)

Panagiotis Chouvardas, panagiotis.chouvardas@dbmr.unibe.ch <br>
Bhavya Dhaka, bhavya.dhaka@ucdconnect.ie  <br>
Rory Johnson*, rory.johnson@dbmr.unibe.ch

### Session Info

```
R version 4.1.2 (2021-11-01)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.1 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0

locale:
 [1] LC_CTYPE=en_IE.UTF-8      
 [2] LC_NUMERIC=C              
 [3] LC_TIME=en_IE.UTF-8       
 [4] LC_COLLATE=en_IE.UTF-8    
 [5] LC_MONETARY=en_IE.UTF-8   
 [6] LC_MESSAGES=en_IE.UTF-8   
 [7] LC_PAPER=en_IE.UTF-8      
 [8] LC_NAME=C                 
 [9] LC_ADDRESS=C              
[10] LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_IE.UTF-8
[12] LC_IDENTIFICATION=C       

attached base packages:
[1] grid      stats     graphics 
[4] grDevices utils     datasets 
[7] methods   base     

other attached packages:
 [1] scales_1.2.1         
 [2] ggh4x_0.2.3          
 [3] reshape_0.8.9        
 [4] VennDiagram_1.7.3    
 [5] futile.logger_1.4.3  
 [6] plyr_1.8.7           
 [7] pCalibrate_0.2-1     
 [8] MCMCpack_1.6-3       
 [9] MASS_7.3-55          
[10] coda_0.19-4          
[11] exact2x2_1.6.6       
[12] exactci_1.4-2        
[13] testthat_3.1.5       
[14] ssanv_1.1            
[15] psych_2.2.5          
[16] RColorBrewer_1.1-3   
[17] pheatmap_1.0.12      
[18] ggpubr_0.5.0         
[19] clusterProfiler_4.2.2
[20] hash_2.2.6.2         
[21] forcats_0.5.2        
[22] stringr_1.4.1        
[23] purrr_0.3.4          
[24] readr_2.1.2          
[25] tidyr_1.2.1          
[26] tibble_3.1.8         
[27] tidyverse_1.3.2      
[28] ggrepel_0.9.1        
[29] ggplot2_3.4.0        
[30] readxl_1.4.1         
[31] data.table_1.14.2    
[32] rlist_0.4.6.2        
[33] dplyr_1.0.10         
[34] rlang_1.0.6          

loaded via a namespace (and not attached):
  [1] shadowtext_0.1.2      
  [2] backports_1.4.1       
  [3] fastmatch_1.1-3       
  [4] igraph_1.3.5          
  [5] lazyeval_0.2.2        
  [6] splines_4.1.2         
  [7] BiocParallel_1.28.3   
  [8] GenomeInfoDb_1.30.1   
  [9] digest_0.6.29         
 [10] yulab.utils_0.0.6     
 [11] GOSemSim_2.20.0       
 [12] viridis_0.6.2         
 [13] GO.db_3.14.0          
 [14] fansi_1.0.3           
 [15] magrittr_2.0.3        
 [16] memoise_2.0.1         
 [17] googlesheets4_1.0.1   
 [18] tzdb_0.3.0            
 [19] Biostrings_2.62.0     
 [20] graphlayouts_0.8.4    
 [21] modelr_0.1.9          
 [22] enrichplot_1.14.2     
 [23] colorspace_2.0-3      
 [24] blob_1.2.3            
 [25] rvest_1.0.3           
 [26] haven_2.5.1           
 [27] crayon_1.5.1          
 [28] RCurl_1.98-1.9        
 [29] jsonlite_1.8.0        
 [30] scatterpie_0.1.8      
 [31] survival_3.2-13       
 [32] ape_5.6-2             
 [33] glue_1.6.2            
 [34] polyclip_1.10-4       
 [35] gtable_0.3.1          
 [36] gargle_1.2.1          
 [37] zlibbioc_1.40.0       
 [38] XVector_0.34.0        
 [39] MatrixModels_0.5-1    
 [40] car_3.1-1             
 [41] BiocGenerics_0.40.0   
 [42] SparseM_1.81          
 [43] abind_1.4-5           
 [44] DOSE_3.20.1           
 [45] futile.options_1.0.1  
 [46] DBI_1.1.3             
 [47] rstatix_0.7.1         
 [48] Rcpp_1.0.9            
 [49] viridisLite_0.4.1     
 [50] gridGraphics_0.5-1    
 [51] tidytree_0.4.2        
 [52] bit_4.0.4             
 [53] stats4_4.1.2          
 [54] httr_1.4.4            
 [55] fgsea_1.20.0          
 [56] ellipsis_0.3.2        
 [57] pkgconfig_2.0.3       
 [58] farver_2.1.1          
 [59] dbplyr_2.2.1          
 [60] utf8_1.2.2            
 [61] labeling_0.4.2        
 [62] ggplotify_0.1.0       
 [63] tidyselect_1.2.0      
 [64] reshape2_1.4.4        
 [65] AnnotationDbi_1.56.2  
 [66] munsell_0.5.0         
 [67] cellranger_1.1.0      
 [68] tools_4.1.2           
 [69] cachem_1.0.6          
 [70] downloader_0.4        
 [71] cli_3.4.0             
 [72] generics_0.1.3        
 [73] RSQLite_2.2.18        
 [74] broom_1.0.1           
 [75] fastmap_1.1.0         
 [76] mcmc_0.9-7            
 [77] ggtree_3.2.1          
 [78] bit64_4.0.5           
 [79] fs_1.5.2              
 [80] tidygraph_1.2.2       
 [81] KEGGREST_1.34.0       
 [82] ggraph_2.1.0          
 [83] nlme_3.1-155          
 [84] quantreg_5.94         
 [85] formatR_1.12          
 [86] aplot_0.1.9           
 [87] DO.db_2.9             
 [88] xml2_1.3.3            
 [89] brio_1.1.3            
 [90] compiler_4.1.2        
 [91] rstudioapi_0.14       
 [92] png_0.1-7             
 [93] ggsignif_0.6.4        
 [94] reprex_2.0.2          
 [95] treeio_1.18.1         
 [96] tweenr_2.0.2          
 [97] stringi_1.7.8         
 [98] lattice_0.20-45       
 [99] Matrix_1.5-1          
[100] vctrs_0.5.1           
[101] pillar_1.8.1          
[102] lifecycle_1.0.3       
[103] bitops_1.0-7          
[104] patchwork_1.1.2       
[105] qvalue_2.26.0         
[106] R6_2.5.1              
[107] gridExtra_2.3         
[108] IRanges_2.28.0        
[109] lambda.r_1.2.4        
[110] assertthat_0.2.1      
[111] withr_2.5.0           
[112] mnormt_2.1.0          
[113] S4Vectors_0.32.4      
[114] GenomeInfoDbData_1.2.7
[115] parallel_4.1.2        
[116] hms_1.1.2             
[117] ggfun_0.0.9           
[118] carData_3.0-5         
[119] googledrive_2.0.0     
[120] ggforce_0.4.1         
[121] Biobase_2.54.0        
[122] lubridate_1.8.0    
```



