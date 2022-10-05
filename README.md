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
* The standalone version is built in R and the main dependencies are R packages: tidyverse, dplyr, rlang, rlist.


### Executing/Installing TransCistor

* For single-case studies we highly recomend to use the webserver version of TransCistor (https://transcistor.unibe.ch/). 
* Standalone version is recommended for analysing multiple datasets. Executing simply requires "sourcing" the TransCistor R script.

``` 
source("transcistor.R")
TransCistor(input.file, 			#input regulation file
	id.type, 				#gene id type: ENSEMBL, SYMBOL or ENTREZ
	species, 				#human or mouse
        lncRNA.name,	 			#only used for visualization purposes
	lncRNA.chr,	  			#chromosome in chr1, chr2, ... format
	lncRNA.tss, 				#IMPORTANT: hg38 or mm10 coordinates	
	lncRNA.strand = "+",			#only used for visualization purposes
        enricher.threshold = 0, 		#currently not used.
	simulations = 1000	 		#analogue number of simulations. 
			)
```

An example regulation file is provided (UMLILO.txt). You may perform the analysis as shown below. It should take a few seconds in a normal computer. \
The results is a list of 8 data frames. The most important ones are [[4]] and [[9]] which included the final pvalues and the list of genes in the same TAD.

``` 
df <- read.file("UMLILO.txt", sep="\t", header = F)
results <- TransCistor(input.file = df, 			
              id.type = "ENSEMBL", 				
              species = "human", 				
              lncRNA.name = "UMLILO", 			
              lncRNA.chr = "chr4",  			
              lncRNA.tss =  73710302, 				
              lncRNA.strand = "+",	
              enricher.threshold = 0, 
              simulations = 1000 		 
  )

```

...or visit the webserver to get publucation ready figures!

## Authors

Main contributors names and contact info. (more people have contributed in this work, see https://www.biorxiv.org/content/10.1101/2022.09.18.508380v1)

Panagiotis Chouvardas, panagiotis.chouvardas@dbmr.unibe.ch
Rory Johnson, rory.johnson@dbmr.unibe.ch

### Session Info

```
R version 4.1.2 (2021-11-01)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044)

Matrix products: default

locale:
[1] LC_COLLATE=English_Switzerland.65001 
[2] LC_CTYPE=C                           
[3] LC_MONETARY=English_Switzerland.65001
[4] LC_NUMERIC=C                         
[5] LC_TIME=English_Switzerland.65001    
system code page: 65001

attached base packages:
[1] stats     graphics  grDevices utils    
[5] datasets  methods   base     

other attached packages:
 [1] psych_2.2.9           RColorBrewer_1.1-2   
 [3] pheatmap_1.0.12       ggpubr_0.4.0         
 [5] clusterProfiler_4.2.2 forcats_0.5.1        
 [7] stringr_1.4.0         purrr_0.3.4          
 [9] readr_2.1.2           tidyr_1.2.0          
[11] tibble_3.1.6          ggplot2_3.3.5        
[13] tidyverse_1.3.1       rlist_0.4.6.2        
[15] dplyr_1.0.8           rlang_1.0.1          

loaded via a namespace (and not attached):
  [1] utf8_1.2.2             reticulate_1.24       
  [3] tidyselect_1.1.1       RSQLite_2.2.10        
  [5] AnnotationDbi_1.56.2   htmlwidgets_1.5.4     
  [7] FactoMineR_2.4         grid_4.1.2            
  [9] BiocParallel_1.28.3    Rtsne_0.15            
 [11] scatterpie_0.1.7       munsell_0.5.0         
 [13] codetools_0.2-18       ica_1.0-2             
 [15] DT_0.22                future_1.27.0         
 [17] miniUI_0.1.1.1         withr_2.4.3           
 [19] spatstat.random_2.1-0  colorspace_2.0-2      
 [21] GOSemSim_2.20.0        Biobase_2.54.0        
 [23] knitr_1.37             rstudioapi_0.13       
 [25] leaps_3.1              Seurat_4.1.0          
 [27] stats4_4.1.2           ROCR_1.0-11           
 [29] ggsignif_0.6.3         tensor_1.5            
 [31] DOSE_3.20.1            listenv_0.8.0         
 [33] GenomeInfoDbData_1.2.7 mnormt_2.0.2          
 [35] polyclip_1.10-0        farver_2.1.0          
 [37] bit64_4.0.5            downloader_0.4        
 [39] rprojroot_2.0.2        treeio_1.18.1         
 [41] parallelly_1.32.1      vctrs_0.3.8           
 [43] generics_0.1.2         xfun_0.29             
 [45] R6_2.5.1               GenomeInfoDb_1.30.1   
 [47] graphlayouts_0.8.0     pals_1.7              
 [49] gridGraphics_0.5-1     bitops_1.0-7          
 [51] spatstat.utils_2.3-0   cachem_1.0.6          
 [53] fgsea_1.20.0           assertthat_0.2.1      
 [55] promises_1.2.0.1       scales_1.1.1          
 [57] ggraph_2.0.5           enrichplot_1.14.1     
 [59] gtable_0.3.0           globals_0.16.0        
 [61] goftest_1.2-3          tidygraph_1.2.0       
 [63] scatterplot3d_0.3-41   splines_4.1.2         
 [65] rstatix_0.7.0          lazyeval_0.2.2        
 [67] dichromat_2.0-0        spatstat.geom_2.3-2   
 [69] broom_0.7.12           yaml_2.2.2            
 [71] reshape2_1.4.4         abind_1.4-5           
 [73] modelr_0.1.8           backports_1.4.1       
 [75] httpuv_1.6.5           qvalue_2.26.0         
 [77] tools_4.1.2            ggplotify_0.1.0       
 [79] ellipsis_0.3.2         spatstat.core_2.4-0   
 [81] BiocGenerics_0.40.0    ggridges_0.5.3        
 [83] Rcpp_1.0.8             plyr_1.8.6            
 [85] zlibbioc_1.40.0        RCurl_1.98-1.6        
 [87] rpart_4.1.16           deldir_1.0-6          
 [89] viridis_0.6.2          pbapply_1.5-0         
 [91] cowplot_1.1.1          S4Vectors_0.32.3      
 [93] zoo_1.8-9              SeuratObject_4.0.4    
 [95] haven_2.4.3            ggrepel_0.9.1         
 [97] cluster_2.1.2          fs_1.5.2              
 [99] factoextra_1.0.7       magrittr_2.0.2        
[101] data.table_1.14.2      scattermore_0.8       
[103] DO.db_2.9              lmtest_0.9-39         
[105] reprex_2.0.1           RANN_2.6.1            
[107] tmvnsim_1.0-2          fitdistrplus_1.1-6    
[109] matrixStats_0.61.0     hms_1.1.1             
[111] patchwork_1.1.1        mime_0.12             
[113] evaluate_0.15          xtable_1.8-4          
[115] readxl_1.3.1           IRanges_2.28.0        
[117] gridExtra_2.3          compiler_4.1.2        
[119] maps_3.4.0             shadowtext_0.1.1      
[121] KernSmooth_2.23-20     crayon_1.5.0          
[123] htmltools_0.5.2        ggfun_0.0.5           
[125] mgcv_1.8-38            later_1.3.0           
[127] tzdb_0.2.0             aplot_0.1.2           
[129] lubridate_1.8.0        DBI_1.1.2             
[131] tweenr_1.0.2           dbplyr_2.1.1          
[133] MASS_7.3-55            car_3.0-12            
[135] Matrix_1.4-0           cli_3.1.1             
[137] parallel_4.1.2         igraph_1.2.11         
[139] pkgconfig_2.0.3        flashClust_1.01-2     
[141] plotly_4.10.0          spatstat.sparse_2.1-0 
[143] xml2_1.3.3             ggtree_3.2.1          
[145] XVector_0.34.0         rvest_1.0.2           
[147] yulab.utils_0.0.4      digest_0.6.29         
[149] sctransform_0.3.3      RcppAnnoy_0.0.19      
[151] spatstat.data_2.1-2    Biostrings_2.62.0     
[153] rmarkdown_2.11         cellranger_1.1.0      
[155] leiden_0.3.9           enrichR_3.0           
[157] fastmatch_1.1-3        tidytree_0.3.8        
[159] uwot_0.1.11            shiny_1.7.1           
[161] rjson_0.2.21           lifecycle_1.0.1       
[163] nlme_3.1-155           jsonlite_1.7.3        
[165] carData_3.0-5          mapproj_1.2.8         
[167] viridisLite_0.4.0      fansi_1.0.2           
[169] pillar_1.7.0           lattice_0.20-45       
[171] KEGGREST_1.34.0        fastmap_1.1.0         
[173] httr_1.4.2             survival_3.2-13       
[175] GO.db_3.14.0           glue_1.6.1            
[177] png_0.1-7              bit_4.0.4             
[179] ggforce_0.3.3          stringi_1.7.6         
[181] blob_1.2.2             memoise_2.0.1         
[183] ape_5.6-1              irlba_2.3.5           
[185] future.apply_1.8.1    
```



