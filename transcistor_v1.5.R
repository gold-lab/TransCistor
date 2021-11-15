TransCistor <- function(input.file, id.type, species, 
                        lncRNA.name, lncRNA.chr, lncRNA.tss, lncRNA.strand = "+",
                        enricher.threshold = 0.05){
  
  #to do
  #1. remake gencode matrix to be sure and to make it more clean
  #2. re-download tad to have a correct link
  #3. add toupper for mouse
  #4. add plots in the script
  #5. increase number of iterations for distance
  #6. run with all data
  #7. update webserver with new version
  #8. add options: to run part of the analysis or all, number of simulations,...
  #9. add step to remove perturbed lncRNA from targets (a couple of FANTOM ones are pc)
  
  library(rlang)
  library(dplyr)
  library(rlist)
  library(tidyverse)
  library(clusterProfiler)
  
  options(scipen = 999)
  `%!in%` = Negate(`%in%`)
  
  #TEST PARAMS - REMOVE WHEN DONE#######
  #input.file = "ASO_C0000830_02.txt"
  #id.type="ensembl"
  #species="human"
  #enricher.threshold = 0.05
  #lncRNA.chr <- "chr1"
  #lncRNA.tss <- 50960745
  #input.file <- read.table(input.file, sep="\t")
  
  
  #############################
  #   Input/Reference files   #
  #############################
  
  #reading gencode and tad files
  species = tolower(species)
  if((species == "human") || (species == "mouse")){
    gencode <- read.table(paste0("References/", species, "_gencode.txt"), sep="\t", header=T, stringsAsFactors=FALSE)
    tad <- read.table(paste0("References/", species, "_tad.txt"), sep="\t", header=F, stringsAsFactors=FALSE)
  }else{
    print("species should be either human or mouse.")
    return()
  }
  
  #REMOVE WHEN YOU FIX GENCODE INPUT##################
  gencode=gencode[,c("ENSEMBL", "SYMBOL","ENTREZID", "CHR", "start","end", "strand", "TSS","TYPE")]
  colnames(gencode) <- toupper(colnames(gencode))
  colnames(tad) <- c("CHR", "START", "END")
  ####################################################
  
  #read input file and merge with gencode metadata
  df.raw <- input.file
  
  id.type <- toupper(id.type)
  colnames(df.raw) <- c(id.type, "Regulation")
  
  if((id.type == "ENSEMBL") || (id.type =="ENTREZ") || (id.type == "SYMBOL")){
    df.raw <- merge(df.raw, gencode, by=id.type)  
  }else{
    print("id.type should be either ENSEMBL, ENTREZ or SYMBOL")
    return()
  }
  
  
  #############################
  #      Pre-processing       #
  #############################
  #remove non protein coding genes
  df.protein <- filter(df.raw, TYPE == "protein_coding")
  df.protein$SYMBOL <- toupper(df.protein$SYMBOL)
  df.protein.targets <- filter(df.protein, Regulation != 0)
  
  if(nrow(df.protein.targets) == 0){
    print(paste0("No targets found. exiting"))
    return()
  }
  
  #secondary target filtering
  TFs <- read.gmt("References/TF_Perturbations_Followed_by_Expression.txt")  
  if(length(which(df.protein.targets$SYMBOL %in% TFs$gene)) > 0){
    TFs <- enricher(df.protein.targets$SYMBOL, TERM2GENE = TFs, qvalueCutoff = enricher.threshold,
                    minGSSize = 1, maxGSSize = 5000)
    
    if(nrow(TFs) > 0){
      Enriched.TFs <- data.frame(TFs=sapply(strsplit(TFs$ID, " "), "[", 1), targets = TFs$geneID)
      Enriched.TFs.DE <- Enriched.TFs[which(Enriched.TFs$TFs %in% df.protein.targets$SYMBOL),]
      if(nrow(Enriched.TFs.DE) > 0){
        Secondary.Targets = paste(Enriched.TFs.DE$targets, collapse="/")
        Secondary.Targets = unique(strsplit(Secondary.Targets, split="/")[[1]])
        df.protein.targets.filtered <- filter(df.protein.targets, SYMBOL %!in% Secondary.Targets)
        print(paste("Filtered", length(Secondary.Targets), "secondary targets"))
      }else{
        print("No secondary targets filtered.")
        df.protein.targets.filtered <- df.protein.targets
      }
    }else{
      print("No secondary targets filtered.")
      df.protein.targets.filtered <- df.protein.targets
    }
  }else{
    print("No secondary targets filtered.")
    df.protein.targets.filtered <- df.protein.targets
  }
  
  #Initializing results tables
  Results <- data.frame(Method=c("TAD","EMPIRICAL", paste0("W",seq(500000, 5000000, 500000))),
                                 Proximal=NA, Distal=NA, Proximal.Targets.Activated=NA, Total.Targets.Activated=NA,
                                 Proximal.Targets.Repressed=NA, Total.Targets.Repressed=NA,
                                 Acticator.Pvalue=NA, Repressor.Pvalue=NA)
  
  #############################
  #            TAD            #
  #############################
  lncRNA.tad <- filter(tad,CHR == lncRNA.chr & START <= lncRNA.tss & END >= lncRNA.tss)
  if(nrow(lncRNA.tad) == 0){
    print("No ovalapping TAD found. NAs reported for TAD window.")
  }else{
    tad.proximal.genes <- filter(df.protein, CHR == lncRNA.tad$CHR & TSS > lncRNA.tad$START & TSS < lncRNA.tad$END)
    tad.proximal.targets <- filter(tad.proximal.genes, SYMBOL %in% df.protein.targets.filtered$SYMBOL)
    
    #estimate numbers for statistical tests
    prox <- nrow(tad.proximal.genes)
    distal <- nrow(df.protein) - prox
    prox.targets <- nrow(tad.proximal.targets)
    prox.targets.activated <- nrow(tad.proximal.targets[which(tad.proximal.targets$Regulation == -1),])
    prox.targets.repressed <- nrow(tad.proximal.targets[which(tad.proximal.targets$Regulation == 1),])
    targets <- nrow(df.protein.targets.filtered)
    targets.activated <- nrow(df.protein.targets.filtered[which(df.protein.targets.filtered$Regulation == -1),])
    targets.repressed <- nrow(df.protein.targets.filtered[which(df.protein.targets.filtered$Regulation == 1),]) 
    
    #hypergeometric tests
    activator <- phyper(prox.targets.activated - 1, prox, distal, targets.activated, lower.tail=F)  
    repressor <- phyper(prox.targets.repressed - 1, prox, distal, targets.repressed, lower.tail=F)  
    
    #Passing results to data.frame
    Results[1,] <- c("TAD", prox, distal, prox.targets.activated, targets.activated, 
                     prox.targets.repressed, targets.repressed, activator, repressor)

  }
  
  #############################
  #          Window           #
  #############################
  windows_vector <- seq(500000, 5000000, 500000)
  for(f in 1:length(windows_vector)){
    half_window <- windows_vector[f]/2
    lncRNA.window <- data.frame(CHR=lncRNA.chr, START=lncRNA.tss-half_window, END=lncRNA.tss+half_window)
    
    proximal.genes <- filter(df.protein, CHR == lncRNA.window$CHR & TSS > lncRNA.window$START & TSS < lncRNA.window$END)
    proximal.targets <- filter(proximal.genes, SYMBOL %in% df.protein.targets.filtered$SYMBOL)
    
    #estimate numbers for statistical tests
    prox <- nrow(proximal.genes)
    distal <- nrow(df.protein) - prox
    prox.targets <- nrow(proximal.targets)
    prox.targets.activated <- nrow(proximal.targets[which(proximal.targets$Regulation == -1),])
    prox.targets.repressed <- nrow(proximal.targets[which(proximal.targets$Regulation == 1),])
    targets <- nrow(df.protein.targets.filtered)
    targets.activated <- nrow(df.protein.targets.filtered[which(df.protein.targets.filtered$Regulation == -1),])
    targets.repressed <- nrow(df.protein.targets.filtered[which(df.protein.targets.filtered$Regulation == 1),]) 
    
    #hypergeometric tests
    activator <- phyper(prox.targets.activated - 1, prox, distal, targets.activated, lower.tail=F)  
    repressor <- phyper(prox.targets.repressed - 1, prox, distal, targets.repressed, lower.tail=F)  
    
    #Passing results to data.frame
    Results[2+f,] <- c(paste0("W",windows_vector[f]), prox, distal, prox.targets.activated, targets.activated, 
                     prox.targets.repressed, targets.repressed, activator, repressor)
    
    
  }
  
  
  #############################
  #         Distance          #
  #############################
  df.protein.chr <- filter(df.protein, CHR == lncRNA.chr)
  df.protein.chr <- data.frame(df.protein.chr, distance = abs(df.protein.chr$TSS - lncRNA.tss))
  df.protein.chr <- df.protein.chr[order(df.protein.chr$distance),]
  
  Activated <- filter(df.protein.chr, Regulation == -1)
  Repressed <- filter(df.protein.chr, Regulation == 1)
  Activated.distance <- mean(Activated$distance)
  Repressed.distance <- mean(Repressed$distance)
  
  x=0
  Activated.distance.random = matrix(nrow = 100, ncol = 1)
  Repressed.distance.random = matrix(nrow = 100, ncol = 1)
  set.seed(23)
  for(x in 1:100){
    df.protein.random <- transform(df.protein.chr, Regulation = sample(Regulation))
    Activated.random <- filter(df.protein.random, Regulation == -1)
    Repressed.random <- filter(df.protein.random, Regulation == 1)
    Activated.distance.random[x,1] <- mean(Activated.random$distance)
    Repressed.distance.random[x,1] <- mean(Repressed.random$distance)  
  }
  
  prox <- nrow(df.protein.chr)
  distal <- nrow(df.protein) - prox
  prox.targets.activated <- nrow(Activated)
  prox.targets.repressed <- nrow(Repressed)
  activator <- round(1-(length(which(Activated.distance.random > Activated.distance))+1)/101,2)
  repressor <- round(1-(length(which(Repressed.distance.random > Repressed.distance))+1)/101,2)
  
  Results[2,] <- c("DistanceMetric", prox, distal, prox.targets.activated, targets.activated,
                   prox.targets.repressed, targets.repressed, activator, repressor)
  
  return(list(
    Results
  ))
}
