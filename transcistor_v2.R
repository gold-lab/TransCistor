TransCistor <- function(input.file, id.type, species, cell = "H1-NPC",
                        lncRNA.name, lncRNA.chr, lncRNA.tss, lncRNA.strand = "+",
                        enricher.threshold = 0, simulations = 1000 ){

  library(rlang)
  library(dplyr)
  library(rlist)
  library(tidyverse)
  library(clusterProfiler)
  library(ggplot2)
  library(ggpubr)
  library(pheatmap)
  library(RColorBrewer)
  library(psych)
  options(scipen = 999)
  `%!in%` = Negate(`%in%`)
  
  
  #############################
  #   Input/Reference files   #
  #############################
  
  #reading gencode and tad files
  species = tolower(species)
  if((species == "human") || (species == "mouse")){
    gencode_file <- paste0("References/Gencode/", species, "_gencode.txt")
    gencode <- read.table(gencode_file, sep="\t", header=T, stringsAsFactors=FALSE)
    tad_file <- paste0("References/TADs/", species,".ESC.txt")
    if(file.exists(tad_file)){
      tad <- read.table(paste0("References/TADs/", species,".ESC.txt"), sep="\t", header=F, stringsAsFactors=FALSE)
    }else{
      print("No tad file could be read. species should be either human or mouse.")
    }
  }else{
    print("species should be either human or mouse.")
    return()
  }
  
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
  df.raw <- unique(select(df.raw, -ENTREZ))
  
  #############################
  #      Pre-processing       #
  #############################
  #remove non protein coding genes
  #df.protein <- filter(df.raw, type == "protein_coding")
  df.raw <- filter(df.raw, !(chr == lncRNA.chr & TSS == lncRNA.tss))
  df.raw$SYMBOL <- toupper(df.raw$SYMBOL)
  df.raw.targets <- filter(df.raw, Regulation != 0)
  

  if(nrow(df.raw.targets) == 0){
    print(paste0("No targets found. exiting"))
    combined_heatmap <- qplot(1:2,1:2)+xlab("This is a dummy plot")+ylab("No targets found")
    barplot_tad <- qplot(1:2,1:2)+xlab("This is a dummy plot")+ylab("No targets found")
    barplot_tad <- qplot(1:2,1:2)+xlab("This is a dummy plot")+ylab("No targets found")
    TADs <- list.files("References/TADs/", species)
    Results_TADs = data.frame(TAD = "",Proximal = rep(0, length(TADs)), 
                              Proximal.Activated = rep(0, length(TADs)), Proximal.Repressed = rep(0, length(TADs)),
                              Total = rep(0, length(TADs)), Total.Activated = rep(0, length(TADs)), Total.Repressed = rep(0, length(TADs)), 
                              Activator = rep(1, length(TADs)), Repressor = rep(1, length(TADs)))
    Results_TADs[,1] = TADs
    windows_vector <- seq(50000, 1000000, 50000)
    Results_window = data.frame(Window = "",Proximal = rep(0, length(windows_vector)), 
                                Proximal.Activated = rep(0, length(windows_vector)), Proximal.Repressed = rep(0, length(windows_vector)),
                                Total = rep(0, length(windows_vector)), Total.Activated = rep(0, length(windows_vector)), Total.Repressed = rep(0, length(windows_vector)), 
                                Activator = rep(1, length(windows_vector)), Repressor = rep(1, length(windows_vector)))
    Results_window[,1] <- windows_vector
    Results_distance = data.frame(Method="" ,Proximal = 0, 
                                  Proximal.Activated = 0, Proximal.Repressed = 0,
                                  Total = 0, Total.Activated = 0, Total.Repressed = 0, 
                                  Activator = 0, Repressor = 0)
    Results_distance[,1] <- "Distance"
    Pvalues <- data.frame(TAD.Activator = 1, TAD.Repressor = 1,
                          Window.Activator = 1, Window.Repressor = 1,
                          Distance.Activator = 1, Distance.Repressor = 1)
    return(list(
      Results_TADs,
      Results_window,
      Results_distance,
      Pvalues))
  }
  
  
  df.raw.targets.filtered <- df.raw.targets
  Pvalues <- data.frame(TAD.Activator = 1, TAD.Repressor = 1,
                        Window.Activator = 1, Window.Repressor = 1,
                        Distance.Activator = 1, Distance.Repressor = 1)
  #############################
  #            TAD            #
  #############################
  TADs = list.files("References/TADs/", species)
  Results_TADs = data.frame(TAD = "",Proximal = rep(0, length(TADs)), 
            Proximal.Activated = rep(0, length(TADs)), Proximal.Repressed = rep(0, length(TADs)),
            Total = rep(0, length(TADs)), Total.Activated = rep(0, length(TADs)), Total.Repressed = rep(0, length(TADs)), 
            Activator = rep(1, length(TADs)), Repressor = rep(1, length(TADs)))

  TAD_bounaries <- c()
  for(f in 1:length(TADs)){
    tad <- read.table(paste0("References/TADs/",TADs[f]))
    colnames(tad) <- c("chr", "start", "end")
    lncRNA.tad <- filter(tad,chr == lncRNA.chr & start <= lncRNA.tss & end >= lncRNA.tss)
    if(nrow(lncRNA.tad) == 0){
      print("No ovalapping TAD found. Pvalues of 1 reported for TAD window.")
      tad.proximal.targets = c()
      LocusPlot.genes = c()
    }else{
      TAD_bounaries <- rbind(TAD_bounaries, lncRNA.tad)
      tad.proximal.genes <- filter(df.raw, chr == lncRNA.tad$chr & TSS > lncRNA.tad$start & TSS < lncRNA.tad$end)
      tad.proximal.targets <- filter(tad.proximal.genes, SYMBOL %in% df.raw.targets.filtered$SYMBOL)
      
      #estimate numbers for statistical tests
      prox <- nrow(tad.proximal.genes)
      distal <- nrow(df.raw) - prox
      total <- nrow(df.raw)
      prox.targets <- nrow(tad.proximal.targets)
      prox.targets.activated <- nrow(tad.proximal.targets[which(tad.proximal.targets$Regulation == -1),])
      prox.targets.repressed <- nrow(tad.proximal.targets[which(tad.proximal.targets$Regulation == 1),])
      targets <- nrow(df.raw.targets.filtered)
      targets.activated <- nrow(df.raw.targets.filtered[which(df.raw.targets.filtered$Regulation == -1),])
      targets.repressed <- nrow(df.raw.targets.filtered[which(df.raw.targets.filtered$Regulation == 1),]) 
      
      #hypergeometric tests
      activator <- phyper(prox.targets.activated - 1, prox, distal, targets.activated, lower.tail=F)  
      repressor <- phyper(prox.targets.repressed - 1, prox, distal, targets.repressed, lower.tail=F)  
      
      #Passing results to data.frame
      Results_TADs[f,2:ncol(Results_TADs)] <- c(prox, prox.targets.activated, prox.targets.repressed,
                                                total, targets.activated, targets.repressed,
                                                activator, repressor)
      Results_TADs[f,1] <- TADs[f]
      
    }
    #gather pvalues
    Pvalues[1:2] <- c(harmonic.mean(Results_TADs$Activator),harmonic.mean(Results_TADs$Repressor))
    TAD_bounaries = data.frame(chr = lncRNA.chr, 
                                            start = min(TAD_bounaries$start),
                                            end = max(TAD_bounaries$end))
    LocusPlot.genes <- filter(df.raw, chr == TAD_bounaries$chr & TSS > TAD_bounaries$start & TSS < TAD_bounaries$end)
    
  }

  #############################
  #          Window           #
  #############################
  windows_vector <- seq(50000, 1000000, 50000)
  Results_window = data.frame(Window = "",Proximal = rep(0, length(windows_vector)), 
                            Proximal.Activated = rep(0, length(windows_vector)), Proximal.Repressed = rep(0, length(windows_vector)),
                            Total = rep(0, length(windows_vector)), Total.Activated = rep(0, length(windows_vector)), Total.Repressed = rep(0, length(windows_vector)), 
                            Activator = rep(1, length(windows_vector)), Repressor = rep(1, length(windows_vector)))
  
  for(f in 1:length(windows_vector)){
    half_window <- windows_vector[f]/2
    lncRNA.window <- data.frame(chr=lncRNA.chr, start=lncRNA.tss-half_window, end=lncRNA.tss+half_window)
    
    proximal.genes <- filter(df.raw, chr == lncRNA.window$chr & TSS > lncRNA.window$start & TSS < lncRNA.window$end)
    proximal.targets <- filter(proximal.genes, SYMBOL %in% df.raw.targets.filtered$SYMBOL)
    
    #estimate numbers for statistical tests
    prox <- nrow(proximal.genes)
    distal <- nrow(df.raw) - prox
    total <- nrow(df.raw)
    prox.targets <- nrow(proximal.targets)
    prox.targets.activated <- nrow(proximal.targets[which(proximal.targets$Regulation == -1),])
    prox.targets.repressed <- nrow(proximal.targets[which(proximal.targets$Regulation == 1),])
    targets <- nrow(df.raw.targets.filtered)
    targets.activated <- nrow(df.raw.targets.filtered[which(df.raw.targets.filtered$Regulation == -1),])
    targets.repressed <- nrow(df.raw.targets.filtered[which(df.raw.targets.filtered$Regulation == 1),]) 
    
    #hypergeometric tests
    activator <- phyper(prox.targets.activated - 1, prox, distal, targets.activated, lower.tail=F)  
    repressor <- phyper(prox.targets.repressed - 1, prox, distal, targets.repressed, lower.tail=F)  
    
    #Passing results to data.frame
    Results_window[f,] <- c(windows_vector[f], prox, prox.targets.activated, prox.targets.repressed,
                          total, targets.activated, targets.repressed,
                          activator, repressor)
  }
  Pvalues[3:4] <- c(harmonic.mean(Results_window$Activator),harmonic.mean(Results_window$Repressor))
  
  #############################
  #         Distance          #
  #############################
  Results_distance = data.frame(Method="" ,Proximal = 0, 
                              Proximal.Activated = 0, Proximal.Repressed = 0,
                              Total = 0, Total.Activated = 0, Total.Repressed = 0, 
                              Activator = 0, Repressor = 0)
  
  df.raw.chr <- filter(df.raw, chr == lncRNA.chr)
  df.raw.chr <- data.frame(df.raw.chr, distance = abs(df.raw.chr$TSS - lncRNA.tss))
  df.raw.chr <- df.raw.chr[order(df.raw.chr$distance),]
  
  Activated <- filter(df.raw.chr, Regulation == -1)
  Repressed <- filter(df.raw.chr, Regulation == 1)
  Activated.distance <- mean(Activated$distance)
  Repressed.distance <- mean(Repressed$distance)
  
  x=0
  Activated.distance.random = matrix(nrow = simulations, ncol = 1)
  Repressed.distance.random = matrix(nrow = simulations, ncol = 1)
  set.seed(23)
  for(x in 1:simulations){
    df.raw.random <- transform(df.raw.chr, Regulation = sample(Regulation))
    Activated.random <- filter(df.raw.random, Regulation == -1)
    Repressed.random <- filter(df.raw.random, Regulation == 1)
    Activated.distance.random[x,1] <- mean(Activated.random$distance)
    Repressed.distance.random[x,1] <- mean(Repressed.random$distance)  
  }
  
  prox <- nrow(df.raw.chr)
  distal <- nrow(df.raw) - prox
  total <- nrow(df.raw)
  prox.targets.activated <- nrow(Activated)
  prox.targets.repressed <- nrow(Repressed)
  activator <- round(1-(length(which(Activated.distance.random > Activated.distance))+1)/(simulations+1),7)
  repressor <- round(1-(length(which(Repressed.distance.random > Repressed.distance))+1)/(simulations+1),7)
  
  DistancePlot.df = data.frame(MeanDistance = rbind(Activated.distance.random/1000, Repressed.distance.random/1000), 
                               Regulation = c(rep("Activated", simulations), rep("Repressed",simulations)))
  
  Results_distance[2:ncol(Results_distance)] <- c(prox, prox.targets.activated, prox.targets.repressed,
                          total, targets.activated, targets.repressed,
                          activator, repressor)
  Results_distance[1] <- "Distance"
  Pvalues[5:6] <- c(Results_distance$Activator, Results_distance$Repressor)
    return(list(
      Results_TADs,
      Results_window,
      Results_distance,
      Pvalues,
      DistancePlot.df,
      Activated.distance,
      Repressed.distance,
      LocusPlot.genes
      
  ))
}
