TransCistor <- function(input.file, id.type, species, cell = "H1-NPC",
                        lncRNA.name, lncRNA.chr, lncRNA.tss, lncRNA.strand = "+",
                        enricher.threshold = 0, simulations = 100000 ){
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
  library(pCalibrate)
  
  options(scipen = 999)
  `%!in%` = Negate(`%in%`)
  
  #TEST PARAMS - REMOVE WHEN DONE#######
  
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
  
  df.raw <- unique(clusterProfiler::select(df.raw, -ENTREZ))
  
  #############################
  #      Pre-processing       #
  #############################
  ### removal of same-strand overlapping genes
  # if(any(df.raw$chr == lncRNA.chr & df.raw$TSS == lncRNA.tss)){
  #   lncRNA.dat <- df.raw[df.raw$chr == lncRNA.chr & df.raw$TSS == lncRNA.tss,][1,]
  #   remove.rows <- (df.raw$start %in% lncRNA.dat$start:lncRNA.dat$end |
  #                              df.raw$end %in% lncRNA.dat$start:lncRNA.dat$end) &
  #     df.raw$strand == lncRNA.strand &
  #     df.raw$chr == lncRNA.chr &
  #     !(df.raw$chr == lncRNA.chr & df.raw$TSS == lncRNA.tss)
  #   if(any(remove.rows)){
  #     print('Removing same strand overlaps:')
  #     print(df.raw[remove.rows,])
  #     df.raw <- df.raw[!remove.rows,]
  #   }
  # }
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
                              Proximal.Acticated = rep(0, length(TADs)), Proximal.Repressed = rep(0, length(TADs)),
                              Total = rep(0, length(TADs)), Total.Activated = rep(0, length(TADs)), Total.Repressed = rep(0, length(TADs)), 
                              Activator = rep(1, length(TADs)), Repressor = rep(1, length(TADs)), Activator.mid = rep(1,length(TADs)),
                              Repressor.mid = rep(1,length(TADs)), Activator.Bayes = rep(1,length(TADs)), Repressor.Bayes = rep(1,length(TADs)))
    Results_TADs[,1] = TADs
    windows_vector <- seq(50000, 1000000, 50000)
    Results_window = data.frame(Window = "",Proximal = rep(0, length(windows_vector)), 
                                Proximal.Acticated = rep(0, length(windows_vector)), Proximal.Repressed = rep(0, length(windows_vector)),
                                Total = rep(0, length(windows_vector)), Total.Activated = rep(0, length(windows_vector)), Total.Repressed = rep(0, length(windows_vector)), 
                                Activator = rep(1, length(windows_vector)), Repressor = rep(1, length(windows_vector)))
    Results_window[,1] <- windows_vector
    Results_distance = data.frame(Method="" ,Proximal = 0, 
                                  Proximal.Acticated = 0, Proximal.Repressed = 0,
                                  Total = 0, Total.Activated = 0, Total.Repressed = 0, 
                                  Activator = 1, Repressor = 1)
    Results_distance[,1] <- "Distance"
    Pvalues <- data.frame(TAD.Activator = 1, TAD.Repressor = 1,
                          TAD.Acti.Mid = 1, TAD.Repr.Mid = 1,
                          TAD.Acti.Bayes = 1, TAD.Repr.Bayes = 1,
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
                        TAD.Acti.Mid = 1, TAD.Repr.Mid = 1,
                        TAD.Acti.Bayes = 1, TAD.Repr.Bayes = 1,
                        Window.Activator = 1, Window.Repressor = 1,
                        Distance.Activator = 1, Distance.Repressor = 1)
  #############################
  #            TAD            #
  #############################
  TADs = list.files("References/TADs/", species)
  Results_TADs = data.frame(TAD = "",Proximal = rep(0, length(TADs)), 
                            Proximal.Acticated = rep(0, length(TADs)), Proximal.Repressed = rep(0, length(TADs)),
                            Total = rep(0, length(TADs)), Total.Activated = rep(0, length(TADs)), Total.Repressed = rep(0, length(TADs)), 
                            Activator = rep(1, length(TADs)), Repressor = rep(1, length(TADs)), Activator.mid = rep(1,length(TADs)),
                            Repressor.mid = rep(1,length(TADs)), Activator.Bayes = rep(1,length(TADs)), Repressor.Bayes = rep(1,length(TADs)))
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
      
      if(prox.targets.activated == 0){
        activator <- list(p.value = list(p.fi = 1, p.mid = 1, p.lie = 1))
      }else{
        actiMat <- matrix(c(prox.targets.activated,prox-prox.targets.activated,
                            targets.activated-prox.targets.activated,total - targets.activated - (prox - prox.targets.activated)), nrow = 2, byrow = T)
        activator <- twoby2Calibrate(actiMat,type = 'one.sided', direction = 'greater')
      }
      
      if(prox.targets.repressed == 0){
        repressor <- list(p.value = list(p.fi = 1, p.mid = 1, p.lie = 1))
      }else{
        reprMat <- matrix(c(prox.targets.repressed,prox-prox.targets.repressed,
                            targets.repressed-prox.targets.repressed,total - targets.repressed - (prox - prox.targets.repressed)), nrow = 2, byrow = T)
        repressor <- twoby2Calibrate(reprMat, type = 'one.sided', direction = 'greater')
      }
      #Passing results to data.frame
      Results_TADs[f,2:ncol(Results_TADs)] <- c(prox, prox.targets.activated, prox.targets.repressed,
                                                total, targets.activated, targets.repressed,
                                                activator$p.value[['p.fi']], repressor$p.value[['p.fi']],
                                                activator$p.value[['p.mid']], repressor$p.value[['p.mid']],
                                                activator$p.value[['p.lie']], repressor$p.value[['p.lie']])
      Results_TADs[f,1] <- TADs[f]
    }
    
    #Only form the harmonic mean over results which are not 1 by default (because there is no overlapping TAD)
    OverlapTADs <- Results_TADs[Results_TADs$Proximal != 0,]
    if(nrow(OverlapTADs) > 0){
      Pvalues[1:6] <- c(harmonic.mean(OverlapTADs$Activator),harmonic.mean(OverlapTADs$Repressor),
                        harmonic.mean(OverlapTADs$Activator.mid),harmonic.mean(OverlapTADs$Repressor.mid),
                        harmonic.mean(OverlapTADs$Activator.Bayes),harmonic.mean(OverlapTADs$Repressor.Bayes))
    }else{Pvalues[1:6] <- rep(1,6)}
    TAD_bounaries = data.frame(chr = lncRNA.chr, 
                               start = min(TAD_bounaries$start),
                               end = max(TAD_bounaries$end))
    LocusPlot.genes <- filter(df.raw, chr == TAD_bounaries$chr & TSS > TAD_bounaries$start & TSS < TAD_bounaries$end)
    
  }
  
  
  windows_vector <- seq(50000, 1000000, 50000)
  Results_window = data.frame(Window = "",Proximal = rep(0, length(windows_vector)),
                              Proximal.Acticated = rep(0, length(windows_vector)), Proximal.Repressed = rep(0, length(windows_vector)),
                              Total = rep(0, length(windows_vector)), Total.Activated = rep(0, length(windows_vector)), Total.Repressed = rep(0, length(windows_vector)),
                              Activator = rep(1, length(windows_vector)), Repressor = rep(1, length(windows_vector)))

  targets <- nrow(df.raw.targets.filtered)
  targets.activated <- nrow(df.raw.targets.filtered[which(df.raw.targets.filtered$Regulation == -1),])
  targets.repressed <- nrow(df.raw.targets.filtered[which(df.raw.targets.filtered$Regulation == 1),])
  Pvalues[7:8] <- c(1,1)
  
  #############################
  #         Distance          #
  #############################
  Results_distance = data.frame(Method="" ,Proximal = 0, 
                                Proximal.Acticated = 0, Proximal.Repressed = 0,
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
  if(is.na(Activated.distance) & is.na(Repressed.distance)){
    activator <- 1
    repressor <- 1
  }else{
    df.raw.random <- df.raw.chr
    ### Currently sampling WITHOUT replacement - otherwise the sample needs to be drawn from the original df every loop
    for(x in 1:simulations){
      df.raw.random$Regulation <- sample(df.raw.random$Regulation)
      Activated.distance.random[x,1] <- mean(df.raw.random[df.raw.random$Regulation == -1,'distance'])
      Repressed.distance.random[x,1] <- mean(df.raw.random[df.raw.random$Regulation == 1,'distance'])
    }
    activator <- round((length(which(Activated.distance.random <= Activated.distance))+1)/(simulations+1),7)
    repressor <- round((length(which(Repressed.distance.random <= Repressed.distance))+1)/(simulations+1),7)
  }
  
  if(is.na(Activated.distance)){activator <- 1}
  if(is.na(Repressed.distance)){repressor <- 1}
  
  prox <- nrow(df.raw.chr)
  distal <- nrow(df.raw) - prox
  total <- nrow(df.raw)
  prox.targets.activated <- nrow(Activated)
  prox.targets.repressed <- nrow(Repressed)
  Results_distance[2:ncol(Results_distance)] <- c(prox, prox.targets.activated, prox.targets.repressed,
                                                  total, targets.activated, targets.repressed,
                                                  activator, repressor)
  Results_distance[1] <- "Distance"
  DistancePlot.df = data.frame(MeanDistance = rbind(Activated.distance.random/1000, Repressed.distance.random/1000), 
                               Regulation = c(rep("Activated", simulations), rep("Repressed",simulations)))
  
  Results_distance[2:ncol(Results_distance)] <- c(prox, prox.targets.activated, prox.targets.repressed,
                                                  total, targets.activated, targets.repressed,
                                                  activator, repressor)
  Pvalues[9:10] <- c(Results_distance$Activator, Results_distance$Repressor)
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


