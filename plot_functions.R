library(tidyr)
library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)
library(ggrepel)
library(VennDiagram)
library(stringr)
library(RColorBrewer)
library(reshape)
library(ggh4x)
library(scales)
plot_digital <- function(Data, title, Pval_act, Pval_rep){
  library(ggplot2, lib.loc = "libraries/")
  TADplot.df.Act <- cbind.data.frame(Regulation = rep("Activated", nrow(Data)), 
                                     Location = rep(c("Proximal","Proximal","Distal","Distal"), nrow(Data)), 
                                     Target = rep(c("Targets","Non-Targets","Targets","Non-Targets"), nrow(Data)), 
                                     COUNT = 0)
  TADplot.df.Rep = cbind.data.frame(Regulation = rep("Repressed", nrow(Data)), 
                                    Location = rep(c("Proximal","Proximal","Distal","Distal"), nrow(Data)), 
                                    Target = rep(c("Targets","Non-Targets","Targets","Non-Targets"), nrow(Data)))
  
  pal <- c(
    "Activated" = "#1B9E77",
    "Repressed" = "#D95F02",
    "ns" = "#666666"
  )
  pal = pal[c(rep(3,2),1:2)]
  names(pal) = c("Activated.Non-Targets","Repressed.Non-Targets",
                 "Activated.Targets","Repressed.Targets")
  
  X = c()
  for(f in 1:nrow(Data)){
    i = ((f-1)*4)+1:4
    X = rbind(X,i)
    TADplot.df.Act[i,"COUNT"] = c(Data$Proximal.Acticated[f], Data$Proximal[f] - Data$Proximal.Acticated[f],
                                  Data$Total.Activated[f] - Data$Proximal.Acticated[f], (Data$Total[f]-Data$Proximal[f])-(Data$Total.Activated[f] - Data$Proximal.Acticated[f]))
    
    TADplot.df.Rep[i,"COUNT"] = c(Data$Proximal.Repressed[f], Data$Proximal[f] - Data$Proximal.Repressed[f],
                                  Data$Total.Repressed[f] - Data$Proximal.Repressed[f], (Data$Total[f]-Data$Proximal[f])-(Data$Total.Repressed[f] - Data$Proximal.Repressed[f]))
    
    
  }
  TADplot.df <- rbind(cbind.data.frame(Regulation = "Activated", Location = c("Proximal","Proximal","Distal","Distal"), Target = c("Targets","Non-Targets","Targets","Non-Targets"), 
                                       COUNT = c(median(TADplot.df.Act[X[,1],"COUNT"]), median(TADplot.df.Act[X[,2],"COUNT"]), 
                                                 median(TADplot.df.Act[X[,3],"COUNT"]), median(TADplot.df.Act[X[,4],"COUNT"]))), 
                      cbind.data.frame(Regulation = "Repressed", Location = c("Proximal","Proximal","Distal","Distal"), Target = c("Targets","Non-Targets","Targets","Non-Targets"), 
                                       COUNT = c(median(TADplot.df.Rep[X[,1],"COUNT"]), median(TADplot.df.Rep[X[,2],"COUNT"]), 
                                                 median(TADplot.df.Rep[X[,3],"COUNT"]), median(TADplot.df.Rep[X[,4],"COUNT"]))))
  
  Regulation.labs <- paste0(paste0("Activator pvalue = ", round(Pval_act, digits = 3)), "  ---  ", paste0("Repressor pvalue = ", round(Pval_rep, digits = 3)))
  digital.plot <- ggplot(TADplot.df, aes(x=Location, y=COUNT, fill=interaction(Regulation,Target), alpha=Target)) + 
    geom_bar(stat="identity", position="fill") + 
    facet_wrap(facets = "Regulation") +
    scale_fill_manual(values = pal, limits = names(pal)) + 
    geom_text(aes(label = COUNT), position = position_fill(vjust = 0.5), alpha=1, size=5) + 
    theme_bw() +  
    theme(axis.text.x = element_text(hjust = 0.5, size=18, angle = 90)) + 
    theme(axis.text.y = element_text(hjust = 0.5, size=18)) +
    theme(axis.title.x = element_text(size=18)) + 
    theme(axis.title.y = element_text(size=18)) + 
    theme(title=element_text(size=15,hjust = 0.5)) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(strip.text.x = element_text(size = 18)) +
    labs(x = "Relative location", y = "Fraction of genes", alpha = "") + guides(fill="none") +
    ggtitle(paste0(title,": [",Regulation.labs,"]"))
  return(digital.plot)
  
}


plot_analogue <- function(Data, title, Activated.distance, Repressed.distance, Pval_act, Pval_rep){
  
  library(ggplot2, lib.loc = "libraries/")
  library(ggpubr, lib.loc = "libraries/")
  Regulation.labs <- paste0(paste0("Activator pvalue = ", round(Pval_act, digits = 3)), "  ---  ", paste0("Repressor pvalue = ", round(Pval_rep, digits = 3)))
  analogue.plot <- gghistogram(Data, x = "MeanDistance", 
                               rug = TRUE,
                               color = "Regulation", 
                               bins=50, 
                               palette = c("#1b9e77", "#d95f02"),
                               facet.by = "Regulation") + 
    geom_vline(data=filter(Data, Regulation=="Activated"), 
               aes(xintercept=Activated.distance/1000), colour="#1b9e77", lwd=2) +
    geom_vline(data=filter(Data, Regulation=="Repressed"), 
               aes(xintercept=Repressed.distance/1000), colour="#d95f02", lwd=2) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Mean Distance (kb)") +
    ylab("Random count") +     ggtitle(paste0(title,": [",Regulation.labs,"]")) +
    theme(axis.text.x = element_text(hjust = 0.5, size=18, angle = 90)) + 
    theme(axis.text.y = element_text(hjust = 0.5, size=18)) +
    theme(axis.title.x = element_text(size=18)) + 
    theme(axis.title.y = element_text(size=18)) + 
    theme(title=element_text(size=15,hjust = 0.5)) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(strip.text.x = element_text(size = 18))  
  
  return(analogue.plot)
  
}

plot_locus <- function(Data){
  library(ggplot2)
  library(ggrepel)
  Data$Regulation[which(Data$Regulation == 1)] <- "Repressed"
  Data$Regulation[which(Data$Regulation == -1)] <- "Activated"
  Data$Regulation[which(Data$Regulation == 0)] <- "ns"
  
  pal <- c(
    "lncRNA" = "#E78AC3",
    "Activated" = "#1B9E77",
    "Repressed" = "#D95F02",
    "ns" = "#666666"
  )
  
  #For plotting umlilo example only -remove later
  Data$RectScale <- 1
  Data[Data$Regulation == "lncRNA","RectScale"] <- 2
  Data[Data$Regulation == "Activated","RectScale"] <- 2
  
  locus_plot <- ggplot(Data, aes(TSS)) +
    geom_hline(yintercept = 0, col="gray") +
    # geom_text(aes(max(Data$TSS) - ((max(Data$TSS)-min(Data$TSS))/2), -1, 
    #               label = paste0(unique(Data$chr),":",min(Data$TSS),"-",max(Data$TSS)), vjust = -1)) +
    geom_rect(aes(fill = Regulation, xmin = TSS-5000, xmax = TSS+5000, ymin=-0.2 * RectScale, ymax=0.2 * RectScale), alpha=0.8) +
 
    scale_fill_manual(values = pal) + theme_void() 
  
  
  return(locus_plot)
}

