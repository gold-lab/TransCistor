source( "plot_functions.R")
source( "transcistor.R")

#####Usage####
TADs = c()
Distance = c()
Pvalues = c()

### Input the Identifier for the lncRNA####
input='lincIRX5'

file_name <- metadata[metadata$SYMBOL==input,]$File_Name
num=1 #### File index number. Can be changes for lncRNAs with multiple experiments###

i=rownames(metadata[metadata$SYMBOL==input,])
print(paste0("current_file : ", file_name[num]))
current_file <- read.table(paste0(data_path, file_name[num]), header= FALSE, sep = '\t', comment.char = '')


  
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

### Store all the results####
TADs <- rbind(TADs, results[[1]])
Distance <- rbind(Distance, results[[3]])
Pvalues <- rbind(Pvalues, results[[4]])
Data = results[[1]]
Data = as.data.frame(Data)
Data = Data[which(Data$Total != 0),]
title = name
Pval_act = Pvalues$TAD.Acti.Mid
Pval_rep = Pvalues$TAD.Repr.Mid

####Plot TransCistor-digital results####

plot_digital(Data, title, Pval_act, Pval_rep)

Data2 = results[[5]]
Data2 = as.data.frame(Data2)
Act.dist <- results[[6]]
Rep.dist <- results[[7]]
Pval_act2 = Pvalues$Distance.Activator
Pval_rep2 = Pvalues$Distance.Repressor

####Plot TransCistor-analogue results####
plot_analogue(Data2, title, Act.dist, Rep.dist, Pval_act2, Pval_rep2)

#### Plot lncRNA chromosome position plot s####

Data3 = as.data.frame(results[[8]])
print(head(Data3))
Data3 = data.frame(ENSEMBL = c(Data3$ENSEMBL,name),
                   TSS = c(Data3$TSS, tss),
                   Regulation = c(Data3$Regulation, "lncRNA"),
                   chr = c(Data3$chr, chr))
plot_locus(Data3)




