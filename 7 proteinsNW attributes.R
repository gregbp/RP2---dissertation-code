library(seqinr)
library(plotly)
library(sets)
library(plyr)



# # ### CREATE nw.species.attrs df
setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged\\output\\checked res")

input_plus <- read.table("S01.TOTAL_with_taxon_levels.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
nw.epitopes <- read.table(file = "nwS01.epitopes_30227.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
nw.proteins <- read.table(file = "nwS01.proteins.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
# nw.proteins <- read.table(file = "nwS01.proteins_queryAccVer.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
dataset <- "S01"


# input_plus <- read.table("S02.TOTAL_with_taxon_levels.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
# nw.epitopes <- read.table(file = "nwS02.epitopes_68476.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
# nw.proteins <- read.table(file = "nwS02.proteins.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
# dataset <- "S02"




## store the NGS counts of the epitopes of epitopesNW
nw.proteins.NGS <- data.frame("node_ID" = character(1), "NGS_counts_log10"= numeric(1),stringsAsFactors=FALSE)

unique.epits <- c(nw.epitopes$source, nw.epitopes$target)
unique.epits <- unique(unique.epits)
head(unique.epits)
tail(unique.epits)
length(unique.epits)

# for (x in unique.epits[1:5]) {
for (x in unique.epits) {
  
  if(length(c(input_plus[input_plus$IDs == x,]$subjectAccVer, log10(input_plus[input_plus$IDs == x,]$NGSCount))) == 2){
    # if(length(c(input_plus[input_plus$IDs == x,]$queryAccVer, log10(input_plus[input_plus$IDs == x,]$NGSCount))) == 2){  
    
    nw.proteins.NGS <- rbind(nw.proteins.NGS, c(input_plus[input_plus$IDs == x,]$subjectAccVer, log10(input_plus[input_plus$IDs == x,]$NGSCount)))
    
    # nw.proteins.NGS <- rbind(nw.proteins.NGS, c(input_plus[input_plus$IDs == x,]$queryAccVer, log10(input_plus[input_plus$IDs == x,]$NGSCount)))
  }
  
}
length(unique.epits)
dim(nw.proteins.NGS)
length(unique(nw.proteins.NGS$node_ID))

nw.proteins.NGS <- nw.proteins.NGS[-c(1), ]

head(nw.proteins.NGS)
tail(nw.proteins.NGS)

# n_occur <-  data.frame(table(nw.proteins.NGS$node_ID))
# n_occur[n_occur$Freq > 1,]
# dim(n_occur[n_occur$Freq > 1,])[1]

# get node attributes for nw.protein

unique.prots <- c(nw.proteins$source, nw.proteins$target)
unique.prots  <- unique(unique.prots )
length(unique.prots )


# nw.proteins.attrs <- data.frame("node_ID" = character(1), "Avg_NGS_counts_log10"= numeric(1), "epitopes_per_protein"=numeric(1), "protein_acc_ver"=numeric(1), stringsAsFactors=FALSE)
nw.proteins.attrs <- data.frame("node_ID" = character(1), "Avg_NGS_counts_log10"= numeric(1), "epitopes_per_protein"=numeric(1), stringsAsFactors=FALSE)
for (x in unique.prots) {
  
  # x = "WP_012314765.1"
  # x = 1031
  
  x <- unique.prots[2]
  # x <-"2A9E_A"
  
  ngsc <- nw.proteins.NGS[nw.proteins.NGS$node_ID == x,]
  ngsc
  avg <- sum(as.numeric(ngsc$NGS_counts_log10))/dim(ngsc)[1]
  avg
  epit.per.prot <- dim(ngsc)[1]
  epit.per.prot
  
  
  # nw.species.attrs <- rbind(nw.species.attrs, c(x, log10(avg), input_plus[input_plus$staxids == x,]$specie[1]))
  
  # nw.proteins.attrs <- rbind(nw.proteins.attrs, c(x, avg, epit.per.prot,input_plus[input_plus$subjectAccVer == x,]$queryAccVer[1]))
  nw.proteins.attrs <- rbind(nw.proteins.attrs, c(x, avg, epit.per.prot))
  
  
}

nw.proteins.attrs <- nw.proteins.attrs[-c(1), ]
head(nw.proteins.attrs)

setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged\\output\\checked res")
write.table(nw.proteins.attrs, file = paste0("nw",dataset,".proteins_nodes_attributes2.txt"), append = F, row.names = F, col.names = TRUE, sep = "\t", quote = F)
# write.table(nw.proteins.attrs, file = paste0("nw",dataset,".proteins_nodes_attributes_queryAccVer.txt"), append = F, row.names = F, col.names = TRUE, sep = "\t", quote = F)


