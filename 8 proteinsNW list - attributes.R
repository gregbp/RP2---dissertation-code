library(seqinr)
library(plotly)
library(sets)
library(plyr)


#### SHORT DESCRITPION: this script creates proteinsNW attributes file ("proteinsNW_list.txt")


## change this to "HG" for the healthy group dataset
dataset <- "DMG"

input_enriched <- read.table(file = paste0(dataset,"_input_enriched.txt"),sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
nw.epitopes.interactions <- read.table(file = paste0(dataset,"_filtered_epitopesNW_interactions.txt"),sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
nw.proteins.interactions <- read.table(file = paste0(dataset,"_proteinsNW_interactions.txt"),sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")


 


## store the NGS counts of the epitopes of epitopesNW
nw.proteins.NGS <- data.frame("node_ID" = character(1), "NGS_counts_log10"= numeric(1),stringsAsFactors=FALSE)

unique.epits <- c(nw.epitopes.interactions$source, nw.epitopes.interactions$target)
unique.epits <- unique(unique.epits)



for (x in unique.epits) {
  
  if(length(c(input_enriched[input_enriched$IDs == x,]$subjectAccVer, log10(input_enriched[input_enriched$IDs == x,]$NGSCount))) == 2){

    nw.proteins.NGS <- rbind(nw.proteins.NGS, c(input_enriched[input_enriched$IDs == x,]$subjectAccVer, log10(input_enriched[input_enriched$IDs == x,]$NGSCount)))
 
  }
  
}

# delete 1st empty row (it's empty)
nw.proteins.NGS <- nw.proteins.NGS[-c(1), ]


# create nw.proteins.attrsdataframe that stores proteins attributes
unique.prots <- c(nw.proteins.interactions$source, nw.proteins.interactions$target)
unique.prots  <- unique(unique.prots )
 
nw.proteins.attrs <- data.frame("node_ID" = character(1), "Avg_NGS_counts_log10"= numeric(1), "epitopes_per_protein"=numeric(1), stringsAsFactors=FALSE)
for (x in unique.prots) {
  
    
  ngsc <- nw.proteins.NGS[nw.proteins.NGS$node_ID == x,]
  avg <- sum(as.numeric(ngsc$NGS_counts_log10))/dim(ngsc)[1]
  epit.per.prot <- dim(ngsc)[1]
  
   
  nw.proteins.attrs <- rbind(nw.proteins.attrs, c(x, avg, epit.per.prot))
  
  
}

# delete 1st empty row (it's empty)
nw.proteins.attrs <- nw.proteins.attrs[-c(1), ]
 

## export nw.proteins.attrs to a file
write.table(nw.proteins.attrs, file = paste0(dataset,"_proteinsNW_list.txt"), append = F, row.names = F, col.names = TRUE, sep = "\t", quote = F)




