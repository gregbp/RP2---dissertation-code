library(seqinr)
library(plotly)
library(sets)
library(plyr)

#### SHORT DESCRITPION: this script creates speciesNW attributes file ("speciesNW_list.txt")


## change this to "HG" for the healthy group dataset
dataset <- "DMG"

input_enriched <- read.table(file = paste0(dataset,"_input_enriched.txt"),sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
nw.epitopes.interactions <- read.table(file = paste0(dataset,"_filtered_epitopesNW_interactions.txt"),sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
nw.species.interactions <- read.table(file = paste0(dataset,"_speciesNW_interactions.txt"),sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")


 
# epitopes with no specie - delete these epitopes
input_enriched <- input_enriched[input_enriched$specie != "",]


## store the NGS counts of the epitopes of epitopesNW
unique.epits <- c(nw.epitopes.interactions$source, nw.epitopes.interactions$target)
unique.epits <- unique(unique.epits)

nw.species.NGS <- data.frame("node_ID" = character(1), "NGS_counts_log10"= numeric(1), stringsAsFactors=FALSE)

n0_specie <- c() #epitopes from epitopesNW that don't classify as species


for (x in unique.epits) {
  
  # filter epitopes from epitopesNW that don't classify as species - input_plus has only the "speciesful" epitopes, while unique.epits leads to "specieless" epitopes
  if(length(c(input_enriched[input_enriched$IDs == x,]$specie, log10(input_enriched[input_enriched$IDs == x,]$NGSCount))) == 2){
    nw.species.NGS <- rbind(nw.species.NGS, c(input_enriched[input_enriched$IDs == x,]$specie, log10(input_enriched[input_enriched$IDs == x,]$NGSCount)))
  }else{
    n0_specie <- c(n0_specie, x)
  }
  
  
}

# delete 1st empty row (it's empty) 
nw.species.NGS <- nw.species.NGS[-c(1), ]

 

# create nw.species.attrs dataframe that stores species attributes
unique.sps <- c(nw.species.interactions$source, nw.species.interactions$target)
unique.sps <- unique(unique.sps)
 
 
nw.species.attrs <- data.frame("node_ID" = character(1), "Avg_NGS_counts_log10"= numeric(1), "epitopes per specie"= numeric(1), stringsAsFactors=FALSE)
 
for (x in unique.sps) {
  
  ngsc <- nw.species.NGS[nw.species.NGS$node_ID == x,]
  avg <- sum(as.numeric(ngsc$NGS_counts_log10))/dim(ngsc)[1]
  epit.per.spec <- dim(ngsc)[1]
  
  nw.species.attrs <- rbind(nw.species.attrs, c(x, avg, epit.per.spec))
  
}

# delete 1st empty row (it's empty)
nw.species.attrs <- nw.species.attrs[-c(1), ]

## export nw.species.attrs to a file
write.table(nw.species.attrs, file = paste0(dataset,"_speciesNW_list.txt"), append = F, row.names = F, col.names = TRUE, sep = "\t", quote = F)
