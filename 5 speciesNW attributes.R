library(seqinr)
library(plotly)
library(sets)
library(plyr)

# ### CREATE nw.species.attrs df
setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged\\output\\checked res")

# input_plus <- read.table("S01.TOTAL_with_taxon_levels.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
# nw.epitopes <- read.table(file = "nwS01.epitopes_30227.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
# nw.species <- read.table(file = "nwS01.species.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
# dataset <- "S01"

input_plus <- read.table("S02.TOTAL_with_taxon_levels.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
nw.epitopes <- read.table(file = "nwS02.epitopes_68476.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
nw.species <- read.table(file = "nwS02.species.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
dataset <- "S02"


# epitopes with no specie - delete these epitopes
length(input_plus[input_plus$specie == "",]$IDs)
length(input_plus[input_plus$specie != "",]$IDs)
dim(input_plus)
input_plus <- input_plus[input_plus$specie != "",]
dim(input_plus)



## store the NGS counts of the epitopes of epitopesNW


unique.epits <- c(nw.epitopes$source, nw.epitopes$target)
unique.epits <- unique(unique.epits)
head(unique.epits)
tail(unique.epits)
length(unique.epits)

nw.species.NGS <- data.frame("node_ID" = character(1), "NGS_counts_log10"= numeric(1), stringsAsFactors=FALSE)

n0_specie <- c() #epitopes from epitopesNW that don't classify as species
# for (x in unique.epits[1:5]) {
for (x in unique.epits) {
  
  ## filter epitopes from epitopesNW that don't classify as species - input_plus has only the "speciesful" epitopes, while unique.epits leads to "specieless" epitopes
  if(length(c(input_plus[input_plus$IDs == x,]$specie, log10(input_plus[input_plus$IDs == x,]$NGSCount))) == 2){
    # print(length(c(input_plus[input_plus$IDs == x,]$staxids, log10(input_plus[input_plus$IDs == x,]$NGSCount))) )
    # print(x)
    nw.species.NGS <- rbind(nw.species.NGS, c(input_plus[input_plus$IDs == x,]$specie, log10(input_plus[input_plus$IDs == x,]$NGSCount)))
  }else{
    n0_specie <- c(n0_specie, x)
  }
  
  
}


# x = "hg465_HYGRKA" 
# length(c(input_plus[input_plus$IDs == x,]$specie, log10(input_plus[input_plus$IDs == x,]$NGSCount)))

nw.species.NGS <- nw.species.NGS[-c(1), ]

length(unique.epits)
length(n0_specie)
dim(nw.species.NGS)
length(unique.epits)-dim(nw.species.NGS)

# n0_specie
# x="hg495_SAYRGG"
# length(c(input_plus[input_plus$IDs == x,]$subjectAccVer, log10(input_plus[input_plus$IDs == x,]$NGSCount)))
# c(input_plus[input_plus$IDs == x,]$subjectAccVer, log10(input_plus[input_plus$IDs == x,]$NGSCount))


head(nw.species.NGS)
tail(nw.species.NGS)


## create nw.species.attrs
unique.sps <- c(nw.species$source, nw.species$target)
unique.sps <- unique(unique.sps)
head(unique.sps)
length(unique.sps)

# nw.species.attrs <- data.frame("node_ID" = character(1), "Avg_NGS_counts_log10"= numeric(1), "specie_name"=character(1), stringsAsFactors=FALSE)
nw.species.attrs <- data.frame("node_ID" = character(1), "Avg_NGS_counts_log10"= numeric(1), "epitopes per specie"= numeric(1), stringsAsFactors=FALSE)
# for (x in unique.sps[1:5]) {
for (x in unique.sps) {
  # x<-"Acidovorax citrulli"
  
  ngsc <- nw.species.NGS[nw.species.NGS$node_ID == x,]
  ngsc
  avg <- sum(as.numeric(ngsc$NGS_counts_log10))/dim(ngsc)[1]
  avg
  epit.per.spec <- dim(ngsc)[1]
  epit.per.spec 
  
  
  # nw.species.attrs <- rbind(nw.species.attrs, c(x, avg, input_plus[input_plus$subjectAccVer == x,]$specie[1]))
  nw.species.attrs <- rbind(nw.species.attrs, c(x, avg, epit.per.spec))
  
}


nw.species.attrs <- nw.species.attrs[-c(1), ]

length(unique.sps)
dim(nw.species.attrs)[1]

head(nw.species.attrs)
tail(nw.species.attrs)

setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged\\output\\checked res")
write.table(nw.species.attrs, file = paste0("nw",dataset,".species_nodes_attributes.txt"), append = F, row.names = F, col.names = TRUE, sep = "\t", quote = F)

