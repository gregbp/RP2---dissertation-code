library(seqinr)
library(plotly)
library(sets)
library(plyr)


#### SHORT DESCRITPION: this script creates file the file that describes speciesNW
#### interactions ("speciesNW_interactions.txt")


# change this to "HG" for the healthy group dataset
dataset <- "DMG"
input_enriched <- read.table(file = paste0(dataset,"_input_enriched.txt"),sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
nw.epitopes.interactions <- read.table(file = paste0(dataset,"_filtered_epitopesNW_interactions.txt"),sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")


 
# epitopes with no specie - delete these epitopes
input_enriched <- input_enriched[input_enriched$specie != "",]

 
# create speciesNW data.frame
nw.species <- data.frame("source" = character(1), "target" = character(1), "num of consecutive AAs"= numeric(1), stringsAsFactors=FALSE)


for (row in 1:nrow(nw.epitopes.interactions)) {
  
  # get the source, target and num of AAs of each row of nw.epitopes
  src <- nw.epitopes.interactions[row,]$source
  trg <- nw.epitopes.interactions[row,]$target
  numAAs <- nw.epitopes.interactions[row,]$num.of.consecutive.AAs
  
   
  # filter deleted epitopes from previously
  if(length(c(input_enriched[input_enriched$IDs == src,]$specie, input_enriched[input_enriched$IDs == trg,]$specie, numAAs)) == 3){
    nw.species <- rbind(nw.species, c(input_enriched[input_enriched$IDs == src,]$specie, input_enriched[input_enriched$IDs == trg,]$specie, numAAs))
  }
}

# delete 1st empty row (it's empty)
nw.species <- nw.species[-c(1), ]

 

# merge duplicate lines of nw.species data.frame

# first, sort nw.species
nw.species <- nw.species[order(nw.species$source, nw.species$target),]

# then merge multiple same interactions - save merged data to nw.species2 dataframe
nw.species2 <- data.frame("source" = character(0), "target" = character(0), "Avg num of consecutive AAs"= numeric(0), "num of interactions"= numeric(0), stringsAsFactors=FALSE)

prvs.src <- ""
prvs.trg <- ""

i <- 0 
for (row in 1:nrow(nw.species)) {
  
  src <- nw.species[row,]$source
  trg <- nw.species[row,]$target
   
  
  ## if src or trg differ from the previous, add a new row to nw.species2
  if(src != prvs.src || trg != prvs.trg){
    tmp.df <- data.frame(src, trg, nw.species[row,]$num.of.consecutive.AAs, 1)
    colnames(tmp.df) <- c("source", "target", "Avg.num.of.consecutive.AAs", "num.of.interactions")
    
    nw.species2 <- rbind(nw.species2, tmp.df)
    i <- i + 1
     
    
  }else{
    nw.species2[i,]$Avg.num.of.consecutive.AAs <- as.numeric(nw.species2[i,]$Avg.num.of.consecutive.AAs) + as.numeric(nw.species[row,]$num.of.consecutive.AAs)
    nw.species2[i,]$num.of.interactions <- as.numeric(nw.species2[i,]$num.of.interactions) + 1
     
  }
  
  
  prvs.src <- nw.species[row,]$source
  prvs.trg <- nw.species[row,]$target
  
  
}
 

# for each row, divide "Avg.num.of.consecutive.AAs" by "num.of.interactions"
# to get the actual average (now "Avg.num.of.consecutive.AAs" is still sum)
for (row in 1:nrow(nw.species2)) {
  nw.species2[row,]$Avg.num.of.consecutive.AAs <- as.numeric(nw.species2[row,]$Avg.num.of.consecutive.AAs) / as.numeric(nw.species2[row,]$num.of.interactions)
}
 


# delete "self-loops"
rows.to.be.deleted <- c()
for (row in 1:nrow(nw.species2)) {
  if(nw.species2[row,]$source == nw.species2[row,]$target){
    rows.to.be.deleted <- c(rows.to.be.deleted, row)
    
  }
}

# delete 1st empty row (it's empty)
nw.species2 <- nw.species2[-c(rows.to.be.deleted), ]

# export nw.species2 to a file 
write.table(nw.species2, file = paste0(dataset,"_speciesNW_interactions.txt"), append = F, row.names = F, col.names = TRUE, sep = "\t", quote = F)



