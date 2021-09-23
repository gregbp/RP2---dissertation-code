library(seqinr)
library(plotly)
library(sets)
library(plyr)

#### SHORT DESCRITPION: this script creates file the file that describes proteinsNW
#### interactions ("proteinsNW_interactions.txt")

## change this to "HG" for the healthy group dataset
dataset <- "DMG"

input_enriched <- read.table(file = paste0(dataset,"_input_enriched.txt"),sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
nw.epitopes.interactions <- read.table(file = paste0(dataset,"_filtered_epitopesNW_interactions.txt"),sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")


 

# create nw.proteins data.frame
nw.proteins <- data.frame("source" = character(1), "target" = character(1), "num of consecutive AAs"= numeric(1), stringsAsFactors=FALSE)
for (row in 1:nrow(nw.epitopes.interactions)) {

  
  src <- nw.epitopes.interactions[row,]$source
  trg <- nw.epitopes.interactions[row,]$target
  numAAs <- nw.epitopes.interactions[row,]$num.of.consecutive.AAs
  
  nw.proteins <- rbind(nw.proteins, c(input_enriched[input_enriched$IDs == src,]$subjectAccVer, input_enriched[input_enriched$IDs == trg,]$subjectAccVer, numAAs))
 
}

# delete 1st empty row (it's empty)
nw.proteins <- nw.proteins[-c(1), ]



# merge duplicate lines of nw.proteins data.frame

# first, sort nw.proteins
nw.proteins <- nw.proteins[order(nw.proteins$source, nw.proteins$target),]


# then merge multiple same interactions - save merged data to nw.proteins2 dataframe
nw.proteins2 <- data.frame("source" = character(1), "target" = character(1), "Avg num of consecutive AAs"= numeric(1), "num of interactions"= numeric(1), stringsAsFactors=FALSE)
 

prvs.src <- ""
prvs.trg <- ""

i <- 1
for (row in 1:nrow(nw.proteins)) {
   
  src <- nw.proteins[row,]$source
  trg <- nw.proteins[row,]$target
   
  
  if(src != prvs.src || trg != prvs.trg){
     
    
    nw.proteins2 <- rbind(nw.proteins2, c(src, trg, nw.proteins[row,]$num.of.consecutive.AAs, 1))
    i <- i + 1
     
    
  }else{
     
    
    nw.proteins2[i,]$Avg.num.of.consecutive.AAs <- as.numeric(nw.proteins2[i,]$Avg.num.of.consecutive.AAs) + as.numeric(nw.proteins[row,]$num.of.consecutive.AAs)
    nw.proteins2[i,]$num.of.interactions <- as.numeric(nw.proteins2[i,]$num.of.interactions) + 1
   
    
  }
  
  
  prvs.src <- nw.proteins[row,]$source
  prvs.trg <- nw.proteins[row,]$target
  
  
  
}

# delete 1st empty row (it's empty)
nw.proteins2 <- nw.proteins2[-c(1), ]


# for each row, divide "Avg.num.of.consecutive.AAs" by "num.of.interactions"
# to get the actual average (now "Avg.num.of.consecutive.AAs" is still sum)
for (row in 1:nrow(nw.proteins2)) {
  nw.proteins2[row,]$Avg.num.of.consecutive.AAs <- as.numeric(nw.proteins2[row,]$Avg.num.of.consecutive.AAs) / as.numeric(nw.proteins2[row,]$num.of.interactions)
}

  
 

# delete "self-loops"
rows.to.be.deleted <- c()
for (row in 1:nrow(nw.proteins2)) {
  if(nw.proteins2[row,]$source == nw.proteins2[row,]$target){
    rows.to.be.deleted <- c(rows.to.be.deleted, row)
    
  }
}
 
# delete 1st empty row (it's empty)
nw.proteins2 <- nw.proteins2[-c(rows.to.be.deleted), ]
 
# export nw.proteins2 to a file 
write.table(nw.proteins2, file = paste0(dataset,"_proteinsNW_interactions.txt"), append = F, row.names = F, col.names = TRUE, sep = "\t", quote = F)
 

