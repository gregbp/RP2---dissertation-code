library(seqinr)
library(plotly)
library(sets)
library(plyr)


# # #### CREATE THE 2 FILES DESCRIBING species NW

setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged\\output\\checked res")
# 
# input_plus <- read.table("S01.TOTAL_with_taxon_levels.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
# dataset <- "S01"
# nw.epitopes <- read.table(file = "nwS01.epitopes_30227.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")

input_plus <- read.table("S02.TOTAL_with_taxon_levels.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
dataset <- "S02"
nw.epitopes <- read.table(file = "nwS02.epitopes_68476.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")


# epitopes with no specie - delete these epitopes
length(input_plus[input_plus$specie == "",]$IDs)
length(input_plus[input_plus$specie != "",]$IDs)
dim(input_plus)
input_plus <- input_plus[input_plus$specie != "",]
dim(input_plus)
# setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged\\output")
# write.table(input_plus, file = paste0("nwS01.epitopes(only completed specie rows).txt"), append = F, row.names = F, col.names = TRUE, sep = "\t", quote = F)
length(unique(input_plus$staxids))
length(unique(input_plus$specie))


# # create speciesNW df
nw.species <- data.frame("source" = character(1), "target" = character(1), "num of consecutive AAs"= numeric(1), stringsAsFactors=FALSE)



# # "translate" peptidesNW df to speciesNW df
# for (row in 1:5) {
for (row in 1:nrow(nw.epitopes)) {
  
  
  
  # row=325
  
  ## get the source, target and num of AAs of each row of nw.epitopes
  src <- nw.epitopes[row,]$source
  trg <- nw.epitopes[row,]$target
  numAAs <- nw.epitopes[row,]$num.of.consecutive.AAs
  
  src
  trg
  numAAs
  
  # filter deleted epitopes from previously
  if(length(c(input_plus[input_plus$IDs == src,]$specie, input_plus[input_plus$IDs == trg,]$specie, numAAs)) == 3){
    nw.species <- rbind(nw.species, c(input_plus[input_plus$IDs == src,]$specie, input_plus[input_plus$IDs == trg,]$specie, numAAs))
  }
}

nw.species <- nw.species[-c(1), ]

head(nw.species)
tail(nw.species)

nw.speciesTMP <- nw.species
nw.speciesTMP

cat("\ninitital nw.species size: ",  dim(nw.speciesTMP)[1])
cat("\nnw,epitopes size: ", dim(nw.epitopes)[1])

## part 2 - merge duplicate lines of speciesNW df

# sort nw.species

head(nw.speciesTMP)
tail(nw.speciesTMP)
nw.species <- nw.species[order(nw.species$source, nw.species$target),]
head(nw.species)
tail(nw.species)


# merging multiple same interactions
nw.species2 <- data.frame("source" = character(0), "target" = character(0), "Avg num of consecutive AAs"= numeric(0), "num of interactions"= numeric(0), stringsAsFactors=FALSE)


prvs.src <- ""
prvs.trg <- ""

i <- 0
# for (row in 1:5) {
for (row in 1:nrow(nw.species)) {
  # row=1
  src <- nw.species[row,]$source
  trg <- nw.species[row,]$target
  # print(src)
  #     print(trg)
  
  ## if src or trg differ from the previous, add a new row to nw.species2
  if(src != prvs.src || trg != prvs.trg){
    tmp.df <- data.frame(src, trg, nw.species[row,]$num.of.consecutive.AAs, 1)
    colnames(tmp.df) <- c("source", "target", "Avg.num.of.consecutive.AAs", "num.of.interactions")
    
    nw.species2 <- rbind(nw.species2, tmp.df)
    i <- i + 1
    # print("eleos1")
    
  }else{
    nw.species2[i,]$Avg.num.of.consecutive.AAs <- as.numeric(nw.species2[i,]$Avg.num.of.consecutive.AAs) + as.numeric(nw.species[row,]$num.of.consecutive.AAs)
    nw.species2[i,]$num.of.interactions <- as.numeric(nw.species2[i,]$num.of.interactions) + 1
    # print("eleos2")
  }
  
  # print(i)
  # cat("\n")
  prvs.src <- nw.species[row,]$source
  prvs.trg <- nw.species[row,]$target
  # prvs.src
  # prvs.trg
  
}
cat("\nnw.species size after merging same interactions: ",  dim(nw.species2 )[1])
cat("\nnw,epitopes size: ", dim(nw.epitopes)[1])
cat("\ndifference: ",dim(nw.epitopes)[1]-dim(nw.species2 )[1])


# head(nw.species2,5)
# #
# head(nw.species,5)

nw.speciesTMP[nw.speciesTMP$source==101510 & nw.speciesTMP$target==1027,]
nw.speciesTMP[nw.speciesTMP$source== 101510   & nw.speciesTMP$target==106590,]
head(nw.species2,10)
tail(nw.species2,10)
#
# head(nw.species,20)

## for each row, divide "Avg.num.of.consecutive.AAs" by "num.of.interactions"
## to get the actual average (now "Avg.num.of.consecutive.AAs" is still sum)
for (row in 1:nrow(nw.species2)) {
  nw.species2[row,]$Avg.num.of.consecutive.AAs <- as.numeric(nw.species2[row,]$Avg.num.of.consecutive.AAs) / as.numeric(nw.species2[row,]$num.of.interactions)
}
# head(nw.species2,10)
# tail(nw.species2,10)
# nw.speciesTMP[nw.speciesTMP$source==986  & nw.speciesTMP$target==80866,]
# nw.species2[nw.species2$source==986  & nw.species2$target==80866,]
#
# nw.species2[nw.species2$source==101510  & nw.species2$target==106590 ,]

nw.species23 <- nw.species2



# delete "self-loops"
rows.to.be.deleted <- c()
for (row in 1:nrow(nw.species2)) {
  if(nw.species2[row,]$source == nw.species2[row,]$target){
    rows.to.be.deleted <- c(rows.to.be.deleted, row)
    
  }
}
rows.to.be.deleted
length(rows.to.be.deleted)
dim(nw.species2[nw.species2$source == nw.species2$target,])
nw.species2TMP <- nw.species2

dim(nw.species2)[1]
nw.species2 <- nw.species2[-c(rows.to.be.deleted), ]
dim(nw.species23)
dim(nw.species2)[1]


cat("\nself loops deleted: ",  length(rows.to.be.deleted))
cat("\nnw.species size after self loops deletion: ",  dim(nw.species2)[1])
cat("\nnw,epitopes size: ", dim(nw.epitopes)[1])
cat("\ndifference between inital nw.epitopes & final nw.species: ",dim(nw.epitopes)[1]-dim(nw.species2)[1])
cat("\nmerging + self loops deleted:", dim(nw.epitopes)[1]-dim(nw.species23)[1] + length(rows.to.be.deleted))



setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged\\output")
write.table(nw.species2, file = paste0("nw",dataset,".species.txt"), append = F, row.names = F, col.names = TRUE, sep = "\t", quote = F)



