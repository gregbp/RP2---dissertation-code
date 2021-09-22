library(seqinr)
library(plotly)
library(sets)
library(plyr)

# # #### CREATE THE 2 FILES DESCRIBING proteins NW

setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged\\output\\checked res")

input_plus <- read.table("S01.TOTAL_with_taxon_levels.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
dataset <- "S01"
nw.epitopes <- read.table(file = "nwS01.epitopes_30227.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")

# input_plus <- read.table("S02.TOTAL_with_taxon_levels.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
# dataset <- "S02"
# nw.epitopes <- read.table(file = "nwS02.epitopes_68476.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
# 



nw.proteins <- data.frame("source" = character(1), "target" = character(1), "num of consecutive AAs"= numeric(1), stringsAsFactors=FALSE)
for (row in 1:nrow(nw.epitopes)) {
  # for (row in 1:5) {
  
  src <- nw.epitopes[row,]$source
  trg <- nw.epitopes[row,]$target
  numAAs <- nw.epitopes[row,]$num.of.consecutive.AAs
  
  src
  trg
  numAAs
  
  # tmp.df <- data.frame(input[input$IDs == src,]$subjectAccVer, input[input$IDs == trg,]$subjectAccVer, numAAs)
  # colnames(tmp.df) <- c("source", "target","num.of.consecutive.AAs")
  
  nw.proteins <- rbind(nw.proteins, c(input_plus[input_plus$IDs == src,]$subjectAccVer, input_plus[input_plus$IDs == trg,]$subjectAccVer, numAAs))
  # nw.proteins <- rbind(nw.proteins, c(input_plus[input_plus$IDs == src,]$queryAccVer, input_plus[input_plus$IDs == trg,]$queryAccVer, numAAs))
  
}

nw.proteins
nw.proteins <- nw.proteins[-c(1), ]
nw.proteins
cat("\ninitital nw.proteins size: ",  dim(nw.proteins)[1])
cat("\nnw,epitopes size: ", dim(nw.epitopes)[1])

head(nw.proteins)
tail(nw.proteins)

#### queryAccVer OR subjectAccVer
# length(unique(input_plus$queryAccVer))
# length(input_plus$queryAccVer)
# length(unique(input_plus$subjectAccVer))
# length(input_plus$subjectAccVer)



# sort nw.protein
nw.proteinsTMP <- nw.proteins
nw.proteins <- nw.proteins[order(nw.proteins$source, nw.proteins$target),]
head(nw.proteins)
tail(nw.proteins)

dim(nw.proteins)
dim(nw.epitopes)


# merging multiple same interactions
# nw.proteins2 <- data.frame("source" = character(0), "target" = character(0), "Avg num of consecutive AAs"= numeric(0), "num of interactions"= numeric(0), stringsAsFactors=FALSE)
nw.proteins2 <- data.frame("source" = character(1), "target" = character(1), "Avg num of consecutive AAs"= numeric(1), "num of interactions"= numeric(1), stringsAsFactors=FALSE)
nw.proteins2


prvs.src <- ""
prvs.trg <- ""

i <- 1
for (row in 1:nrow(nw.proteins)) {
  # row <-50
  # nw.proteins[48:50,]
  src <- nw.proteins[row,]$source
  trg <- nw.proteins[row,]$target
  
  src
  trg
  
  if(src != prvs.src || trg != prvs.trg){
    
    # tmp.df <- data.frame(src, trg, nw.proteins[row,]$num.of.consecutive.AAs, 1)
    # colnames(tmp.df) <- c("source", "target", "Avg.num.of.consecutive.AAs", "num.of.interactions")
    # nw.proteins2 <- rbind(nw.proteins2, tmp.df)
    
    nw.proteins2 <- rbind(nw.proteins2, c(src, trg, nw.proteins[row,]$num.of.consecutive.AAs, 1))
    i <- i + 1
    # print("eleos1")
    
  }else{
    # print(row)
    
    nw.proteins2[i,]$Avg.num.of.consecutive.AAs <- as.numeric(nw.proteins2[i,]$Avg.num.of.consecutive.AAs) + as.numeric(nw.proteins[row,]$num.of.consecutive.AAs)
    # nw.proteins2[i,]$Avg.num.of.consecutive.AAs <- nw.proteins2[i,]$Avg.num.of.consecutive.AAs + nw.proteins[row,]$num.of.consecutive.AAs
    
    nw.proteins2[i,]$num.of.interactions <- as.numeric(nw.proteins2[i,]$num.of.interactions) + 1
    # print("eleos2")
    
    
  }
  
  
  prvs.src <- nw.proteins[row,]$source
  prvs.trg <- nw.proteins[row,]$target
  
  
  
}


nw.proteins2 <- nw.proteins2[-c(1), ]



head(nw.proteins2,10)
tail(nw.proteins2,10)

head(nw.proteins,20)
nw.proteins2[48:52,]
nw.proteins[48:52,]


# for each row, divide "Avg.num.of.consecutive.AAs" by "num.of.interactions"
# to get the actual average (now "Avg.num.of.consecutive.AAs" is still sum)
for (row in 1:nrow(nw.proteins2)) {
  
  nw.proteins2[row,]$Avg.num.of.consecutive.AAs <- as.numeric(nw.proteins2[row,]$Avg.num.of.consecutive.AAs) / as.numeric(nw.proteins2[row,]$num.of.interactions)
}

head(nw.proteins2,30)
tail(nw.proteins2,10)
nw.proteins2[48:52,]
nw.proteins[48:52,]



nw.proteins23 <- nw.proteins2

cat("\nnw.proteins size after merging same interactions: ",  dim(nw.proteins23)[1])
cat("\nnw,epitopes size: ", dim(nw.epitopes)[1])
cat("\ndifference: ",dim(nw.epitopes)[1]-dim(nw.proteins23)[1])

# delete "self-loops"
rows.to.be.deleted <- c()
for (row in 1:nrow(nw.proteins2)) {
  if(nw.proteins2[row,]$source == nw.proteins2[row,]$target){
    rows.to.be.deleted <- c(rows.to.be.deleted, row)
    
  }
}
rows.to.be.deleted
length(rows.to.be.deleted)
## check
dim(nw.proteins2[nw.proteins2$source == nw.proteins2$target,])

cat("\nself loops deleted: ",  length(rows.to.be.deleted))
cat("\nnw.proteins size after self loops deletion: ",  dim(nw.proteins2)[1])
cat("\nnw,epitopes size: ", dim(nw.epitopes)[1])
cat("\ndifference between inital nw.epitopes & final nw.species: ",dim(nw.epitopes)[1]-dim(nw.proteins2)[1])
cat("\nmerging + self loops deleted:", dim(nw.epitopes)[1]-dim(nw.proteins23)[1] + length(rows.to.be.deleted))


nw.proteins2 <- nw.proteins2[-c(rows.to.be.deleted), ]
dim(nw.proteins23)
dim(nw.proteins2)

dim(nw.proteins2[nw.proteins2$num.of.interactions == 2,])[1]*2 +length(rows.to.be.deleted)
dim(nw.epitopes)[1]-dim(nw.proteins2)[1]

# nw.pr <- read.table("output\\nwS01.proteins2.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
# dim(nw.pr[nw.pr$num.of.interactions == 2,])
# dim(nw.epitopes)-dim(nw.pr)


setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged\\output\\checked res")
write.table(nw.proteins2, file = paste0("nw",dataset,".proteins.txt"), append = F, row.names = F, col.names = TRUE, sep = "\t", quote = F)
# write.table(nw.proteins2, file = paste0("nw",dataset,".proteins_queryAccVer.txt"), append = F, row.names = F, col.names = TRUE, sep = "\t", quote = F)



