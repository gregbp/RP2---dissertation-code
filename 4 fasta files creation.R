library(seqinr)
library(plotly)
library(plyr) 

#### SHORT DESCRITPION: this script converts subnetwork csv files to fasta 
#### in order for the fasta to be used in CLC (create cladograms)

## change this to "HG" for the healthy group dataset
dataset <- "DMG"

for(sb in 1:10){
  
  # read the input files (.csv)
  if (dataset == "HG"){
    input <- read.table(paste0("HG_subNW", sb, ".csv"),sep=",",header=TRUE,fill=T, comment.char = "")
  }else if (dataset == "DMG"){
    input <- read.table(paste0("DMG_subNW", sb, ".csv"),sep=",",header=TRUE,fill=T, comment.char = "")
  }  
  
  nms <- c()
  seqs <- c()
  
  # convert csv to fasta format
  for(i in 1:nrow(input)){
    
    input[i,]$name
    
    nm <- input[i,]$name
    seq <- strsplit(input[i,]$name,"_")[[1]][[2]]
    
    
    nms <- c(nms, nm)
    seqs <- c(seqs, seq)
  }
  
  # write output file (.fasta)
  if (dataset == "HG"){
    write.fasta(sequences = as.list(seqs), names = nms, file.out = paste0("HG_subNW", sb, ".fasta"))
  }else if (dataset == "DMG"){
    write.fasta(sequences = as.list(seqs), names = nms, file.out = paste0("DMG_subNW", sb, ".fasta"))
  } 
}




