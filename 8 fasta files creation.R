library(seqinr)
library(plotly)
library(plyr) 


#### CREATE FASTA FILES FOR CLC 

dataset <- "S02"
# sb<-1
for(sb in 1:10){
  # read the input file
  
  if (dataset == "S01"){
    setwd("C:\\Users\\grigo\\Desktop\\total S01 NWs\\subNWs_S01")
    input <- read.table(paste0("nwS01.epitopes_30227_subNW", sb, ".csv"),sep=",",header=TRUE,fill=T, comment.char = "")
    # head(input)
  }else if (dataset == "S02"){
    setwd("C:\\Users\\grigo\\Desktop\\total S02 NWs\\subNWs_S02")
    input <- read.table(paste0("nwS02.epitopes_68476_subNW", sb, ".csv"),sep=",",header=TRUE,fill=T, comment.char = "")
    # head(input)
  }  
  
  nms <- c()
  seqs <- c()
  
  # for(i in 1:5){
  for(i in 1:nrow(input)){
    # i=1
    input[i,]$name
    
    nm <- input[i,]$name
    seq <- strsplit(input[i,]$name,"_")[[1]][[2]]
    
    nm
    seq
    
    nms <- c(nms, nm)
    seqs <- c(seqs, seq)
  }
  nms
  seqs 
  
  
  
  if (dataset == "S01"){
    setwd("C:\\Users\\grigo\\Desktop\\total S01 NWs")
    write.fasta(sequences = as.list(seqs), names = nms, file.out = paste0("fasta_files_subNWs_S01\\nwS01.epitopes_sub", sb, ".fasta"))
  }else if (dataset == "S02"){
    setwd("C:\\Users\\grigo\\Desktop\\total S02 NWs")
    write.fasta(sequences = as.list(seqs), names = nms, file.out = paste0("fasta_files_subNWs_S02\\nwS02.epitopes_sub", sb, ".fasta"))
  } 
}




