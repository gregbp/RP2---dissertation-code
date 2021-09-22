
## INSTALL THIS PACKAGE IF IT IS NOT ALREADY INSTALLED
install.packages("NameNeedle")

library(NameNeedle)



start_time = Sys.time()

print("Program started running at: ")
print(start_time) 

## calculates the longest continuous match between 2 aligned sequences
longest_continuous_match_score <- function(seq1, seq2){
  
  # current continuous match
  curr_match <- 0
  
  # top continuous match
  max_match <- 0
  
  # comparison of the two aligned sequences character by character
  for(i in 1:nchar(seq1)){
    
    # if the compared sequences have the same character, increase curr_match by 1
    if(substr(seq1, i, i) == substr(seq2, i, i)){
      curr_match <- curr_match + 1
    }
    else{ # otherwise check if the curr_match is the top found so far and then set curr_match to 0
      
      if(curr_match > max_match){
        max_match <- curr_match
      }
      
      curr_match <- 0
      
    }
    
    
  }
  
  ## for identical seqs !!!!!
  if(curr_match > max_match){
    max_match <- curr_match
  }
  
  return(max_match)
}





# print("Working directory: ")
# print(getwd()) 

## CHANGE THIS TO YOUR FILESYSTEM
## INPUT FILE  
input <- read.table("30.S02.FilteredDF.AASeq.ByTaxID.FromStaxids.10-12AA_10Sig99.Min100.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")



size <- dim(input)[1]
cmnt <- ""
dataset <- "S02"





# Needles algorithm is going to use these default parameters
# defaultNeedleParams <- list(MATCH=1, MISMATCH=-1, GAP=-1, GAPCHAR="*")

# create a single ID for each peptide
IDs <- list()
for (i in 1:length(input$qseq)) {
  if(dataset == "S01"){
    IDs <- c(IDs, paste0("hg",i,"_",input$qseq[i]))
  }else if(dataset == "S02"){
    IDs <- c(IDs, paste0("m",i,"_",input$qseq[i]))
  }
  
}
input$IDs <-IDs

# Output matrices with nxn dimensions: m1 = Alignment score, m2 = Longest continuous pairwise matches
# and m3 = Network binary matrix

m1 <- matrix(nrow = length(input$qseq), ncol=length(input$qseq), dimnames = list(IDs, IDs))
m2 <- matrix(nrow = length(input$qseq), ncol=length(input$qseq), dimnames = list(IDs, IDs))
m3 <- matrix(nrow = length(input$qseq), ncol=length(input$qseq), dimnames = list(IDs, IDs))


# dataframe that stores the source and target nodes of the network (same info with m3, different representation)
nw.epitopes <- data.frame("source" = character(1), "target" = character(1), "num of consecutive AAs"= numeric(1), stringsAsFactors=FALSE)




# # threshold for filling m3 and df3
thres = 5



for (i in 1:nrow(input)) {
  
  # sequence to be matched
  pattern <- as.character(input$qseq[i])
  
  # cat("\n input$qseq[i] = ",input$qseq[i])
  # cat("\ti = ", i,"\n\n")
  
  for (j in 1:nrow(input)) {
    
    # avoid comparing a seq with itself
    if (j==i){
      break
    }
    # cat("\nj = ", j,"\t input$qseq[j] = ",input$qseq[j],"\n")
    
    
    # sequence to be matched against "pattern"
    subject <- as.character(input$qseq[j])
    
    
    # run Needles algorithm
    nd <- needles(pattern, subject, params=defaultNeedleParams)
    
    # store alignment score in m1
    m1[i,j] <- nd$score
    
    # store the longest continuous match score between the 2 aligned sequences in m2
    m2[i,j] <- longest_continuous_match_score(nd$align1,nd$align2)
    
    # if longest continuous match score between the 2 aligned sequences is greater than threshold,
    # then put 1 in m3
    if (m2[i,j]>thres){
      m3[i,j] <- 1
      
      
      nw.epitopes <- rbind(nw.epitopes, c(as.character(input$IDs[i]), as.character(input$IDs[j]), m2[i,j]))
      
      
      
      # else put 0 in m3
    }else{
      m3[i,j] <- 0
      
    }
    
  }
}


# delete 1st row (it's empty)
nw.epitopes <- nw.epitopes[-c(1), ]

## export the 3 matrices to 3 files
size=dim(input)[1]
write.table(m1, file = paste0("m1_", dataset, "_" , size, ".txt"), append = F, row.names = T, col.names = TRUE, sep = ",", quote = T)
write.table(m2, file = paste0("m2_", dataset, "_" , size, ".txt"), append = F, row.names = T, col.names = TRUE, sep = ",", quote = T)
write.table(m3, file = paste0("m3_", dataset, "_" , size, ".txt"), append = F, row.names = T, col.names = TRUE, sep = ",", quote = T)


## export nw.epitopes to a file
write.table(nw.epitopes, file = paste0("nw",dataset,".epitopes_", size, cmnt,".txt"), append = F, row.names = F, col.names = TRUE, sep = "\t", quote = F)

