library(NameNeedle)


#### SHORT DESCRITPION: this script filters the "input_enriched.txt" file in order to find the epitopes that have
#### more than 5 identical consecutive amino acids with other epitopes 
#### It creates the file that describes filtered epitopesNW interactions ("filtered_epitopesNW_interactions.txt")
#### and the 3 matrices files.



# this function calculates the longest continuous match between 2 aligned sequences
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




## MAIN PROGRAM
start_time = Sys.time()

print("Program started running at: ")
print(start_time) 


## change this to "HG" for the healthy group dataset
dataset <- "DMG"
input_enriched <- read.table(file = paste0(dataset,"_input_enriched.txt"),sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")


 

# output matrices with nxn dimensions: m1 = Alignment score, m2 = Longest continuous pairwise matches and m3 = Network binary matrix
m1 <- matrix(nrow = length(input_enriched$qseq), ncol=length(input_enriched$qseq), dimnames = list(IDs, IDs))
m2 <- matrix(nrow = length(input_enriched$qseq), ncol=length(input_enriched$qseq), dimnames = list(IDs, IDs))
m3 <- matrix(nrow = length(input_enriched$qseq), ncol=length(input_enriched$qseq), dimnames = list(IDs, IDs))


# dataframe that stores the source and target nodes of the network (same info with m3, different representation)
nw.epitopes.interactions <- data.frame("source" = character(1), "target" = character(1), "num of consecutive AAs"= numeric(1), stringsAsFactors=FALSE)




# threshold for filling m3 and nw.epitopes.interactions -> more than 5 identical consecutive amino acids
thres = 5


# compare each epitope sequence with every other sequence
for (i in 1:nrow(input_enriched)) {
  
  # sequence to be matched
  pattern <- as.character(input_enriched$qseq[i])
  
  
  for (j in 1:nrow(input_enriched)) {
    
    # avoid comparing a seq with itself
    if (j==i){
      break
    }
    
    # sequence to be matched against "pattern"
    subject <- as.character(input_enriched$qseq[j])
    
    # run Needles algorithm
    nd <- needles(pattern, subject, params=defaultNeedleParams)
    
    # store alignment score in m1
    m1[i,j] <- nd$score
    
    # store the longest continuous match score between the 2 aligned sequences in m2
    m2[i,j] <- longest_continuous_match_score(nd$align1,nd$align2)
    
    # if longest continuous match score between the 2 aligned sequences is greater than threshold, then put 1 in m3
    if (m2[i,j]>thres){
      m3[i,j] <- 1
      
      nw.epitopes.interactions <- rbind(nw.epitopes.interactions, c(as.character(input_enriched$IDs[i]), as.character(input_enriched$IDs[j]), m2[i,j]))
      
      # else put 0 in m3
    }else{
      m3[i,j] <- 0
      
    }
    
  }
}


# delete 1st row (it's empty)
nw.epitopes.interactions <- nw.epitopes.interactions[-c(1), ]

## export the 3 matrices to 3 files
write.table(m1, file = paste0("m1_", dataset, ".txt"), append = F, row.names = T, col.names = TRUE, sep = ",", quote = T)
write.table(m2, file = paste0("m2_", dataset, ".txt"), append = F, row.names = T, col.names = TRUE, sep = ",", quote = T)
write.table(m3, file = paste0("m3_", dataset, ".txt"), append = F, row.names = T, col.names = TRUE, sep = ",", quote = T)


## export nw.epitopes.interactions to a file
write.table(nw.epitopes.interactions, file = paste0(dataset,"_filtered_epitopesNW_interactions.txt"), append = F, row.names = F, col.names = TRUE, sep = "\t", quote = F)

