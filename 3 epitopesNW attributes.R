library(seqinr)
library(plotly)
library(sets)
library(plyr)

## CREATE epitopes NW attributes file


setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged")


# input_plus <- read.table("S01.TOTAL_with_taxon_levels.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
# dataset <- "S01"
# # nw.epitopes <- read.table(file = "nwS01.epitopes_merged.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
# nw.epitopes <- read.table(file = "nwS01.epitopes_30227.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")

input_plus <- read.table("S02.TOTAL_with_taxon_levels.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
dataset <- "S02"
# nw.epitopes <- read.table(file = "nwS02.epitopes_merged.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
nw.epitopes <- read.table(file = "nwS02.epitopes_68476.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")

## create nw.species.attrs df
unique.IDs <- c(nw.epitopes$source, nw.epitopes$target)
unique.IDs <- unique(unique.IDs)
head(unique.IDs)
tail(unique.IDs)
length(unique.IDs)


nw.epitopes_nodes_attrs <- data.frame("node_ID"=character(1), "qseq"=character(1), "NGS_counts_log10"=numeric(1), "queryAccVer"=numeric(1), "subjectAccVer"=character(1),
                                      "org_taxID"=numeric(1), "org_name"=character(1), "specie"=character(1), "genus"=character(1),
                                      "family"=character(1), "phylum"=character(1), "superkingdom"=character(1),stringsAsFactors=FALSE)

# for (x in unique.IDs[1:5]){ 
for (x in unique.IDs) {
  
  
  node_ID <- x
  qseq <- input_plus[input_plus$IDs == x,]$qseq
  NGS_counts_log10 <- log10(input_plus[input_plus$IDs == x,]$NGSCount)
  queryAccVer <- input_plus[input_plus$IDs == x,]$queryAccVer
  subjectAccVer <- input_plus[input_plus$IDs == x,]$subjectAccVer
  org_taxID <- input_plus[input_plus$IDs == x,]$staxids
  org_name <- input_plus[input_plus$IDs == x,]$Organism
  
  specie <- input_plus[input_plus$IDs == x,]$specie
  genus <- input_plus[input_plus$IDs == x,]$genus
  family <- input_plus[input_plus$IDs == x,]$family
  phylum <- input_plus[input_plus$IDs == x,]$phylum
  superkingdom <- input_plus[input_plus$IDs == x,]$superkingdom
  
  nw.epitopes_nodes_attrs <- rbind(nw.epitopes_nodes_attrs, c(node_ID, qseq, NGS_counts_log10, queryAccVer, subjectAccVer, org_taxID, org_name, specie, genus, family, phylum, superkingdom))
  
  
}
nw.epitopes_nodes_attrs <- nw.epitopes_nodes_attrs[-c(1), ]
head(nw.epitopes_nodes_attrs)
tail(nw.epitopes_nodes_attrs)



setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged\\output")
write.table(nw.epitopes_nodes_attrs, file = paste0("nw",dataset,".epitopes_nodes_attributes.txt"), append = F, row.names = F, col.names = TRUE, sep = "\t", quote = F)