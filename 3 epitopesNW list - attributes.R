library(seqinr)
library(plotly)
library(sets)
library(plyr)

#### SHORT DESCRITPION: this script creates epitopesNW attributes file ("filtered_epitopesNW_list.txt")


# change this to "HG" for the healthy group dataset
dataset <- "DMG"
  
input_enriched <- read.table(file = paste0(dataset,"_input_enriched.txt"),sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
nw.epitopes.interactions <- read.table(file = paste0(dataset,"_filtered_epitopesNW_interactions.txt"),sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")



# keep the the distinct IDs of the filtered epitopes
unique.IDs <- c(nw.epitopes.interactions$source, nw.epitopes.interactions$target)
unique.IDs <- unique(unique.IDs)
 

# create a data.frame that stores the distinct filtered epitopes and their attributes
nw.epitopes_nodes_attrs <- data.frame("node_ID"=character(1), "qseq"=character(1), "NGS_counts_log10"=numeric(1), "queryAccVer"=numeric(1), "subjectAccVer"=character(1),
                                      "org_taxID"=numeric(1), "org_name"=character(1), "specie"=character(1), "genus"=character(1),
                                      "family"=character(1), "phylum"=character(1), "superkingdom"=character(1),stringsAsFactors=FALSE)

for (x in unique.IDs) {
  
  
  node_ID <- x
  qseq <- input_enriched[input_enriched$IDs == x,]$qseq
  NGS_counts_log10 <- log10(input_enriched[input_enriched$IDs == x,]$NGSCount)
  queryAccVer <- input_enriched[input_enriched$IDs == x,]$queryAccVer
  subjectAccVer <- input_enriched[input_enriched$IDs == x,]$subjectAccVer
  org_taxID <- input_enriched[input_enriched$IDs == x,]$staxids
  org_name <- input_enriched[input_enriched$IDs == x,]$Organism
  
  specie <- input_enriched[input_enriched$IDs == x,]$specie
  genus <- input_enriched[input_enriched$IDs == x,]$genus
  family <- input_enriched[input_enriched$IDs == x,]$family
  phylum <- input_enriched[input_enriched$IDs == x,]$phylum
  superkingdom <- input_enriched[input_enriched$IDs == x,]$superkingdom
  
  nw.epitopes_nodes_attrs <- rbind(nw.epitopes_nodes_attrs, c(node_ID, qseq, NGS_counts_log10, queryAccVer, subjectAccVer, org_taxID, org_name, specie, genus, family, phylum, superkingdom))
  
  
}

# delete 1st empty row (it's empty)
nw.epitopes_nodes_attrs <- nw.epitopes_nodes_attrs[-c(1), ]


# export nw.epitopes_nodes_attrs to a file
write.table(nw.epitopes_nodes_attrs, file = paste0(dataset,"_filtered_epitopesNW_list.txt"), append = F, row.names = F, col.names = TRUE, sep = "\t", quote = F)



