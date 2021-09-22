
## ADD taxonomic columns to the original datasets
setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged")


# input <- read.table("30.S01.FilteredDF.AASeq.ByTaxID.FromStaxids.10-12AA.10Sig99.Min100.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
# dataset <- "S01"
# org.IDs.NCBI <- read.table("S01_tax_report_NCBI.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
# head(org.IDs.NCBI)
# tax_groups <- read.table("S01_lineage_TaxonToolkit.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
# head(tax_groups)





input <- read.table("30.S02.FilteredDF.AASeq.ByTaxID.FromStaxids.10-12AA_10Sig99.Min100.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
dataset <- "S02"
org.IDs.NCBI <- read.table("S02_tax_report_NCBI.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
head(org.IDs.NCBI)
tax_groups <- read.table("S02_lineage_TaxonToolkit.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
head(tax_groups)

##create a single ID for each peptide
IDs <- c()
for (i in 1:length(input$qseq)) {
  if(dataset == "S01"){
    IDs <- c(IDs, paste0("hg",i,"_",input$qseq[i]))
  }else if(dataset == "S02"){
    IDs <- c(IDs, paste0("m",i,"_",input$qseq[i]))
  }
  
}
input$IDs <-IDs


## correct organisms names 
for(i in 1:nrow(input)){
  org.name <- org.IDs.NCBI[org.IDs.NCBI$taxid == input[i,]$staxids ,]$taxname
  input[i,]$Organism <- org.name
  
}
# head(input$staxids)
# head(input$Organism)
# tail(input$staxids)
# tail(input$Organism)


# ## add the new "taxonomic" columns
specie <- c()
genus <- c()
family <- c()
phylum <- c()
superkingdom <- c()

for (i in 1:nrow(input)) {
  
  staxID <- input[i,]$staxids
  
  specie <- c(specie, tax_groups[tax_groups$taxIDs == staxID,]$species)
  genus <- c(genus, tax_groups[tax_groups$taxIDs == staxID,]$genus)
  family <- c(family, tax_groups[tax_groups$taxIDs == staxID,]$family)
  phylum <- c(phylum, tax_groups[tax_groups$taxIDs == staxID,]$phylum)
  superkingdom <- c(superkingdom, tax_groups[tax_groups$taxIDs == staxID,]$superkingdom)
  
  
}
# head(input$staxids,4)
# specie[1:4]
# genus[1:4]
# family[1:4]
# phylum[1:4]
# superkingdom[1:4]

input$specie <- specie
input$genus <- genus
input$family <- family
input$phylum <- phylum
input$superkingdom <- superkingdom

head(input)
tail(input)

setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged\\output")
# write.table(input, file = paste0(dataset,".TOTAL_with_taxon_levels.txt"), append = F, row.names = F, col.names = TRUE, sep = "\t", quote = F)

