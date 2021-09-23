
#### SHORT DESCRITPION: this script adds taxonomic columns to the original datasets



 
## uncomment this to use "HG" dataset for the healthy group dataset

# input <- read.table("30.S01.FilteredDF.AASeq.ByTaxID.FromStaxids.10-12AA.10Sig99.Min100.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
# dataset <- "HG"
# tax_groups <- read.table("HG_lineage_TaxonToolkit.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")


input <- read.table("30.S02.FilteredDF.AASeq.ByTaxID.FromStaxids.10-12AA_10Sig99.Min100.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
dataset <- "DMG"
tax_groups <- read.table("DMG_lineage_TaxonToolkit.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
head(tax_groups)



# create a single ID for each epitope
IDs <- c()
for (i in 1:length(input$qseq)) {
  if(dataset == "HG"){
    IDs <- c(IDs, paste0("hg",i,"_",input$qseq[i]))
  }else if(dataset == "DMG"){
    IDs <- c(IDs, paste0("dmg",i,"_",input$qseq[i]))
  }
  
}
input$IDs <-IDs


# add the new "taxonomic" columns
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

input$specie <- specie
input$genus <- genus
input$family <- family
input$phylum <- phylum
input$superkingdom <- superkingdom
 

# create input_enriched.txt file
write.table(input, file = paste0(dataset,"_input_enriched.txt"), append = F, row.names = F, col.names = TRUE, sep = "\t", quote = F)

