library(VennDiagram)     
library(seqinr)
library(plotly)
library(sets)
library(plyr) 
library(RColorBrewer)
 


#### SHORT DESCRITPION: this script creates: 
#### 1) Venn diagrams
#### 2) the great table that describes HG and DMG original and filtered datasets
#### 3) the files for creating the distributions in original and filtered epitopes datasets
#### 4) the files for creating comparative distributions original DMG VS HG and filtered DMG VS HG

 



### VENN DIAGRAMS

library(VennDiagram)   
library(RColorBrewer)

datasets <-c("original","filtered")

for(dataset in datasets){
  
  if(dataset == "original"){
    input1 <- read.table(file = "HG_input_enriched.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    input2 <- read.table(file = "DMG_input_enriched.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    
  }else if(dataset == "filtered"){
    input1 <- read.table(file = "HG_filtered_epitopesNW_list.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    input2 <- read.table(file = "DMG_filtered_epitopesNW_list.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    
    
  }
  
   
  ## VennDiagram 1 - epitopes
  HG.epitopes <- input1$qseq
  DMG.epitopes <- input2$qseq
  
  
  
  myCol <- brewer.pal(3, "Pastel2")
  
  venn.diagram(
    x = list(HG.epitopes, DMG.epitopes),
    category.names = c(paste0(dataset," HG - epitopes") , paste0(dataset, " DMG - epitopes")),
    filename = paste0(dataset,'_1_epitopes_venn_diagramm.png'),
    output=TRUE,
    imagetype="png" ,
    resolution = 500,
    compression = "lzw",
    height = 2500 ,
    width = 2500,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = c("green","red"),
    
    cex = 1.3,
    fontfamily = "sans",
    cat.cex = 1.3,
    cat.default.pos = "outer",
    cat.pos = c(-00.1, 00.1)
    
    
  )
  
  
  ## VennDiagram 2 - proteins
  HG.proteins <- input1$subjectAccVer
  DMG.proteins <- input2$subjectAccVer
  
  venn.diagram(
    x = list(HG.proteins, DMG.proteins),
    category.names = c(paste0(dataset," HG - proteins") , paste0(dataset, " DMG - proteins")),
    filename = paste0(dataset,'_2_proteins_venn_diagramm.png'),
    output=TRUE,
    
    imagetype="png" ,
    resolution = 500,
    compression = "lzw",
    height = 2500 ,
    width = 2500,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = c("green","red"),
    
    cex = 1.3,
    fontfamily = "sans",
    cat.cex = 1.3,
    cat.default.pos = "outer",
    cat.pos = c(-00.1, 00.1)
    
  ) 
  
  
  ## VennDiagram 3 - species
  HG.species <- input1$specie
  DMG.species <- input2$specie
  
  
  venn.diagram(
    x = list(HG.species, DMG.species),
    category.names = c(paste0(dataset," HG - species") , paste0(dataset, " DMG - species")),
    filename = paste0(dataset,'_3_species_venn_diagramm.png'),
    output=TRUE,
    imagetype="png" ,
    resolution = 500,
    compression = "lzw",
    height = 2500 ,
    width = 2500,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = c("green","red"),
    
    cex = 1.3,
    fontfamily = "sans",
    cat.cex = 1.3,
    cat.default.pos = "outer",
    cat.pos = c(-20, 20),
    cat.dist = c(0.05, 0.05)
    
    
  ) 
  
  ## VennDiagram 4 - genera
  HG.genera <- input1$genus
  DMG.genera <- input2$genus
  
  
  venn.diagram(
    x = list(HG.genera, DMG.genera),
    category.names = c(paste0(dataset," HG - genera") , paste0(dataset, " DMG - genera")),
    filename = paste0(dataset,'_4_genera_venn_diagramm.png'),
    output=TRUE,
    imagetype="png" ,
    resolution = 500,
    compression = "lzw",
    height = 2500 ,
    width = 2500,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = c("green","red"),
    
    cex = 1.3,
    fontfamily = "sans",
    cat.cex = 1.3,
    cat.default.pos = "outer",
    cat.pos = c(-20, 20),
    cat.dist = c(0.05, 0.05)
    
    
  ) 
  
  
  ## VennDiagram 5 - families
  HG.families <- input1$family
  DMG.families <- input2$family
  
  
  venn.diagram(
    x = list(HG.families, DMG.families),
    category.names = c(paste0(dataset," HG - families") , paste0(dataset, " DMG - families")),
    filename = paste0(dataset,'_5_families_venn_diagramm.png'),
    output=TRUE,
    imagetype="png" ,
    resolution = 500,
    compression = "lzw",
    height = 2500 ,
    width = 2500,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = c("green","red"),
    
    cex = 1.3,
    fontfamily = "sans",
    cat.cex = 1.3,
    cat.default.pos = "outer",
    cat.pos = c(-20, 20),
    cat.dist = c(0.05, 0.05)
    
    
  ) 
  
  ## VennDiagram 6 - phyla
  HG.phyla <- input1$phylum
  DMG.phyla <- input2$phylum
  
  venn.diagram(
    x = list(HG.phyla, DMG.phyla),
    category.names = c(paste0(dataset," HG - phyla") , paste0(dataset, " DMG - phyla")),
    filename = paste0(dataset,'_6_phyla_venn_diagramm tmp.png'),
    output=TRUE,
    imagetype="png" ,
    resolution = 500,
    compression = "lzw",
    height = 2500 ,
    width = 2500,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = c("green","red"),
    
    cex = 1.3,
    fontfamily = "sans",
    cat.cex = 1.3,
    cat.default.pos = "outer",
    cat.pos = c(-20, 20),
    cat.dist = c(0.05, 0.05)
    
    
  ) 
  
  
  ## VennDiagram 7 - superkingdoms
  HG.superkingdoms <- input1$superkingdom
  DMG.superkingdoms <- input2$superkingdom
  
  
  venn.diagram(
    x = list(HG.superkingdoms, DMG.superkingdoms),
    category.names = c(paste0(dataset," HG - superkingdoms") , paste0(dataset, " DMG - superkingdoms")),
    filename = paste0(dataset,'_7_superkingdoms_venn_diagramm.png'),
    output=TRUE,
    imagetype="png" ,
    resolution = 500,
    compression = "lzw",
    height = 2500 ,
    width = 2500,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = c("green","red"),
    
    cex = 1.15,
    fontfamily = "sans",
    cat.cex = 1.15,
    cat.default.pos = "outer",
    cat.pos = c(-27, 27),
    cat.dist = c(0.05, 0.05)
    
    
  ) 
  
}  


 
 

  


### GREAT TABLE that describes HG and DMG original and filtered datasets

descr.df <- data.frame("dataset.size"=numeric(1),
                       "classified.as.specie"=numeric(1), "unclassified.as.specie"=numeric(1),
                       "classified.as.specie.perc"=numeric(1), "unclassified.as.specie.perc"=numeric(1), "distinct.species"=numeric(1),
                       
                       "classified.as.genus"=numeric(1), "unclassified.as.genus"=numeric(1),
                       "classified.as.genus.perc"=numeric(1), "unclassified.as.genus.perc"=numeric(1), "distinct.genera"=numeric(1),
                       
                       "classified.as.family"=numeric(1), "unclassified.as.family"=numeric(1),
                       "classified.as.family.perc"=numeric(1), "unclassified.as.family.perc"=numeric(1), "distinct.families"=numeric(1),
                       
                       "classified.as.phylum"=numeric(1), "unclassified.as.phylum"=numeric(1),
                       "classified.as.phylum.perc"=numeric(1), "unclassified.as.phylum.perc"=numeric(1), "distinct.phyla"=numeric(1),
                       
                       "classified.as.s.kingdom"=numeric(1), "unclassified.as.s.kingdom"=numeric(1),
                       "classified.as.s.kingdom.perc"=numeric(1), "unclassified.as.s.kingdom.perc"=numeric(1), "distinct.s.kingdoms"=numeric(1),
                       
                       stringsAsFactors=FALSE)



datasets <-c("HG.original","DMG.original", "HG.filtered","DMG.filtered")


for(dataset in datasets){
   
  if(dataset == "HG.original"){
    input <- read.table(file = "HG_input_enriched.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    
  }else if(dataset == "DMG.original"){
    input <- read.table(file = "DMG_input_enriched.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    
  }else if(dataset == "HG.filtered"){
    input <- read.table(file = "HG_filtered_epitopesNW_list.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    
  }else if(dataset == "DMG.filtered"){
    input <- read.table(file = "DMG_filtered_epitopesNW_list.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    
  }
  
  
  head(input)
  size <- dim(input)[1]
  print(size)
  
  ## species
  species.classfd <- dim(input[input$specie != "",])[1]
  species.unclassfd <- dim(input[input$specie == "",])[1]
  species.classfd.perc <- format(round(species.classfd*100/size, 5), nsmall = 2) 
  species.unclassfd.perc <- format(round(species.unclassfd*100/size, 5), nsmall = 2)  
  distinct.species <- length(unique(input$specie))
  
  ## genuses
  gen.classfd <- dim(input[input$genus != "",])[1]
  gen.unclassfd <- dim(input[input$genus == "",])[1]
  gen.classfd.perc <- format(round(gen.classfd*100/size, 5), nsmall = 2)  
  gen.unclassfd.perc <- format(round(gen.unclassfd*100/size, 5), nsmall = 2)  
  distinct.gen <- length(unique(input$genus))
  
  
  ## families
  fam.classfd <- dim(input[input$family != "",])[1]
  fam.unclassfd <- dim(input[input$family == "",])[1]
  fam.classfd.perc <- format(round(fam.classfd*100/size, 5), nsmall = 2)  
  fam.unclassfd.perc <- format(round(fam.unclassfd*100/size, 5), nsmall = 2)  
  distinct.fam <- length(unique(input$family))
  
  ## phyla
  phyla.classfd <- dim(input[input$phylum != "",])[1]
  phyla.unclassfd <- dim(input[input$phylum == "",])[1]
  phyla.classfd.perc <- format(round(phyla.classfd*100/size, 5), nsmall = 2)  
  phyla.unclassfd.perc <- format(round(phyla.unclassfd*100/size, 5), nsmall = 2)  
  distinct.phyla <- length(unique(input$phylum))
  
  
  ## superkingdoms
  spking.classfd <- dim(input[input$superkingdom != "",])[1]
  spking.unclassfd <- dim(input[input$superkingdom == "",])[1]
  spking.classfd.perc <- format(round(spking.classfd*100/size, 5), nsmall = 2)  
  spking.unclassfd.perc <- format(round(spking.unclassfd*100/size, 5), nsmall = 2)  
  distinct.spking <- length(unique(input$superkingdom))
  

  
  descr.df <- rbind(descr.df, c(size,species.classfd, species.unclassfd, species.classfd.perc, species.unclassfd.perc, distinct.species,
                                gen.classfd, gen.unclassfd, gen.classfd.perc, gen.unclassfd.perc, distinct.gen,
                                fam.classfd, fam.unclassfd, fam.classfd.perc, fam.unclassfd.perc, distinct.fam, 
                                phyla.classfd, phyla.unclassfd, phyla.classfd.perc, phyla.unclassfd.perc, distinct.phyla,
                                spking.classfd, spking.unclassfd, spking.classfd.perc, spking.unclassfd.perc, distinct.spking
                                
  ))
  
  
  
}  

descr.df <- descr.df[c(-1),]
write.table(descr.df, file = paste0("original_and_filtered_descriptions.txt"), append = F, row.names = F, col.names = T, sep = "\t", quote = F)






### DISTRIBUTIONs IN ORIGINAL AND FILTERED EPITOPES DATASETS
 
datasets <-c("HG.original","DMG.original", "HG.filtered","DMG.filtered")
 

for(dataset in datasets){
  
  
  if(dataset == "HG.original"){
    input <- read.table(file = "HG_input_enriched.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    
  }else if(dataset == "DMG.original"){
    input <- read.table(file = "DMG_input_enriched.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    
  }else if(dataset == "HG.filtered"){
    input <- read.table(file = "HG_filtered_epitopesNW_list.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    
  }else if(dataset == "DMG.filtered"){
    input <- read.table(file = "DMG_filtered_epitopesNW_list.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    
  }
  
 
   
  ## specie
  tax_group <- "specie"
  distr <- sort(table(input$specie), decreasing = T)
  write.table(distr, file = paste0(dataset,"_distr.",tax_group,".txt"), append = F, row.names = F, col.names = T, sep = "\t", quote = F)
  
  ## genus
  tax_group <- "genus"
  distr <- sort(table(input$genus), decreasing = T)
  write.table(distr, file = paste0(dataset,"_distr.",tax_group,".txt"), append = F, row.names = F, col.names = T, sep = "\t", quote = F)
  
  ## family
  tax_group <- "family"
  distr <- sort(table(input$family), decreasing = T)
  write.table(distr, file = paste0(dataset,"_distr.",tax_group,".txt"), append = F, row.names = F, col.names = T, sep = "\t", quote = F)
  
  ## phylum
  tax_group <- "phylum"
  distr <- sort(table(input$phylum), decreasing = T)
  write.table(distr, file = paste0(dataset,"_distr.",tax_group,".txt"), append = F, row.names = F, col.names = T, sep = "\t", quote = F)
  
  ## superkingdom
  tax_group <- "superkingdom"
  distr <- sort(table(input$superkingdom), decreasing = T)
  write.table(distr, file = paste0(dataset,"_distr.",tax_group,".txt"), append = F, row.names = F, col.names = T, sep = "\t", quote = F)
}




 

 

### COMPRATIVE DISTRIBUTIONs original DMG VS HG

taxgps <- c("specie", "genus", "family", "phylum", "superkingdom")

for(tax_group in taxgps){
  
 
  input2 <- read.table(file = paste0("HG.original_distr.",tax_group,".txt"),sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
  input1 <- read.table(file = paste0("DMG.original_distr.",tax_group,".txt"),sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
  
  head(input1)
  head(input2)
  
  
  freq2 <- c()
  for(i in 1:nrow(input1)){
    x <- input1[i,]$Var1 
    
    
    # if the name from S02 exists in S01, add the freq value from S02
    if (length(input2[input2$Var1 == x,]$Freq)){
      freq2 <- c(freq2, input2[input2$Var1 == x,]$Freq)
      
    }else{ 
      freq2 <- c(freq2, 0)
      
    }
    
  }
 
  input1$freq2 <- freq2
  
  write.table(input1, file = paste0("DMG.vs.HG_original__distr.",tax_group,".txt"), append = F, row.names = F, col.names = T, sep = "\t", quote = F)
}







### COMPRATIVE DISTRIBUTIONs filtered DMG VS HG

taxgps <- c("specie", "genus", "family", "phylum", "superkingdom")

for(tax_group in taxgps){
  
  
  
  input2 <- read.table(file = paste0("HG.filtered_distr.",tax_group,".txt"),sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
  input1 <- read.table(file = paste0("DMG.filtered__distr.",tax_group,".txt"),sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
  
  
  
  freq2 <- c()
  for(i in 1:nrow(input1)){
    x <- input1[i,]$Var1 
     
    # if the name from S02 exists in S01, add the freq value from S02
    if (length(input2[input2$Var1 == x,]$Freq)){
      freq2 <- c(freq2, input2[input2$Var1 == x,]$Freq)
      
    }else{ 
       
      freq2 <- c(freq2, 0)
      
    }
    
    
  }
   
  input1$freq2 <- freq2
   
  
  write.table(input1, file = paste0("DMG.vs.HG_filtered__distr.",tax_group,".txt"), append = F, row.names = F, col.names = T, sep = "\t", quote = F)
}



























 


