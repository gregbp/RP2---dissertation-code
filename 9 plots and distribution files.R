library(VennDiagram)     
library(seqinr)
library(plotly)
library(sets)
library(plyr) 
library(RColorBrewer)
 



### VennDiagrams
library(VennDiagram)   
library(RColorBrewer)

# datasets <-c("original","epitopesNW")
# datasets <-c("original")
# datasets <-c("epitopesNW")
# dataset<-c("original")
datasets <-c("original","filtered")

for(dataset in datasets){
  
  setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged\\output\\checked res")
  if(dataset == "original"){
    input1 <- read.table(file = "S01.TOTAL_with_taxon_levels.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    input2 <- read.table(file = "S02.TOTAL_with_taxon_levels.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    
    
    
  }else if(dataset == "filtered"){
    input1 <- read.table(file = "nwS01.epitopes_nodes_attributes.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    input2 <- read.table(file = "nwS02.epitopes_nodes_attributes.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    
    
    
  }
  
  
  setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged\\output\\checked res\\venn diagrams\\last")
  cat("\n\ndataset =",dataset,"\n\n")
  
  
  
  
  ## VennDiagram 1 - epitopes
  S01.epitopes <- input1$qseq
  S02.epitopes <- input2$qseq
  
  length(S01.epitopes)
  length(unique(S01.epitopes))
  length(S02.epitopes)
  length(unique(S02.epitopes))
  
  myCol <- brewer.pal(3, "Pastel2")
  
  venn.diagram(
    x = list(S01.epitopes, S02.epitopes),
    category.names = c(paste0(dataset," HG - epitopes") , paste0(dataset, " DMG - epitopes")),
    filename = paste0(dataset,'_1_epitopes_venn_diagramm.png'),
    output=TRUE,
    # fill = c("green","red"),
    # Output features
    imagetype="png" ,
    # height = 1980 ,
    # width = 1980 ,
    resolution = 500,
    compression = "lzw",
    height = 2500 ,
    width = 2500,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = c("green","red"),
    
    cex = 1,
    fontfamily = "sans",
    cat.cex = 1,
    cat.default.pos = "outer",
    cat.pos = c(-00.1, 00.1)
    # cat.dist = c(0.005, 0.005)
    
    
  )
  # break
  
  ## checks
  # length(unique(S01.epitopes))
  # length(unique(S02.epitopes))
  
  
  
  ## VennDiagram 2 - proteins
  S01.proteins <- input1$subjectAccVer
  S02.proteins <- input2$subjectAccVer
  
  venn.diagram(
    x = list(S01.proteins, S02.proteins),
    category.names = c(paste0(dataset," HG - proteins") , paste0(dataset, " DMG - proteins")),
    filename = paste0(dataset,'_2_proteins_venn_diagramm.png'),
    output=TRUE,
    # fill = c("green","red"),
    # Output features
    imagetype="png" ,
    # height = 1980 ,
    # width = 1980 ,
    resolution = 500,
    compression = "lzw",
    height = 2500 ,
    width = 2500,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = c("green","red"),
    
    cex = 1,
    fontfamily = "sans",
    cat.cex = 1,
    cat.default.pos = "outer",
    cat.pos = c(-00.1, 00.1)
    # cat.dist = c(0.005, 0.005)
    
  ) 
  
  ## checks
  # length(unique(S01.proteins))
  # length(unique(S02.proteins))
  
  
  ## VennDiagram 3 - species
  S01.species <- input1$specie
  S02.species <- input2$specie
  
  
  venn.diagram(
    x = list(S01.species, S02.species),
    category.names = c(paste0(dataset," HG - species") , paste0(dataset, " DMG - species")),
    filename = paste0(dataset,'_3_species_venn_diagramm.png'),
    output=TRUE,
    # fill = c("green","red"),
    # Output features
    imagetype="png" ,
    # height = 1980 ,
    # width = 1980 ,
    resolution = 500,
    compression = "lzw",
    height = 2500 ,
    width = 2500,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = c("green","red"),
    
    cex = 1,
    fontfamily = "sans",
    cat.cex = 1,
    cat.default.pos = "outer",
    cat.pos = c(-20, 20),
    cat.dist = c(0.05, 0.05)
    
    
  ) 
  
  ## checks
  # length(unique(S01.species))
  # length(unique(S02.species))
  
  ## VennDiagram 4 - genera
  S01.genera <- input1$genus
  S02.genera <- input2$genus
  
  
  venn.diagram(
    x = list(S01.genera, S02.genera),
    category.names = c(paste0(dataset," HG - genera") , paste0(dataset, " DMG - genera")),
    filename = paste0(dataset,'_4_genera_venn_diagramm.png'),
    output=TRUE,
    # fill = c("green","red"),
    # Output features
    imagetype="png" ,
    # height = 1980 ,
    # width = 1980 ,
    resolution = 500,
    compression = "lzw",
    height = 2500 ,
    width = 2500,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = c("green","red"),
    
    cex = 1,
    fontfamily = "sans",
    cat.cex = 1,
    cat.default.pos = "outer",
    cat.pos = c(-20, 20),
    cat.dist = c(0.05, 0.05)
    
    
  ) 
  
  
  ## checks
  # length(unique(S01.genera))
  # length(unique(S02.genera))
  
  
  ## VennDiagram 5 - families
  S01.families <- input1$family
  S02.families <- input2$family
  
  
  venn.diagram(
    x = list(S01.families, S02.families),
    category.names = c(paste0(dataset," HG - families") , paste0(dataset, " DMG - families")),
    filename = paste0(dataset,'_5_families_venn_diagramm.png'),
    output=TRUE,
    # fill = c("green","red"),
    # Output features
    imagetype="png" ,
    # height = 1980 ,
    # width = 1980 ,
    resolution = 500,
    compression = "lzw",
    height = 2500 ,
    width = 2500,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = c("green","red"),
    
    cex = 1,
    fontfamily = "sans",
    cat.cex = 1,
    cat.default.pos = "outer",
    cat.pos = c(-20, 20),
    cat.dist = c(0.05, 0.05)
    
    
  ) 
  
  
  
  ## checks
  # length(unique(S01.families))
  # length(unique(S02.families))
  
  ## VennDiagram 6 - phyla
  S01.phyla <- input1$phylum
  S02.phyla <- input2$phylum
  
  venn.diagram(
    x = list(S01.phyla, S02.phyla),
    category.names = c(paste0(dataset," HG - phyla") , paste0(dataset, " DMG - phyla")),
    filename = paste0(dataset,'_6_phyla_venn_diagramm tmp.png'),
    output=TRUE,
    # fill = c("green","red"),
    # Output features
    imagetype="png" ,
    # height = 1980 ,
    # width = 1980 ,
    resolution = 500,
    compression = "lzw",
    height = 2500 ,
    width = 2500,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = c("green","red"),
    
    cex = 1,
    fontfamily = "sans",
    cat.cex = 1,
    cat.default.pos = "outer",
    cat.pos = c(-20, 20),
    cat.dist = c(0.05, 0.05)
    
    
  ) 
  
  
  ## checks
  # length(unique(S01.phyla))
  # length(unique(S02.phyla))
  
  ## VennDiagram 7 - superkingdoms
  S01.superkingdoms <- input1$superkingdom
  S02.superkingdoms <- input2$superkingdom
  
  
  venn.diagram(
    x = list(S01.superkingdoms, S02.superkingdoms),
    category.names = c(paste0(dataset," HG - superkingdoms") , paste0(dataset, " DMG - superkingdoms")),
    filename = paste0('z ',dataset,'_7_superkingdoms_venn_diagramm.png'),
    output=TRUE,
    # fill = c("green","red"),
    # Output features
    imagetype="png" ,
    # height = 1980 ,
    # width = 1980 ,
    resolution = 500,
    compression = "lzw",
    height = 2500 ,
    width = 2500,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = c("green","red"),
    
    cex = 1,
    fontfamily = "sans",
    cat.cex = 1,
    cat.default.pos = "outer",
    # cat.pos = c(-0.01, 0.01)
    cat.pos = c(-27, 27),
    cat.dist = c(0.05, 0.05)
    
    
  ) 
  
  
  
  ## checks
  # length(unique(S01.superkingdoms))
  # length(unique(S02.superkingdoms))
  
  
  
}  


  


### GREAT TABLE
 



##### test
setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged\\output\\checked res")

dataset <- "S01"
nw.epitopes <- read.table(file = "nwS01.epitopes_nodes_attributes.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")

length(unique(nw.epitopes$specie))
length(unique(nw.epitopes$queryAccVer))
length(unique(nw.epitopes$subjectAccVer))
length(unique(nw.epitopes$node_ID))



# ## DESCRIPTION INITIAL DATASET

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



datasets <-c("S01.initial","S02.initial", "S01.epitopes","S02.epitopes")
# datasets <-c("S01.initial","S02.initial")
# datasets <-c("S01.epitopes","S02.epitopes")
# datasets <-c("S01.initial")


for(dataset in datasets){
  
  setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged\\output\\checked res")
  if(dataset == "S01.initial"){
    input <- read.table(file = "S01.TOTAL_with_taxon_levels.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    
  }else if(dataset == "S02.initial"){
    input <- read.table(file = "S02.TOTAL_with_taxon_levels.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    
  }else if(dataset == "S01.epitopes"){
    input <- read.table(file = "nwS01.epitopes_nodes_attributes.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    
  }else if(dataset == "S02.epitopes"){
    input <- read.table(file = "nwS02.epitopes_nodes_attributes.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    
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
  
  # size
  # species.classfd
  # species.unclassfd
  # species.classfd.perc
  # species.unclassfd.perc
  # distinct.species
  
  
  
  descr.df <- rbind(descr.df, c(size,species.classfd, species.unclassfd, species.classfd.perc, species.unclassfd.perc, distinct.species,
                                gen.classfd, gen.unclassfd, gen.classfd.perc, gen.unclassfd.perc, distinct.gen,
                                fam.classfd, fam.unclassfd, fam.classfd.perc, fam.unclassfd.perc, distinct.fam, 
                                phyla.classfd, phyla.unclassfd, phyla.classfd.perc, phyla.unclassfd.perc, distinct.phyla,
                                spking.classfd, spking.unclassfd, spking.classfd.perc, spking.unclassfd.perc, distinct.spking
                                
  ))
  
  
  
}  

descr.df <- descr.df[c(-1),]
descr.df
setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged\\distributions")
write.table(descr.df, file = paste0("initial_and_epitopesNW_descriptions.txt"), append = F, row.names = F, col.names = T, sep = "\t", quote = F)






# ## DISTRIBUTIONs IN INITIAL DATASETS & EPITOPES NWS
setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged\\output\\checked res")

datasets <-c("S01.initial","S02.initial", "S01.epitopes","S02.epitopes")
# datasets <-c("S01.initial","S02.initial")
# datasets <-c("S01.epitopes","S02.epitopes")
# dataset <-c("S01.initial")

for(dataset in datasets){
  
  setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged\\output\\checked res")
  
  
  if(dataset == "S01.initial"){
    input <- read.table(file = "S01.TOTAL_with_taxon_levels.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    
  }else if(dataset == "S02.initial"){
    input <- read.table(file = "S02.TOTAL_with_taxon_levels.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    
  }else if(dataset == "S01.epitopes"){
    input <- read.table(file = "nwS01.epitopes_nodes_attributes.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    
  }else if(dataset == "S02.epitopes"){
    input <- read.table(file = "nwS02.epitopes_nodes_attributes.txt",sep="\t",header=TRUE, quote = "",fill=T, comment.char = "")
    
  }
  
  
  head(input)
  size <- dim(input)[1]
  print(size)
  
  
  setwd("D:\\Research Project 2 RP2\\RP2 code\\final res\\merged")
  
  ## specie
  tax_group <- "specie"
  distr <- sort(table(input$specie), decreasing = T)
  
  # names(distr[2]) <- "NA"
  head(distr)
  print("unclassified")
  print(tax_group)
  print(distr[[]])
  
  
  
  
  
  
  write.table(distr, file = paste0("distributions\\", dataset,"_distr.",tax_group,".txt"), append = F, row.names = F, col.names = T, sep = "\t", quote = F)
  
  for(i in 1:nrow(distr)){
    distr[[i]] <- distr[[i]]/size*100
    
  }
  
  # head(distr)
  write.table(distr, file = paste0("distributions\\", dataset,"_distr.",tax_group,"_perc.txt"), append = F, row.names = F, col.names = T, sep = "\t", quote = F)
  
  
  
  ## genus
  tax_group <- "genus"
  distr <- sort(table(input$genus), decreasing = T)
  
  # names(distr[2]) <- "NA"
  head(distr)
  print("unclassified")
  print(tax_group)
  print(distr[[]])
  
  write.table(distr, file = paste0("distributions\\", dataset,"_distr.",tax_group,".txt"), append = F, row.names = F, col.names = T, sep = "\t", quote = F)
  
  for(i in 1:nrow(distr)){
    distr[[i]] <- distr[[i]]/size*100
    
  }
  
  # head(distr)
  write.table(distr, file = paste0("distributions\\", dataset,"_distr.",tax_group,"_perc.txt"), append = F, row.names = F, col.names = T, sep = "\t", quote = F)
  
  
  
  
  ## family
  tax_group <- "family"
  distr <- sort(table(input$family), decreasing = T)
  
  # names(distr[2]) <- "NA"
  # head(distr)
  print("unclassified")
  print(tax_group)
  print(distr[[]])
  
  write.table(distr, file = paste0("distributions\\", dataset,"_distr.",tax_group,".txt"), append = F, row.names = F, col.names = T, sep = "\t", quote = F)
  
  for(i in 1:nrow(distr)){
    distr[[i]] <- distr[[i]]/size*100
    
  }
  
  # head(distr)
  write.table(distr, file = paste0("distributions\\", dataset,"_distr.",tax_group,"_perc.txt"), append = F, row.names = F, col.names = T, sep = "\t", quote = F)
  
  
  
  
  ## phylum
  tax_group <- "phylum"
  distr <- sort(table(input$phylum), decreasing = T)
  
  # names(distr[2]) <- "NA"
  # head(distr)
  print("unclassified")
  print(tax_group)
  print(distr[[]])
  
  write.table(distr, file = paste0("distributions\\", dataset,"_distr.",tax_group,".txt"), append = F, row.names = F, col.names = T, sep = "\t", quote = F)
  
  for(i in 1:nrow(distr)){
    distr[[i]] <- distr[[i]]/size*100
    
  }
  
  # head(distr)
  write.table(distr, file = paste0("distributions\\", dataset,"_distr.",tax_group,"_perc.txt"), append = F, row.names = F, col.names = T, sep = "\t", quote = F)
  
  
  
  ## superkingdom
  tax_group <- "superkingdom"
  distr <- sort(table(input$superkingdom), decreasing = T)
  
  # names(distr[2]) <- "NA"
  # head(distr)
  print("unclassified")
  print(tax_group)
  print(distr[[]])
  print("")
  
  write.table(distr, file = paste0("distributions\\", dataset,"_distr.",tax_group,".txt"), append = F, row.names = F, col.names = T, sep = "\t", quote = F)
  
  for(i in 1:nrow(distr)){
    distr[[i]] <- distr[[i]]/size*100
    
  }
  
  # head(distr)
  write.table(distr, file = paste0("distributions\\", dataset,"_distr.",tax_group,"_perc.txt"), append = F, row.names = F, col.names = T, sep = "\t", quote = F)
  
  
  
  
}



































 


