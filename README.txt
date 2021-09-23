The files correspond to the Scripts 1-8 of the dissertation figures 8-10, so that "1 original datasets enrichment.R" corresponds to Script 1 and so forth.

R script "9 plots and distribution files.R" is not mentioned in these figures. It refers to "Dataset comparisons analysis" and includes all the R processing 
needed to create the various plots (venn diagrams, line plots etc).


Short description of each R script:

- the "1 original datasets enrichment.R" script adds taxonomic columns to the original datasets.


- the "2 original datasets filtering.R" script filters the "input_enriched.txt" file in order to find the epitopes that have
  more than 5 identical consecutive amino acids with other epitopes. It creates the file that describes filtered epitopesNW interactions 
  ("filtered_epitopesNW_interactions.txt") and the 3 matrices files.


- the "3 epitopesNW list - attributes.R" script creates epitopesNW attributes file ("filtered_epitopesNW_list.txt").


- the "4 fasta files creation.R" script converts subnetwork csv files to fasta in order for the fasta to be used in CLC (create cladograms).


- the "5 speciesNW interactions.R" script creates file the file that describes speciesNW
  interactions ("speciesNW_interactions.txt").


- the "6 speciesNW list - attributes.R" script creates speciesNW attributes file ("speciesNW_list.txt").


- the "7 proteinsNW interactions.R" script creates file the file that describes proteinsNW
  interactions ("proteinsNW_interactions.txt").


- the "8 proteinsNW list - attributes.R" script creates proteinsNW attributes file ("proteinsNW_list.txt").


- the "9 plots and distribution files.R" script creates: 
  1) Venn diagrams.
  2) the great table (Table 2 in dissertation) that describes HG and DMG original and filtered datasets.
  3) the files for creating the distributions in original and filtered epitopes datasets.
  4) the files for creating comparative distributions original DMG VS HG and filtered DMG VS HG.









