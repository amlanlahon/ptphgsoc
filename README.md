# ptphgsoc

### RUN THE CODE ACCORDING TO THE LIST ###

#1 GEO2R script.R

#2 Microarray Analysis.R
  #Use covdesc.txt in your working directory for the treatment value. It is case sensitive.
  #After getting the .annot file from GPL10348, create a annotation.csv file and replicate the dataframe in the working directory before the Annotation section.

#3 Removal of Non-annotated DEGs.R
  #filter the non-annotated genes from the annotation.csv file and create a new .csv file called annotated_genes.csv

#4 Heatmap.R

#5 volcano_plot.R
  #Use the two_fold_change.csv file as your base file for this analysis by making a copy and renaming it along the lines of the code.
  
#6 WGCNA.R
  #You will need the list of annotated DEGs and then average every duplicate. Conditional Filter the duplicates out and average their expression values.
  
#7 PTP Normalization and Heatmap.R
  #Create a list of PTP genes using https://febs.onlinelibrary.wiley.com/doi/10.1111/febs.13748, table 1.
