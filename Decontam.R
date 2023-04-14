# Introduction to decontam
Decontam package facilitates the identification and visualization of contaminants, enabling their elimination and creation of more precide representation
of sampled communites and metagenomic data

# Neccessary Input files
The primary element of deconatm is a feature table obtained from the original data, containing the relative abundace of sequence features in each
sample. These features could encompass a broad range of types, such as amplicon sequence variants (ASVs), operational taxonomic units (OTUs), 
taxonomic groups or phylotypes (such as genera) or metagenome-assembled genomes (MAGs). The second file required is group of "negative control" 
sampple, where the sequencing was conducted on blanks devoid of any biogical specimen.

# Setting up
```
Create phyloseq object and install packages in R
library(phyloseq); packageVersion("phyloseq") ## '1.34.0'
library(ggplot2); packageVersion("ggplot2") ## '3.3.3' 
library(decontam); packageVersion("decontam") ## '1.8.0'    #BiocManager::install("decontam")
library(plyr); packageVersion("plyr")
library(metagMisc); packageVersion("metagMisc")
```
#Import Biom table
```
data<-import_biom("table.biom")
```
#Import metadata
```
metadata<-import_qiime_sample_data("metadata.txt")
```
#Create phyloseq object
Phyloseq objects can contain biom tables and metadata
```
phyloseq.object<-merge_phyloseq(data,metadata)
```
Check phyloseq object
```
phyloseq.object
```
#Identify Contaminants- prevelance
We will use the "prevelance'method  that look for absence/presence
```
sample_data(phyloseq.object)$is.neg<-sample_data(phyloseq.object)$SampleType== "Control"
```
Run the isContaminant () function to identify which species/taxa are classified as negative control
```
contamdf.prev<-isContaminant(phyloseq.object, method = "prevalence", neg = "is.neg")
```
#Check the table to prevelance of species/taxa identified as contaminants
```
table(contamdf.prev$contaminant)
```
#The threshhold to identify contaminant is raised from default ( 0.1) to 0.5
```
contamdf.prev05<-isContaminant(phyloseq.object, method = "prevalence", neg = "is.neg", threshold = 0.5)
```
#Result Visualization
Make phyloseq object of presence-absence in negative controls and true samples 
```
ps.pa<-transform_sample_counts(phyloseq.object, function(abund) 1*(abund>0))
ps.pa.neg<-prune_samples(sample_data(ps.pa)$SampleType == "Control", ps.pa)
ps.pa.pos<-prune_samples(sample_data(ps.pa)$SampleType == "Calculus", ps.pa)
```
# Make data.frame of prevalence in positive and negative samples
```
df.pa.05<-data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg), contaminant=contamdf.prev05$contaminant) 
ggplot(data=df.pa.05, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (Calculus Samples)")
  ```
# Based on the plot, we see a split in the blue and orange dots.

#Prune taxa
# Now that we have identified which taxa are likely contaminants, we can remove those taxa from our biom table. 
# This will provide a 'cleaner' table and should improve the robustness of our downstream analyses. 
  ```
phyloseq.decontam.05<-prune_taxa(!df.pa.05$contaminant, phyloseq.object)
phyloseq.decontam.05
  ``` 

#Convert phyloseq objects into txt file 

# Now that you have a clean dataset. You may want to analyze this table in another program (e.g., QIIME2)
# To do that, you may want to export the table as a txt file 
# First, convert the biom table in the phyloseq object to a dataframe  
  ```
species_postdecontam.05.df<-phyloseq_to_df(phyloseq.decontam.05, addtax = FALSE, addtot = FALSE, addmaxrank = FALSE, sorting = "NULL")
View(species_postdecontam.05.df)
  ```
# Second, write the table to a txt file 
  ```
write.table(species_postdecontam.05.df, "data-postdecontam-prev05.txt", sep = "\t", row.names = FALSE)
  ```
# check to see if the text file was created 

# You can now proceed to perform downstream analyses 

