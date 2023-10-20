###Introduction to decontam ### 

##Using data stored in decontam package ###

##https://benjjneb.github.io/decontam/vignettes/decontam_intro.html#identifying-contaminants-in-marker-gene-and-metagenomics-data


###RAN on R 4.2.3 ##

###packages versions decontam 1.18.0 ; ggplot2 3.4.4, phyloseq 1.42.0

#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("decontam")

#library(decontam)

#library(phyloseq); packageVersion("phyloseq")
## [1] '1.42.0'
#library(ggplot2); packageVersion("ggplot2")

## [1] '3.4.4'
#library(decontam); packageVersion("decontam")
## [1] '1.18.0'


###Install and load packages ####

list.of.packages <-
  c("decontam",
    "phyloseq",
    "ggplot2"
  )

#Install new packages and load all packages 
new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]

if (length(new.packages))
  install.packages(new.packages)

packages_load <-
  lapply(list.of.packages, require, character.only = TRUE)


####Check for problems loading the packages 


if (any(as.numeric(packages_load) == 0)) {
  warning(paste("Package/s: ", paste(list.of.packages[packages_load != TRUE], sep = ", "), "not loaded!"))
} else {
  print("All packages were successfully loaded.")
}


#Remove packages object from global environment
rm(list.of.packages, new.packages, packages_load)



###Read in 16S phyloseq object from decontam package
ps <- readRDS(system.file("extdata", "MUClite.rds", package="decontam"))
ps


###Phyloseq object has 1951 ASVs that were inferred using DADA2 from 16S amplicon based sequencing 

head(sample_data(ps))

###This lets you look at the metadata component of your phyloseq object. It is important here to have a column for quantitative reading which is your DNA concentration determined by fluorescent intensity measurement and a Sample_or_Control column which notes which of your samples are negative controls and which are samples and which are negative controls


###Look at the number of reads you have in each sample 

df <- as.data.frame(sample_data(ps)) # Puts information from phyloseq object into a dataframe
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()


##Number of reads from your samples is between 10,000 and 40,000 reads 
###You do not want to toss or filter out your low read samples 
###Want to use the negative controls to identify contaminants especially in the raw reads

##Frequency Method to Identify Contamination ####
###The distribution of the frequency of each sequence is used as function of the input DNA concentration to identify contamimants 

contamdf.freq <- isContaminant(ps, method="frequency", conc="quant_reading")
head(contamdf.freq)


##Your output is a data frame the column called p shows the probability that was used to classify if this particular sequence is a contamimant 

###The contaminant column has TRUE or FALSE in it 
##  A designation of TRUE means that statistically the sequence is a contaminant that exceeds the threshold
##The default threshold is threshold = 0.1 which means that in order for the contaminant column to be labelled TRUE 
#p < 0.1

table(contamdf.freq$contaminant)

###This means that 58 of the ASVS are classified as contaminants

head(which(contamdf.freq$contaminant))

##Tells us which of our sequences are contaminants for example sequence 3, 30, 53 etc 


###Let's plot a contaminant (Seq3) next to a non contaminant (Seq1)

plot_frequency(ps, taxa_names(ps)[c(1,3)], conc="quant_reading") + 
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")


###In this plot the black line is the model of noncontaminant sequences which is the frequency of which we would expect the sequence to be independent of the DNA concentration 
##Red line is the model of contaminantion the frequency is exected to be inversely proportion to the input concentration of DNA (contaminating DNA makes up a larger amount of DNA in samples that have have smaller amounts of total DNA)

##Seq3 fits the red line 


##Inspecting some more reads to ensure they look like what we would expect to see 

set.seed(100)
plot_frequency(ps, taxa_names(ps)[sample(which(contamdf.freq$contaminant),3)], conc="quant_reading") +
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

###These look like what we would expect to see for contaminants

###Now that we have likely identified contaminants let's remove them from our phyloseq object. 

ps ###our original phyloseq object 

ps.noncontam <- prune_taxa(!contamdf.freq$contaminant, ps) ###prune the phyloseq object get rid of the contaminants 
ps.noncontam

###This step allows you to remove the contaminants from your phyloseq object. You can then export this phyloseq object or use the ps.noncontam object for downstream analysis 



#### Prevalance Method to ID contaminants####

##Another way to identify contaminants is to look at prevalence
##This looks at the presence or absence across samples 
##Each sequence in a sample is compared to the prevalence of that sequence in the negative controls to identify contaminants

sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

head(which(contamdf.prev$contaminant))

##This method identified more potential contaminants than the frequency method but it also missed some abundant contaminant sequences like Seq3
##the default threshold for this method is 0.1 
##Let's try again with threshold = 0.5
##This will ID all sequences that are more prevalent in negative controls than samples 

contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

##Visualizing the times these taxa were observed in controls and samples

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

##Branching of taxa that are on contaminants and noncontaminants

