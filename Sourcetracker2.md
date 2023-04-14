# Sourcetracker 2
Sourceracker2 is a tool to study microbial source using the Gibb's sampler. As a simple explanation to the Gibb's sampling, it randomly draws an instance from the distribution of each variable, and this is conditional to the other values in order to estimate the complex joint distributions of the entire sample. Sourcetracker2 uses  the Gibb's sample as a underlying mechanism to determine the original source to a sink. It works in 4 steps: 
- 1. randomly assign sequences in a sink to source environment so that the random assignment represents the source that a given sequence sink came from. 
- 2. Select one of the sequences from step 1, calculate the actual probabilities of that sequence having come from any of the source environments, and update the assigned source environment of the sequence based on a random draw with the calculated probabilities. Repeat many times.
- At intervals in the repeats of step 2, take the source environment assingments of all the sequences in a sink and record them.
- After doing step 3 a certain number of times (i.e. recording the full assignments of source environments for each sink sequence), terminate the iteration and move to the next sink.
Here's a fuller [explanation of the machinary and the raw code](https://github.com/biota/sourcetracker2/blob/master/ipynb/Sourcetracking%20using%20a%20Gibbs%20Sampler.ipynb) for sourcetracker2. This page will go through the installation and usage of sourcetracker2 to track the source of each sequences in order to assess the contamination proportion. 

# Installation
```
conda create -n st2 -c biocore python=3.8 numpy scipy scikit-bio biom-format h5py hdf5 seaborn
```
Restart terminal
```
conda activate st2
```
```
pip install sourcetracker
```
To test if the installation is successful: 
```
sourcetracker2 gibbs --help
```

# Download data files from previous sessions (MEGAN result)
```
cd scratch
mkdir sourcetracker2
cd sourcetracker2
wget --no-check-certificate --content-disposition "https://raw.githubusercontent.com/SusanTian/DAWG_workshop/main/Sourcetracker-Georgia-Filtered-MappingFile28MAR23.txt" -O Sourcetracker-Georgia-Filtered-MappingFile28MAR23.txt
wget --no-check-certificate --content-disposition "https://github.com/SusanTian/DAWG_workshop/blob/main/feature-table.biom?raw=true" -O feature-table.biom
```

# Sourcetracker2 Script
```
nano sourcetracker2.sh
```
```
#run with plaque and modern calculus as Oral (combine the two)
sourcetracker2 gibbs -i feature-table.biom \
-m Sourcetracker-Georgia-Filtered-MappingFile28MAR23.txt --alpha2 1.000 \
-o sourcetracker-results-10000-alpha2-ORAL --sink_rarefaction_depth 1000 --source_rarefaction_depth 1000
```
```
bash sourcetracker2.sh
```

# Plotting using R
```
wget --no-check-certificate --content-disposition "https://raw.githubusercontent.com/SusanTian/DAWG_workshop/main/sourcetracker-results.csv" -O sourcetracker-results.csv
```
```
# summarize sourcetracker results 

#### set working directory ####
setwd("/storage/home/s/sbt5355/scratch/sourcetracker2")

#### load packages ####
library(ggplot2)
library(dplyr)
library(tidyr)

#### Read in Tables ####
result <- read.csv("sourcetracker-results.csv")


#### Plot the data ####
result %>% 
  pivot_longer(!SampleID, names_to = "Environment", values_to = "Proportion") %>% 
  ggplot()+
  geom_bar(aes(y=Proportion, x=SampleID, fill = Environment), stat = "identity") +
  labs(x = "Samples", y = "Percentage") +
  theme_bw() +
  scale_fill_manual(name = "Environment", values = c("lightskyblue", "thistle", "firebrick", "darkseagreen","peachpuff", "gray80"))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14)
  )
