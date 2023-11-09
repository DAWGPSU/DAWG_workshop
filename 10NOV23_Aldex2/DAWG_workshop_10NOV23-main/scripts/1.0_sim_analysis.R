###Analysis file for simulated data

#####Starts on slide 19#####

library(ALDEx2)
library(tidyverse)
library(ggplot2)
library(ggpattern)
library(cowplot)

set.seed(12345)
##Reading in data (see "data_simulation.R" for details.)
rdat <- read.csv(file.path("data/simulation/sim_seq_dat.csv"))

##Reading in the simulated flow data for building a scale model
flow_data <- read.csv(file.path("data/simulation/sim_flow_dat.csv"))

##Inspecting elements

## "Y" represents the OTU table
Y <- t(rdat[,-1])
head(Y)

##Vector denoting whether samples was in pre- or post- antibiotic condition.
conds <- as.character(rdat[,1])
head(conds)


## Fitting and analyzing the original ALDEx2 model
mod.base <- aldex(Y, conds) ##equivalent to `gamma = NULL`
mod.base %>% filter(we.eBH < 0.05)


## Recreating ALDEx2
mod.clr <- aldex(Y, conds, gamma = 1e-3)
mod.clr %>% filter(we.eBH < 0.05)

##Checking for concordance in effect sizes
plot(mod.base$effect, mod.clr$effect, xlab = "Original ALDEx2 
  Effect Size", ylab = "CLR Scale Model Effect Size")
abline(a=0,b=1, col = "red", lty = "dashed")


## Adding noise via the default scale model
mod.ss <- aldex(Y, conds, gamma = .25)
mod.ss %>% filter(we.eBH < 0.05)


## Adding more noise via the default scale model
mod.ss.high <- aldex(Y, conds, gamma = 1)
mod.ss.high %>% filter(we.eBH < 0.05)


##Creating an informed model using biological reasoning
scales <- c(rep(1, 50), rep(0.9, 50))
scale_samps <- aldex.makeScaleMatrix(.15, scales, conds, log=FALSE)

mod.know <- aldex(Y, conds, gamma = scale_samps)
mod.know %>% filter(we.eBH < 0.05)

##Now creating an informed model using the flow data
head(flow_data)

flow_data_collapse <- flow_data %>%
  group_by(sample) %>%
  mutate(mean = mean(flow)) %>%
  mutate(stdev = sd(flow)) %>%
  dplyr::select(-flow) %>%
  ungroup() %>%
  unique()

scale_samps <- matrix(NA, nrow = nrow(flow_data_collapse), ncol = 128)
for(i in 1:nrow(scale_samps)){
  scale_samps[i,] <- rnorm(128, flow_data_collapse$mean[i], flow_data_collapse$stdev[i])
}

mod.flow <- aldex(Y, conds, gamma = scale_samps)
mod.flow %>% filter(we.eBH < 0.05)


####Plotting and bringing it all together

##Reading in the true data
dat <- read.csv(file.path("data/simulation/sim_obs_dat.csv"))


## Helper functions
##Function to label True/false positive/negatives
sig_code <- function(sig, Taxa, truth){
  out <- rep("TN", length(Taxa))
  out[sig &(Taxa %in% truth)] <- "TP" # True Positives
  out[sig & (out!="TP")] <- "FP" # False Positives
  out[!sig & (Taxa %in% truth)] <- "FN" # False Negatives
  return(out)
}

##Function to summarize aldex2 output
summary_aldex2 <- function(fit, pval = 0.05){
  fit %>%
    as.data.frame() %>%
    rownames_to_column("category") %>%
    dplyr::select(category, effect, we.ep, we.eBH) %>%
    mutate(padj=we.eBH) %>%
    mutate(mean=effect) %>%
    mutate(low=NA, high=NA) %>%
    mutate(sig = ifelse(padj <= pval, TRUE, FALSE))
}

##Function to create the grid plot
plot_sig2 <- function(rrs, truth, ...){
  names(rrs) <- model.names[names(rrs)]
  bind_rows(rrs, .id="Model") %>%
    dplyr::select(Model, category, sig) %>%
    mutate(Taxa = category) %>%
    mutate(Taxa=as.numeric(sub("Taxa", "", Taxa))) %>%
    mutate(sigcode = sig_code(sig, Taxa, truth)) %>%
    mutate(Taxa=factor(Taxa), sigcode=factor(sigcode,
                                             levels=c("TP", "TN",
                                                      "FP", "FN"))) %>%
    mutate(Model=factor(Model, levels=model.name.levels)) %>%
    ggplot(aes(x=Taxa, y=Model)) +
    geom_tile_pattern(aes(fill=sigcode, pattern = sigcode), color="darkgrey",pattern_fill = 'grey',pattern_colour  = 'grey', pattern_density = 0.015) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          legend.title=element_blank(),
          text = element_text(size=16),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_pattern_manual(values = c(TP = "none", TN = "none", FP = "none", FN = "stripe")) +
    scale_fill_manual(values= c("black", "white", "grey", "white"))
}

##Plotting the results
##Pvalue at default of 0.05

p1 <- gather(dat, Taxa, Count, -Condition) %>%
  mutate(Taxa=as.numeric(sub("Taxa", "", Taxa))) %>%
  mutate(Taxa=factor(Taxa)) %>%
  ggplot(aes(x=Taxa, y=Count)) +
  geom_boxplot(aes(fill = Condition, color = Condition), position=position_dodge(width=1),
               size=1)+
  scale_y_log10() +
  theme_bw() +
  scale_fill_manual(values = c("#fdae61", "#2b83ba")) +
  scale_color_manual(values = c("#fdae61", "#2b83ba")) +
  labs(color='Antibiotic\nTreatment') +
  labs(fill='Antibiotic\nTreatment') +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size=16))

truth <- c(3,4,15,20)##Locations of the differences

model.names <- c("mod.base"="ALDEx2 (Original)",
                 "mod.clr" = "Default Model (Gamma = 1e-3)",
                 "mod.ss"= "Default Model (Gamma = 0.25)","mod.ss.high"= "Default Model (Gamma = 1.00)",
                 "mod.know" = "Knowledged-Based Model",
                 "mod.flow" = "Flow-Based Model")
model.name.levels <- c("Flow-Based Model", "Knowledged-Based Model", "Default Model (Gamma = 1.00)", "Default Model (Gamma = 0.25)",  "Default Model (Gamma = 1e-3)", "ALDEx2 (Original)")

rrs <- list(mod.base=summary_aldex2(mod.base), 
            mod.clr = summary_aldex2(mod.clr),
            mod.ss = summary_aldex2(mod.ss),
            mod.ss.high = summary_aldex2(mod.ss.high),
            mod.know = summary_aldex2(mod.know),
            mod.flow = summary_aldex2(mod.flow))

p2 <- plot_sig2(rrs, truth=truth)
p <- plot_grid(p1, p2, nrow=2, align="v", rel_heights=c(1.7, 1))
p


##Now running a sensitivity analysis over the default scale model

##First, specifying different values for the noise in the scale
gamma_to_test <- c(1e-3, .1, .25, .5, .75, 1, 2, 3, 4, 5)

##Run the CLR function
clr <- aldex.clr(Y, conds)

##Run sensitivity analysis function
sen_res <- aldex.senAnalysis(clr, gamma = gamma_to_test)

##Inspecting results. Note that it is a list with length equal to the number of gammas tested.
##Each element of the list gives the ALDEx2 output.
length(sen_res)
length(gamma_to_test)

head(sen_res[[1]])

##Plotting the sensitivity results.
plotGamma(sen_res, thresh = .1)
