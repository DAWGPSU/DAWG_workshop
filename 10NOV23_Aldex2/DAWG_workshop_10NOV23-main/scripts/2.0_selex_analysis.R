###Code for analyzing the selex data
library(ALDEx2)
library(tidyverse)
library(ggrepel)

##Specifying the conditions. 7 in each the nonselected "NS" and selected "S"
conds <- c(rep("NS", 7), rep("S", 7))
conds

##Examining the selex data
data(selex)
dim(selex)
head(selex)

##Ground truth, from outside experiments
true.pos <- read.delim(file.path("data","selex","selex-truth.txt"))$X
true.pos

##ALDEx2
aldex.mod <- aldex(selex,conds)
aldex.pos <- aldex.mod %>% filter(we.eBH < 0.05)
table(row.names(aldex.pos) %in% true.pos)

##Adding scale
aldex.scale <- aldex(selex, conds, gamma = .5)
aldex.scale.pos <- aldex.scale %>% filter(we.eBH < 0.05)
table(row.names(aldex.scale.pos) %in% true.pos)

##Plotting TP and FP over different levels of gamma
gamma.to.test <- c(1e-3,.25,.5,.75,1,2,3, 4, 5)

TP <- rep(NA, length(gamma.to.test))
FP <- rep(NA, length(gamma.to.test))

#Calculates true positives and false positives for each level of gamma.to.test
for(i in 1:length(gamma.to.test)){
  mod <- aldex(selex,conds,gamma = gamma.to.test[i])
  mod.pos <- row.names(mod %>% filter(we.eBH < 0.05))
  
  TP[i] <- sum(mod.pos %in% true.pos)
  FP[i] <- sum(!(mod.pos %in% true.pos))
  print(i)
}


##Plotting
graph.df <- data.frame("Gamma" = rep(gamma.to.test, 2), "Counts" = c(TP, FP), "Group" = c(rep("TP", length(gamma.to.test)), rep("FP", length(gamma.to.test))))

ggplot(graph.df, aes(x = Gamma, y = Counts, group = Group, color = Group)) + 
  geom_line(aes(linetype = Group), lwd = 1.25) +
  scale_color_manual(values = c("red3", "royalblue3")) +
  scale_linetype_manual(values=c("twodash", "longdash")) +
  theme_bw() +
  ylab("Counts") 

##Built-in sensitivity analysis
clr <- aldex.clr(selex, conds)

##Run sensitivity analysis function
sen_res <- aldex.senAnalysis(clr, gamma = gamma.to.test)

##Plotting the sensitivity results.
plotGamma(sen_res, thresh = .1)

