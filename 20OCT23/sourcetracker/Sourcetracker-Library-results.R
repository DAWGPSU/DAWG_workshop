
# Sourcetracker analysis

library(ggplot2)
library(li)
setwd("/Users/svw5689/Desktop/microArch/RESEARCH_ASSISTANTSHIP/Library_Project/Sourcetracker/")

data<- read.csv("sourcetracker-results-Hungary-Gobero.csv", header = TRUE)
ST_results<-melt(data, id=c("SampleID"))
#ST_results$variable<-factor(ST_results$variable, levels = ST_results$variable[order(ST_results$value)])

ggplot() + geom_bar(aes(y = ST_results$value, x = ST_results$SampleID, fill = ST_results$variable),
                    stat = "identity") +
  labs(x = "Samples", y = "Percentage") +
  theme_bw() +
  scale_fill_manual(name = "Environment", values = c("tan4", "lightskyblue", "mediumvioletred", "darkseagreen","lightsalmon3", "gray44")) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 70, hjust = 1, size =12),
        text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.background = element_blank())


#### Gobero ####
data<- read.csv("Gobero-ST-RESULTS.csv", header = TRUE)
ST_results<-melt(data, id=c("SampleID"))
#ST_results$variable<-factor(ST_results$variable, levels = ST_results$variable[order(ST_results$value)])

plot1<-ggplot() + geom_bar(aes(y = ST_results$value, x = ST_results$SampleID, fill = ST_results$variable),
                    stat = "identity") +
  labs(x = "Samples", y = "Proportion") +
  theme_bw() +
  scale_fill_manual(name = "Environment", values = c("tan4", "lightskyblue", "mediumvioletred", "darkseagreen","lightsalmon3", "gray44")) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 70, hjust = 1, size =12),
        text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.background = element_blank())

#### Hungary ####
data2<- read.csv("sourcetracker-results-Hungary.csv", header = TRUE)
ST_results2<-melt(data2, id=c("SampleID"))
#ST_results$variable<-factor(ST_results$variable, levels = ST_results$variable[order(ST_results$value)])

plot2<-ggplot() + geom_bar(aes(y = ST_results2$value, x = ST_results2$SampleID, fill = ST_results2$variable),
                           stat = "identity") +
  labs(x = "Samples", y = "Proportion") +
  theme_bw() +
  scale_fill_manual(name = "Environment", values = c("tan4", "lightskyblue", "mediumvioletred", "darkseagreen","lightsalmon3", "gray44")) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 70, hjust = 1, size =12),
        text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.background = element_blank())

ggarrange(plot1, plot2,
          labels = c("A", "B"),
          common.legend = TRUE,
          legend = "bottom",
          ncol=2)
