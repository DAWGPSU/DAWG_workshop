########################
###Data_simulation.R####
###Michelle Nixon#######
###October 30, 2023#####
########################

###This file can be used to create the simulated data set used in the workshop.
set.seed(12345)
##Loading libraries
library(tidyverse)

##Function to create the true abundances via Poisson resampling
##d = mean vector of abundances in the form c(mean of all taxa for condition 1, mean of all taxa for condition 2)
##n = number of samples per condition
##Note that the number of taxa is d/2
create_true_abundances <- function(d, n){
  dd <- length(d)/2
  dat <- d %>%
    sapply(function(x) rpois(n, lambda=x)) %>%
    t() %>%
    as.data.frame() %>%
    split(rep(1:2, each=dd)) %>%
    purrr::map(~`rownames<-`(.x, paste0("Taxa", 1:dd))) %>%
    purrr::map(t) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    cbind(Condition=factor(rep(c("Pre", "Post"), each=n), 
                           levels = c("Pre", "Post")), .) %>% `rownames<-`(., NULL)
  return(dat)
}

##Function to resample data to an arbitrary size
##dat = true abundances from `create_true_abundances()`
##seq.depth = sequencing depth
resample_data <- function(dat, seq.depth){
  ddat <- as.matrix(dat[,-1])/rowSums(as.matrix(dat[,-1]))
  for (i in 1:nrow(dat)){
    dat[i,-1] <- rmultinom(1, size=seq.depth, prob=ddat[i,])
  }
  return(dat)
}

###Setting the data parameters for the simulation
##Denotes the mean for the 20 taxa
##Note only taxa 3, 4, 15, and 20 change
d.pre <- c(4000, 4000, 4000, 4000, 4000,
           400,400,400,400,4000,400,500,500,500,
           400,400,400,400,400,400)
d.post <- d.pre
d.post[c(3,4,15,20)] <- c(3000, 2000, 200, 100)

#Combining together
d <- c(d.pre, d.post)

##Create true abundances. 50 samples per condition
dat <- create_true_abundances(d, n=50)
##Resampling data to an arbitrary sequencing depth
rdat <- resample_data(dat, seq.depth=5000)
##Saving
write.csv(rdat, file.path("data/simulation/sim_seq_dat.csv"), row.names = FALSE)
write.csv(dat, file.path("data/simulation/sim_obs_dat.csv"), row.names = FALSE)

##Now we are going to create mock flow cytometry data
##Creating three flow measurements for each sample
##Based on the true totals (which we know)
## Finding sample totals
totals <- rowSums(dat[,-1])

##Function to create our mock flow data. Note we use a normal distribution with standard deviation 5e2.
flow_cytometry <- function(totals, samp.names, replicates = 3){
  samp.names <- rep(samp.names, each = replicates)
  flow_vals <- sapply(totals, FUN = function(totals,   
                                             replicates){rnorm(replicates,totals,5e2)}, 
                      replicates = replicates, simplify = TRUE)
  flow_data <- data.frame("sample" = samp.names, "flow" = c(flow_vals))
  return(flow_data)
}

##Sampling for flow data
flow_data <- flow_cytometry(totals, rownames(dat))

##Saving
write.csv(flow_data, file.path("data/simulation/sim_flow_dat.csv"), row.names = FALSE)
