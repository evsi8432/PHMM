library(momentuHMM)
library(circular)
library(CircStats)
library(runner)
library(ggplot2)
library(tidyr)
library(dplyr)
library(mclust)
library(mixreg)
library(readxl)
library(diveMove)
library(boot)
library(data.table)
library(party)
library(numDeriv)
library(data.tree)

setwd("~/Documents/Research/PHMM/src")

#sex <- "Male"
sex <- "Female"

rand_seed <- 1

if(sex == "Male"){
  whales <- c("D21","I107","L87","L88")
} else {
  whales <- c("A100","A113","D26","I129","I145","R48","R58")
}

### run data on test whales ###
conf_matrix <- matrix(rep(0,N_coarse*N_coarse),
                      nrow=N_coarse,
                      ncol=N_coarse)

for(holdout_whale in whales){
  conf_matrix_whale <- read.csv(paste0("../params/",date,"/conf_matrix_",
                                        holdout_whale,"_",
                                        sex,"_",
                                        rand_seed,".csv"))
  
  
  conf_matrix <- conf_matrix + conf_matrix_whale[1:3,2:4]
}

rownames(conf_matrix) <- c("True Resting", "True Travelling", "True Foraging")
colnames(conf_matrix) <- c("Predicted Resting", "Predicted Travelling", "Predicted Foraging")

write.csv(conf_matrix,paste0("../params/",date,"/conf_matrix_",
                             "full","_",
                             sex,"_",
                             rand_seed,".csv"))



### make confusion matrices for each behaviour ###

behaviours <- c("Resting", "Travelling", "Foraging")

for(n in c(1,2,3)){
  
  conf_matrix_behaviour <- matrix(c(0,0,0,0),nrow=2,ncol=2)
  conf_matrix_behaviour[1,1] <- sum(conf_matrix[ n, n])
  conf_matrix_behaviour[1,2] <- sum(conf_matrix[ n,-n])
  conf_matrix_behaviour[2,1] <- sum(conf_matrix[-n, n])
  conf_matrix_behaviour[2,2] <- sum(conf_matrix[-n,-n])
  
  rownames(conf_matrix_behaviour) <- c("True Positive","True Negative")
  colnames(conf_matrix_behaviour) <- c("Predicted Positive", "Predicted Negative")
  
  write.csv(conf_matrix_behaviour,paste0("../params/",date,"/conf_matrix_",
                                         "full","_",
                                          behaviours[n],"_",
                                          sex,"_",
                                          rand_seed,".csv"))
}
