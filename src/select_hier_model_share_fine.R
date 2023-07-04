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

sex <- "Female"
date <- "2023-04-03"
holdout_whale <- "None"

best_hhmm <- NULL
max_like <- -Inf
best_seed <- NULL

for(rand_seed in 1:5){
  hhmm <- readRDS(paste0("../params/",date,"/hierHmm_", 
                         holdout_whale,"_",
                         sex,"_",
                         rand_seed,".rds"))
  if (-hhmm$mod$minimum > max_like){
    best_hhmm <- hhmm
    max_like <- -hhmm$mod$minimum
    best_seed <- rand_seed
  }
}

