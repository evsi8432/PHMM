#.rs.restartR()
library(Rcpp)
library(tools)
#compileAttributes("/Users/evsi8432/Documents/Research/momentuHMM")
#package_native_routine_registration_skeleton("/Users/evsi8432/Documents/Research/momentuHMM")
#install.packages("/Users/evsi8432/Documents/Research/momentuHMM",
#                 repos=NULL,
#                 type="source")
#library(momentuHMM)
#detach("package:momentuHMM", unload=TRUE)
library(momentuHMM)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(GGally)
library(signal)
library(oce)
library(pROC)
library(latex2exp)

directory <- "/Users/evsi8432/Documents/Research/PHMM/case_study_2"
setwd(directory)

# set options

plot <- T # whether to plot results
load_raw <- F # whether to load raw Data from scratch

args = commandArgs(trailingOnly=TRUE)

K <- 1 # number of cross-validations (one means just do all the data)
lambda <- 0.01 # lambda for paper
num_seeds <- 10 # number of random seeds

# set distributions
dist <- list()
dist[["delt_d"]] <- "norm"
#dist[["rajp"]] <- "vm"
dist[["htv"]] <- "gamma"
dist[["jp_normed"]] <- "gamma"

print("fitting PHMM...")
best_hmm <- NULL
max_ll <- -Inf

for(rand_seed in 1:num_seeds){
  
  directory <- "/Users/evsi8432/Documents/Research/PHMM/case_study_2"
  
  model_name <- paste0("hmm_",
                       paste0(setdiff(names(dist),"knownState"),collapse = "_"),"_",
                       log10(lambda),"_",
                       1,"_",1,"_",
                       rand_seed,".rds")
  
  hmm <- readRDS(paste0(directory,"/params/",model_name))
  
  print(-hmm$mod$minimum)
  
  if(-hmm$mod$minimum > max_ll){
    best_hmm <- hmm
    max_ll <- -hmm$mod$minimum
    print("new best hmm")
  }
}
hmm <- best_hmm

print(hmm)
