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
 
K <- as.numeric(args[1]) # number of cross-validations (one means just do all the data)
lambda <- as.numeric(args[2]) # lambda for paper
num_seeds <- as.numeric(args[3]) # number of random seeds

# create directories
dir.create(directory, showWarnings = FALSE)
dir.create(paste0(directory,"/params"), showWarnings = FALSE)
dir.create(paste0(directory,"/plt"), showWarnings = FALSE)

# set seed
set.seed(1)

# load in data
if(load_raw){
  print("loading fine scale data...")
  source("src/load_data_fine.R") # load in Data_fine
} else {
  Data_fine <- data.frame(fread("../../dat/Final_Data_fine1.csv"))
}

# load in coarse data
print("loading coarse scale data...")
source("src/load_data_coarse.R") # load in Data

print("labeling dives and prey captures...")
source("src/label_data.R") # label dives, foraging events

# plot data before fitting model
if(plot & load_raw){
  print("plotting data with labels...")
  source("src/EDA.R")
}

# create cross-validation groups
source("src/make_test_train.R")

# initialize model lists
models_base <- list()
models_PHMM <- list()

probs_PHMM <- list()
probs_base <- list()

AUCs_base <- rep(0,K)
AUCs_PHMM <- rep(0,K)

conf_matrices_base <- list()
conf_matrices_PHMM <- list()

for(k in 1:K){
  
  train_dives <- train_sets[[k]]
  test_dives <- test_sets[[k]]
  
  # fit baseline
  print("fitting baseline...")
  source("src/fit_base.R") 
  models_base[[k]] <- base_model
  
  # evaluate baseline
  print("evaluating baseline...")
  source("src/eval_base.R")  
  
  # fit the PHMM
  print("fitting PHMM...")
  best_hmm <- NULL
  max_ll <- -Inf
  for(rand_seed in 1:num_seeds){
    print(paste("working on fold",k,"of",K,"and seed",rand_seed,"of",num_seeds))
    source("src/fit_PHMM_coarser.R")
    if(-hmm$mod$minimum > max_ll){
      best_hmm <- hmm
      max_ll <- -hmm$mod$minimum
      print("new best hmm")
    }
  }
  hmm <- best_hmm
  models_PHMM[[k]] <- hmm

  # evaluate PHMM
  print("evaluating PHMM...")
  source("src/eval_PHMM_coarser.R")
  
  # plot hmm results
  if(plot){
    print("plotting PHMM...")
    source("src/plot_PHMM_coarser.R")
  }
}

# summarize cross-validation results
print("summarizing results...")
print(lambda)
print(test_sets)
print(probs_PHMM)
source("src/summarize_results.R")

# plot AUCs
source("src/plot_AUCs.R")