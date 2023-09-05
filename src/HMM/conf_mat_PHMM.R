library(momentuHMM)
library(dplyr)
library(mclust)
library(data.table)

#setwd("/Users/evsi8432/Documents/Research/PHMM/src/bash")

# get command-line arguments
args <- commandArgs(trailingOnly=TRUE)

opt_file <- args[1]
args <- as.integer(args[2])

# get options
source(paste0('../opt/',opt_file))

# do the hierarchical thing
if(hier){
  source("../HHMM/conf_mat_PHHMM.R")
  quit()
}

# set seed
set.seed(1)

sind <- 0
if(is.na(args)){
  args_list <- sind:54
} else {
  args_list <- c(args)
}

for(args in args_list){

# Select Model
model_ind <- (args[1] %% 5) + 1
models <- c("no","fixed_1","fixed_2","half_random","random")
model <- models[model_ind]

# select holdout whale
whale_ind <- floor(args[1] / 5) + 1
whales <- c("A100","A113","D21","D26","I107","I129","I145","L87","L88","R48","R58")
holdout_whale <- whales[whale_ind]

print(model)
print(holdout_whale)

### BEGIN COMPUTATION ###

source("../HMM/load_data.R")

# get (un)labelled Data for held out whale
Data_labeled <- Data[Data$ID %in% holdout_whale,]
Data_unlabeled <- Data[Data$ID %in% holdout_whale,]

Data_unlabeled$label <- 7
Data_unlabeled <- prepData(Data_unlabeled,coordNames=NULL)

# load in hmm
make_title <- function(start,end){
  title <- paste0(start,statesPerBehaviour[1])
  for(nstates in statesPerBehaviour[2:3]){
    title <- paste0(title,"-",nstates)
  }
  for(feature in features1){
    title <- paste0(title,"-",feature)
  }
  if(length(sex) > 1){
    title <- paste0(title,"_all")
  } else {
    title <- paste0(title,"_",sex)
  }
  title <- paste0(title,"_",end)
  return(title)
}

# load in best hmm
files <- Sys.glob(make_title(paste0(directory,"/params/"),
                             paste0(model,"-",holdout_whale,"-*-hmm.rds")))

best_hmm <- NULL
best_nll <- Inf
for(file in files){
  hmm <- readRDS(file)
  if(hmm$mod$minimum < best_nll){
    best_hmm <- hmm
    best_nll <- hmm$mod$minimum
  }
}
hmm <- best_hmm

# get new hmm with held-out whale
Par0 <- getPar0(hmm)
eps <- 1e-8

if(model == "no"){
  hmm0 <- fitHMM(data=Data_unlabeled,
                 nbStates=N,
                 dist=hmm$conditions$dist,
                 DM=hmm$conditions$DM,
                 beta0=Par0$beta,
                 delta0=(1-eps)*(Par0$delta)+eps*rep(1/N,N),
                 Par0=Par0$Par,
                 userBounds=userBounds,
                 workBounds=workBounds,
                 nlmPar = list('stepmax'=1e-100,
                               'iterlim'=1))
} else {
  hmm0 <- fitHMM(data=Data_unlabeled,
               nbStates=N,
               dist=hmm$conditions$dist,
               fixPar=list(label = hmm$conditions$fixPar$label),
               DM=hmm$conditions$DM,
               beta0=Par0$beta,
               delta0=(1-eps)*(Par0$delta)+eps*rep(1/N,N),
               Par0=Par0$Par,
               userBounds=userBounds,
               workBounds=workBounds,
               nlmPar = list('stepmax'=1e-100,
                             'iterlim'=1))
}
# add pairs to hmm0
hmm0$pairs <- hmm$pairs

# decode states
probs <- stateProbs(hmm0)
inds <- cumsum(statesPerBehaviour)

rest_inds <- hmm0$pairs[1:inds[1],2]
Data_labeled$prob_resting <- rowSums(probs[,rest_inds,drop=FALSE])

trav_inds <- hmm0$pairs[(inds[1]+1):inds[2],2]
Data_labeled$prob_travelling <- rowSums(probs[,trav_inds,drop=FALSE])

forg_inds <- hmm0$pairs[(inds[2]+1):inds[3],2]
Data_labeled$prob_foraging <- rowSums(probs[,forg_inds,drop=FALSE])

# just get held-out whale
Data_labeled_small <- Data_labeled[Data_labeled$knownState != 4,c("prob_resting",
                                                                  "prob_travelling",
                                                                  "prob_foraging",
                                                                  "knownState")]

### run data on test whales ###
conf_matrix <- matrix(c(0,0,0,
                        0,0,0,
                        0,0,0),
                      nrow=3,ncol=3)
for(i in 1:3){
  for(j in 1:3){
    conf_matrix[i,j] = sum(Data_labeled_small[Data_labeled_small$knownState == i,j])
  }
}

rownames(conf_matrix) <- c("True Resting", "True Travelling", "True Foraging")
colnames(conf_matrix) <- c("Predicted Resting", "Predicted Travelling", "Predicted Foraging")

# Save hmm
make_title <- function(start,end){
  title <- paste0(start,statesPerBehaviour[1])
  for(nstates in statesPerBehaviour[2:3]){
    title <- paste0(title,"-",nstates)
  }
  for(feature in features1){
    title <- paste0(title,"-",feature)
  }
  if(length(sex) > 1){
    title <- paste0(title,"_all")
  } else {
    title <- paste0(title,"_",sex)
  }
  title <- paste0(title,"_",end)
  return(title)
}

write.csv(conf_matrix,
          make_title(paste0(directory,"/params/"),
                     paste0(model,"-",
                            holdout_whale,"-",
                            "confusion_matrix.csv")))

}
