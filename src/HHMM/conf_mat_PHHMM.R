library(momentuHMM)
library(dplyr)
library(mclust)
library(data.table)
library(RcppHungarian)

#setwd("/Users/evsi8432/Documents/Research/PHMM/src/bash")

args <- commandArgs(trailingOnly=TRUE)

opt_file <- args[1]
args <- as.integer(args[2])

# get options
source(paste0('../opt/',opt_file))

# define way to make titles
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

sind <- 0

if(is.na(args)){
  args_list <- sind:(19*6-1)
} else {
  args_list <- c(args)
}

for(args in args_list){

# set seed
set.seed(1)

# Select Model
model_ind <- (args[1] %% 6) + 1
models <- c("no","fixed_1","fixed_2","half_random","random","random_2")
model <- models[model_ind]

# select holdout whale
whale_ind <- floor(args[1] / 6) + 1
whales <- c("A100a","A100b","A113a","A113b",
            "D21a","D21b","D26a","D26b",
            "I107a","I107b","I129","I145a","I145b",
            "L87","L88","R48a","R48b","R58a","R58b")
holdout_whale <- whales[whale_ind]

print(model)
print(holdout_whale)

### BEGIN COMPUTATION ###

# get data
source("../HMM/load_data.R")

# get (un)labelled Data for held out whale
Data_labeled <- Data[Data$ID %in% holdout_whale,]
Data_unlabeled <- Data[Data$ID %in% holdout_whale,]
Data_unlabeled$label[Data_unlabeled$level %in% 2] <- 7
Data_unlabeled <- prepData(Data_unlabeled,
                           coordNames=NULL,
                           hierLevels=c("1","2i","2"))

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

# add distribution
if(model != "no"){
  dist$label = paste0("cat",nlevels(Data$label))
}

# update user and work bounds and fixPar
userBounds <- list()
workBounds <- list()
fixPar <- list()
for(feature in names(dist)){
  userBounds[[feature]] <- hmm$conditions$userBounds[[feature]]
  workBounds[[feature]] <- hmm$conditions$workBounds[[feature]]
  fixPar[[feature]] <- hmm$conditions$fixPar[[feature]]
}

# fit new HMM
hmm0 <- fitHMM(data=Data_unlabeled,
               hierStates=hmm$conditions$hierStates,
               hierDist=hmm$conditions$hierDist,
               hierBeta=getPar0(hmm)$hierBeta,
               hierDelta=getPar0(hmm)$hierDelta,
               Par0=getPar0(hmm)$Par,
               fixPar=fixPar,
               userBounds=userBounds,
               workBounds=workBounds,
               DM=hmm$conditions$DM,
               nlmPar = list('stepmax'=1e-100,
                             'iterlim'=1))

# get new indices (match best indices with each index set)
inds <- c(0,cumsum(statesPerBehaviour))
rest_inds <- (inds[1]+1):inds[2]
trav_inds <- (inds[2]+1):inds[3]
forg_inds <- (inds[3]+1):inds[4]

# decode states
probs0 <- stateProbs(hmm0)
Data_labeled$prob_resting <- rowSums(probs0[,rest_inds,drop=FALSE])
Data_labeled$prob_travelling <- rowSums(probs0[,trav_inds,drop=FALSE])
Data_labeled$prob_foraging <- rowSums(probs0[,forg_inds,drop=FALSE])

# just get held-out whale
Data_labeled_small <- Data_labeled[Data_labeled$knownState != 4,c("prob_resting",
                                                                  "prob_travelling",
                                                                  "prob_foraging",
                                                                  "boutnum",
                                                                  "knownState")]
Data_labeled_small <- Data_labeled_small[!is.na(Data_labeled_small$boutnum),]
Data_labeled_small <- Data_labeled_small %>% distinct(boutnum,.keep_all=TRUE)

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
write.csv(conf_matrix,
          make_title(paste0(directory,"/params/"),
                     paste0(model,"-",
                            holdout_whale,"-",
                            "confusion_matrix.csv")))
}
