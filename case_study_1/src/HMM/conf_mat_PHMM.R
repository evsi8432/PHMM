library(momentuHMM)
library(dplyr)
library(mclust)
library(data.table)
library(RcppHungarian)

opt_file <- "logMDDD_1-1-1_dd-30_2023-10-23.R"

# get options
source(paste0('../opt/',opt_file))

# define whales
whales <- c("A100a","A100b","A113a","A113b",
            "D21a","D21b","D26a","D26b",
            "I107a","I107b","I129","I145a","I145b",
            "L87","L88","R48a","R48b","R58a","R58b")
n_whales <- length(whales)

# define models
Data <- data.frame(fread("../../dat/case_study_1_data.csv"))

ratio <- sum(!(Data$knownState %in% 4)) / nrow(Data)
ratio <- round(ratio,3)
models <- list()
models[[1]] <- c("fixed",  0.0)       # no weight
models[[2]] <- c("fixed",  0.5*ratio) 
models[[3]] <- c("fixed",  ratio)     # equal weight
models[[4]] <- c("fixed",  0.5 + 0.5*ratio)
models[[5]] <- c("fixed",  1.0)       # natural weight
n_models <- length(models) 

# set seed
set.seed(1)

sind <- 0
args_list <- sind:(n_whales*n_models-1)

for(args in args_list){

# Set Model
model_ind <- (floor(args[1]) %% n_models) + 1
model <- models[[model_ind]][1]
lamb  <- as.numeric(models[[model_ind]][2])

# Select Holdout Whale
whale_ind <- (floor(args[1] / n_models) %% n_whales) + 1
holdout_whale <- whales[whale_ind]

print(model)
print(lamb)
print(holdout_whale)

### BEGIN COMPUTATION ###

Data <- data.frame(fread("../../dat/case_study_1_data.csv"))

# get (un)labelled Data for held out whale
Data_labeled <- Data[Data$ID %in% holdout_whale,]
Data_unlabeled <- Data[Data$ID %in% holdout_whale,]
Data_unlabeled$label <- 4
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
                             paste0(model,"-",lamb,"-",holdout_whale,"-*-hmm.rds")))

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

# get new indices (match best indices with each index set)
inds <- c(0,cumsum(statesPerBehaviour))
rest_inds <- (inds[1]+1):inds[2]
trav_inds <- (inds[2]+1):inds[3]
forg_inds <- (inds[3]+1):inds[4]

# decode states with new HMM
probs0 <- stateProbs(hmm0)
Data_labeled$prob_resting <- rowSums(probs0[,rest_inds,drop=FALSE])
Data_labeled$prob_travelling <- rowSums(probs0[,trav_inds,drop=FALSE])
Data_labeled$prob_foraging <- rowSums(probs0[,forg_inds,drop=FALSE])

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

# Save confusion matrix
write.csv(conf_matrix,
          make_title(paste0(directory,"/params/"),
                     paste0(model,"-",
                            lamb,"-",
                            holdout_whale,"-",
                            "confusion_matrix.csv")))

write.csv(Data_labeled_small,
          make_title(paste0(directory,"/params/"),
                     paste0(model,"-",
                            lamb,"-",
                            holdout_whale,"-",
                            "probs_labs.csv")))
}
