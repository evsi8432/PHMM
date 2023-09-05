library(momentuHMM)
library(dplyr)
library(mclust)
library(data.table)
library(RcppHungarian)

#setwd("/Users/evsi8432/Documents/Research/PHMM/src/bash")

# get command-line arguments
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

# do the hierarchical thing
if(hier){
  source("../HHMM/fit_model_PHHMM.R")
  quit()
}

sind <- 0

if(is.na(args)){
  args_list <- sind:(60*n_retries-1)
} else {
  args_list <- c(args)
}

for(args in args_list){

# Set Seed
rand_seed <- (args[1] %% n_retries) + 1
set.seed(rand_seed)

# select holdout whale
whale_ind <- (floor(args[1] / n_retries) %% 12) + 1
whales <- c("none","A100","A113","D21","D26",
            "I107","I129","I145","L87","L88","R48","R58")
holdout_whale <- whales[whale_ind]

# Select Model
model_ind <- floor(args[1] / (12*n_retries)) + 1
models <- c("no","fixed_1","fixed_2","half_random","random")
model <- models[model_ind]

print(rand_seed)
print(holdout_whale)
print(model)

# check if we already have a model
file <- make_title(paste0(directory,"/params/"),
                   paste0(model,"-",
                          holdout_whale,"-",
                          rand_seed,"-",
                          "hmm.rds"))
if(file.exists(file)){
  next
}

# get data
source("../HMM/load_data.R")

# get rid of held out whale
Data <- Data[!(Data$ID %in% holdout_whale),]

# set up label distribution
if(model != "no"){
  dist$label = paste0("cat",nlevels(Data$label))
}

# set up DM matrix
DM <- list()

if(!(model %in% c("no"))){
  old_hmm <- readRDS(make_title(paste0(directory,"/params/"),
                                paste0("no","-",
                                       holdout_whale,"-",
                                       rand_seed,"-",
                                       "hmm.rds")))
  inds <- cumsum(statesPerBehaviour)
  b_inds <- list()
  b_inds$resting <- old_hmm$pairs[1:inds[1],2]
  b_inds$travelling <- old_hmm$pairs[(inds[1]+1):inds[2],2]
  b_inds$foraging <- old_hmm$pairs[(inds[2]+1):inds[3],2]
  behaviours <- c("resting","travelling","foraging")
}

if (model %in% c("fixed_1")){

  labelDM <- c()
  for(b in 1:3){
    for(i in 1:N){
      if(i %in% b_inds[[behaviours[b]]]){
        labelDM <- c(labelDM, c(1,0))
      } else {
        labelDM <- c(labelDM, c(0,1))
      }
    }
    labelDM <- c(labelDM,rep(c(0,0),N)) # pseudolabels
  }
  labelDM <- matrix(labelDM,6*N,2,byrow = TRUE)
  DM$label <- labelDM

} else if (model %in% c("fixed_2")){

  labelDM <- c()
  for(b in 1:3){
    for(temp in 1:2){ # do for pseudolabels and real labels
      for(i in 1:N){
        if(i %in% b_inds[[behaviours[b]]]){
          labelDM <- c(labelDM, c(1,0))
        } else {
          labelDM <- c(labelDM, c(0,1))
        }
      }
    }
  }
  labelDM <- matrix(labelDM,6*N,2,byrow = TRUE)
  DM$label <- labelDM

} else if (model %in% c("half_random","random")){

  labelDM <- c()
  for(b in 1:3){
    for(i in 1:N){ # real labels
      if(i %in% b_inds[[behaviours[b]]]){
        labelDM <- c(labelDM, c(1,0,0,0))
      } else {
        labelDM <- c(labelDM, c(0,0,0,1))
      }
    }
    for(i in 1:N){ # pseudolabels
      if(i %in% b_inds[[behaviours[b]]]){
        labelDM <- c(labelDM, c(0,1,0,0))
      } else {
        labelDM <- c(labelDM, c(0,0,1,0))
      }
    }
  }
  labelDM <- matrix(labelDM,6*N,4,byrow = TRUE)
  DM$label <- labelDM
}

for(feature in features1){
  if(dist[[feature]] == "mvnorm2"){
    DM0 <- diag(5*N)
    for(i in 1:N){
      DM0[3*N+i,2*N+i] <- 0.5
      DM0[3*N+i,4*N+i] <- 0.5
    }
    DM[[feature]] <- DM0
  } else if (dist[[feature]] == "norm"){
    DM[[feature]] <- list(mean = ~1, sd = ~1)
  } else if (substring(dist[[feature]], 1,3) == "cat"){
    ncats <- as.integer(substring(dist[[feature]], 4))
    DM[[feature]] <- list()
    for(i in 1:(ncats-1)){
      DM[[feature]][[paste0("prob",i)]] = ~1
    }
  }
}

# Set fixed Parameters
fixPar <- list()
if(model %in% c("fixed_1","fixed_2")){
  fixPar$label <- c(0,-100)
} else if (model == "half_random"){
  fixPar$label <- c(NA,NA,NA,-100)
} else if (model == "random"){
  fixPar$label <- c(NA,NA,NA,NA)
}

# Set initial Parameters
Par0 <- list()

for(feature in features1){
  if(dist[[feature]] == "norm"){
    Par0[[feature]] <- matrix(rep(NA,2*N), nrow = N)
  } else if (dist[[feature]] == "mvnorm2"){
    Par0[[feature]] <- matrix(rep(NA,5*N), nrow = N)
  } else if (substring(dist[[feature]], 1,3) == "cat") {
    ncats <- as.integer(substring(dist[[feature]], 4))
    Par0[[feature]] <- matrix(rep(NA,(ncats-1)*N), nrow = N)
  } else {
    print("feature distribtuion not recognized")
  }
}

if(model == "no"){
  for(i in 1:N){
    for(feature in features1){
      if(dist[[feature]] == "norm"){
        Par0[[feature]][i,1] <- mean(Data[,feature]) + rnorm(1)*sd(Data[,feature]) # mean
        Par0[[feature]][i,2] <- exp(0.5*rnorm(1))/N * sd(Data[,feature])           # sd
      } else if (dist[[feature]] == "mvnorm2") {
        featurex <- paste0(feature,".x")
        featurey <- paste0(feature,".y")
        Par0[[feature]][i,1] <- mean(Data[,featurex]) + rnorm(1)*sd(Data[,featurex]) # mean
        Par0[[feature]][i,2] <- mean(Data[,featurey]) + rnorm(1)*sd(Data[,featurey]) # mean
        Par0[[feature]][i,3] <- log((exp(0.5*rnorm(1))/N) * sd(Data[,featurex]))     # sd x
        Par0[[feature]][i,4] <- 0.0                                                  # cov xy
        Par0[[feature]][i,5] <- log((exp(0.5*rnorm(1))/N) * sd(Data[,featurey]))     # sd y
      } else if (substring(dist[[feature]], 1,3) == "cat") {
        ncats <- as.integer(substring(dist[[feature]], 4))
        for(j in 1:(ncats-1)){
          Par0[[feature]][i,j] <- mean(Data[,feature] == j , na.rm=T)
        }
        noise <- rexp(ncats)
        noise <- noise / sum(noise)
        Par0[[feature]][i,] = 0.9*Par0[[feature]][i,] + 0.1*head(noise,-1)
      }
    }
  }
} else {
  Par0 <- getPar0(old_hmm)$Par
}

if(model %in% c("fixed_1","fixed_2")){
  Par0$label = fixPar$label
} else if(model == "half_random"){
  Par0$label <- c(qlogis(0.1),qlogis(0.1),qlogis(0.001),-100)
} else if(model == "random"){
  Par0$label <- c(qlogis(0.1),qlogis(0.1),qlogis(0.001),qlogis(0.001))
}

# Set initial TPM
if(model %in% "no"){
  beta0 <- matrix(rnorm(N*(N-1),mean=-2,sd=1),nrow=1)
} else {
  beta0 <- getPar0(old_hmm)$beta
}

# Set initial delta
if(model %in% "no"){
  delta0 <- matrix(rep(1/N,N),nrow=1)
} else {
  delta0 <- getPar0(old_hmm)$delta
}

# prep data
Data0 <- prepData(Data,coordNames=NULL)

eps <- 1e-8

# fit HMM
hmm <- fitHMM(data=Data0,
              nbStates=N,
              dist=dist,
              DM=DM,
              beta0=beta0,
              delta0=(1-eps)*delta0 + eps*rep(1/N,N),
              Par0=Par0,
              fixPar=fixPar,
              userBounds=userBounds,
              workBounds=workBounds,
              nlmPar = list('steptol'=1e-4,
                            'iterlim'=1000))

print(hmm)

# get labels switching values
Data <- cbind(Data,stateProbs(hmm))

accs <- matrix(nrow = N, ncol = N)
behaviours <- c("resting","travelling","foraging")

# make accuracy matrix
for(i in 1:N){
  for(j in 1:N){
    behaviour_ind <- which(i <= cumsum(statesPerBehaviour))[1]
    behaviour_data <- Data[Data$knownState %in% behaviour_ind,]
    accs[i,j] <- sum(behaviour_data[,paste("state",j)])
  }
}

# get pairs
if (model %in% "no"){
  hs <- HungarianSolver(-accs)
  pairs <- hs$pairs
} else {
  pairs <- old_hmm$pairs
}

# order pairs by increasing mean of first feature
feature <- features1[1]
inds <- c(0,cumsum(statesPerBehaviour))
new_pairs <- matrix(nrow = N, ncol = 2)
new_pairs[,1] <- seq(1,N)

for (i in 1:3){
  behaviour_inds <- seq(inds[i]+1,inds[i+1])
  old_behaviour_pairs <- pairs[behaviour_inds,2]
  mus <- c()
  for (j in old_behaviour_pairs) {
    if (dist[[feature]] == "mvnorm2"){
      mus <- c(mus,hmm$mle[[feature]]["mean.x",j])
    } else if (dist[[feature]] == "norm"){
      mus <- c(mus,hmm$mle[[feature]]["mean",j])
    } else if (substring(dist[[feature]], 1,3) == "cat"){
      mus <- c(mus,hmm$mle[[feature]]["prob1",j])
    }
  }
  new_pairs[behaviour_inds,2] <- old_behaviour_pairs[order(mus)]
}
hmm$pairs <- new_pairs

print(hmm$pairs)

# Save hmm
saveRDS(hmm,
        make_title(paste0(directory,"/params/"),
                   paste0(model,"-",
                          holdout_whale,"-",
                          rand_seed,"-",
                          "hmm.rds")))

}
