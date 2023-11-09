#.rs.restartR()
library(Rcpp)
library(tools)
compileAttributes("/Users/evsi8432/Documents/Research/momentuHMM")
package_native_routine_registration_skeleton("/Users/evsi8432/Documents/Research/momentuHMM")
install.packages("/Users/evsi8432/Documents/Research/momentuHMM",
                 repos=NULL,
                 type="source")
library(momentuHMM)
detach("package:momentuHMM", unload=TRUE)
library(momentuHMM)

library(dplyr)
library(mclust)
library(data.table)
library(RcppHungarian)
library(mvtnorm)

setwd("/Users/evsi8432/Documents/Research/PHMM/src/bash")

# get command-line arguments
args <- commandArgs(trailingOnly=TRUE)
args <- c("logMDDD_1-1-1_dd-30_2023-10-23.R",NA)

opt_file <- args[1]
args <- as.integer(args[2])

# get options
source(paste0('../opt/',opt_file))

# do the hierarchical thing
if(hier){
  source("../HHMM/fit_model_PHHMM.R")
  quit()
}

# define whales
whales <- c("none","A100a","A100b","A113a","A113b",
            "D21a","D21b","D26a","D26b",
            "I107a","I107b","I129","I145a","I145b",
            "L87","L88","R48a","R48b","R58a","R58b")
n_whales <- length(whales)

# define models
source("../preprocessing/load_data.R")
ratio <- sum(!(Data$knownState %in% 4)) / nrow(Data)
ratio <- round(ratio,3)
models <- list()
models[[1]] <- c("no",     1.0)
models[[2]] <- c("random", 1.0)
models[[3]] <- c("fixed",  0.0)       # no weight
models[[4]] <- c("fixed",  0.5*ratio) 
models[[5]] <- c("fixed",  ratio)     # equal weight
models[[6]] <- c("fixed",  0.5 + 0.5*ratio)
models[[7]] <- c("fixed",  1.0)       # natural weight
n_models <- length(models) 

sind <- 0

if(is.na(args)){
  args_list <- sind:(n_retries*n_whales*n_models-1)
} else {
  args_list <- c(args)
}

for(args in args_list){

# Set Model
model_ind <- (floor(args[1]) %% n_models) + 1
model <- models[[model_ind]][1]
lamb  <- as.numeric(models[[model_ind]][2])

# Select Holdout Whale
whale_ind <- (floor(args[1] / n_models) %% n_whales) + 1
holdout_whale <- whales[whale_ind]
  
# Select Seed
rand_seed <- (floor(args[1] / (n_models*n_whales)) %% n_retries) + 1
set.seed(rand_seed)

print(rand_seed)
print(holdout_whale)
print(model)
print(lamb)

# check if we already have a model
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

file <- make_title(paste0(directory,"/params/"),
                   paste0(model,"-",
                          lamb,"-",
                          holdout_whale,"-",
                          rand_seed,"-",
                          "hmm.rds"))
if(file.exists(file)){
  print(paste(file,"already exists. continuing..."))
  next
}

# get data and remove heldout whale
source("../preprocessing/load_data.R")
Data <- Data[!(Data$ID %in% holdout_whale),]

# define states
bhavs <- c("rest","trav","forg")
behaviours <- c("resting","travelling","foraging")

inds <- c(0,cumsum(statesPerBehaviour))
b_inds <- list()
b_inds$resting <- (inds[1]+1):inds[2]
b_inds$travelling <- (inds[2]+1):inds[3]
b_inds$foraging <- (inds[3]+1):inds[4]

# set up label distribution
dist$label <- NULL
if(model != "no"){
  dist$label = paste0("cat",nlevels(Data$label))
}

# set up DM matrix
DM <- list()

if (model %in% c("fixed")){

  labelDM <- c()
  for(b in 1:3){
    for(i in 1:N){
      if(i %in% b_inds[[behaviours[b]]]){
        labelDM <- c(labelDM, c(1,0))
      } else {
        labelDM <- c(labelDM, c(0,1))
      }
    }
  }
  labelDM <- matrix(labelDM,3*N,2,byrow = TRUE)
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
if(model %in% c("fixed")){
  fixPar$label <- c(0,-100)
}

# Set initial Parameters
Par0 <- list()

for(feature in features1){
  if(dist[[feature]] == "norm"){
    Par0[[feature]] <- matrix(rep(NA,2*N), nrow = N)
    Par0[[feature]][,1] <- mean(Data[,feature]) + 0.1*rnorm(N)*sd(Data[,feature]) # mean
    Par0[[feature]][,2] <- exp(0.1*rnorm(N))/sqrt(N) * sd(Data[,feature])           # sd
  } else if (dist[[feature]] == "mvnorm2") {
    Par0[[feature]] <- matrix(rep(NA,5*N), nrow = N)
    featurex <- paste0(feature,".x")
    featurey <- paste0(feature,".y")
    Par0[[feature]][,1] <- mean(Data[,featurex]) + 0.1*rnorm(N)*sd(Data[,featurex]) # mean
    Par0[[feature]][,2] <- mean(Data[,featurey]) + 0.1*rnorm(N)*sd(Data[,featurey]) # mean
    Par0[[feature]][,3] <- 0.1*rnorm(N) + log(sd(Data[,featurex])) - 0.1*log(N) # sd x
    Par0[[feature]][,4] <- 0.0                                                  # cov xy
    Par0[[feature]][,5] <- 0.1*rnorm(N) + log(sd(Data[,featurey])) - 0.1*log(N) # sd y
  } else if (substring(dist[[feature]], 1,3) == "cat") {
    ncats <- as.integer(substring(dist[[feature]], 4))
    Par0[[feature]] <- matrix(rep(NA,(ncats-1)*N), nrow = N)
    for(j in 1:(ncats-1)){
      Par0[[feature]][,j] <- mean(Data[,feature] == j , na.rm=T)
    }
    noise <- matrix(rexp(ncats*N), nrow = N)
    noise <- noise / rowSums(noise)
    Par0[[feature]] = 0.9*Par0[[feature]] + 0.1*noise[,-1]
  }
}

# get accuracy of mixture model for each behaviour and state
accs <- matrix(rep(0,N^2),nrow=N,ncol=N)
for(i in 1:N){

  # get data associated with the behaviour in question
  behaviour_ind <- which(i <= cumsum(statesPerBehaviour))[1]
  behaviour_data <- Data[Data$knownState %in% behaviour_ind,]

  # get p(X,Y)
  pXY <- matrix(rep(0,N*nrow(behaviour_data)),
                nrow=nrow(behaviour_data),ncol=N)
  for(j in 1:N){
    for(feature in features1){
      if(dist[[feature]] == "norm"){
        pXY[,j] <- pXY[,j] + dnorm(behaviour_data[,feature],
                                   mean=Par0[[feature]][j,1],
                                   sd=exp(Par0[[feature]][j,2]) + 0.01)
      } else if (dist[[feature]] == "mvnorm2") {
        featurex <- paste0(feature,".x")
        featurey <- paste0(feature,".y")
        meanx <- Par0[[feature]][j,1]
        meany <- Par0[[feature]][j,2]
        sigxx <- exp(Par0[[feature]][j,3]) + 0.01
        sigyy <- exp(Par0[[feature]][j,5]) + 0.01
        sigxy <- exp(0.5*Par0[[feature]][j,3] +
                     0.5*Par0[[feature]][j,5] +
                    -exp(-Par0[[feature]][j,4]))
        pXY[,j] <- pXY[,j] + dmvnorm(behaviour_data[,c(featurex,featurey)],
                                     mean=c(meanx,meany),
                                     sigma=matrix(c(sigxx,sigxy,sigxy,sigyy),nrow=2))
      } else if (substring(dist[[feature]], 1,3) == "cat") {
        probs <- Par0[[feature]][j,]
        probs <- c(probs,1.0-sum(probs))
        pXY[,j] <- pXY[,j] + probs[behaviour_data[,feature]]
      }
    }
  }
  # get p(X|Y) and add to accuracies
  p_X_given_Y <- pXY / rowSums(pXY)
  accs[i,] <- colSums(p_X_given_Y) / nrow(behaviour_data)
}

hs <- HungarianSolver(-accs)
pairs <- hs$pairs

for(feature in features1){
  Par0[[feature]] <- Par0[[feature]][pairs[,2],]
}

m_prob <- 0.001

if (model %in% c("fixed")){
  Par0$label = fixPar$label
} else if (model == "random") {
  Par0$label <- c(3*mean(Data[,"label"] == 1, na.rm=T), m_prob, m_prob, 
                  m_prob, 3*mean(Data[,"label"] == 2, na.rm=T), m_prob,  
                  m_prob, m_prob, 3*mean(Data[,"label"] == 3, na.rm=T))
}

# Set initial gamma and delta
beta0 <- matrix(rnorm(N*(N-1),mean=-2,sd=1),nrow=1)
delta0 <- matrix(rep(1/N,N),nrow=1)

# prep data
Data0 <- prepData(Data,coordNames=NULL)

eps <- 1e-8

# define the known states
knownStates <- Data$knownState
knownStates[knownStates == 4] <- NA
if (model %in% c("no","random")){
  knownStates <- rep(NA,length(knownStates))
}

# fit HMM
hmm <- fitHMM(data=Data0,
              nbStates=N,
              dist=dist,
              DM=DM,
              beta0=beta0,
              delta0=(1-eps)*delta0 + eps*rep(1/N,N),
              Par0=Par0,
              fixPar=fixPar,
              knownStates=knownStates,
              lambda=lamb,
              userBounds=userBounds,
              workBounds=workBounds,
              nlmPar = list('iterlim'=iterlim,
                            'print.level'=0))

# find the best pairings and refit the model
probs <- stateProbs(hmm)
hmm$data <- cbind(hmm$data,probs)
accs <- matrix(nrow = N, ncol = N)
behaviours <- c("resting","travelling","foraging")
for(i in 1:N){
  for(j in 1:N){
    behaviour_ind <- which(i <= cumsum(statesPerBehaviour))[1]
    behaviour_data <- hmm$data[hmm$data$knownState %in% behaviour_ind,]
    accs[i,j] <- sum(behaviour_data[,paste("state",j)])
  }
}
hs <- HungarianSolver(-accs)
pairs <- hs$pairs

# fix Par0
oldPar0 <- getPar0(hmm)
newPar0 <- list()
for(feature in features1){
  if(dist[feature] == "norm"){
    oldPar0feature <- matrix(oldPar0$Par[[feature]],nrow=N)
    newPar0[[feature]] <- matrix(rep(NA,2*N), nrow = N)
    for(i in 1:N){
      newPar0[[feature]][pairs[i,1],] <- oldPar0feature[pairs[i,2],]
    }
  } else if (dist[feature] == "mvnorm2") {
    oldPar0feature <- matrix(oldPar0$Par[[feature]],nrow=N)
    newPar0[[feature]] <- matrix(rep(NA,5*N), nrow=N)
    for(i in 1:N){
      newPar0[[feature]][pairs[i,1],] <- oldPar0feature[pairs[i,2],]
    }
  } else if (substring(dist[[feature]], 1,3) == "cat") {
    ncats <- as.integer(substring(dist[[feature]], 4))
    oldPar0feature <- matrix(oldPar0$Par[[feature]],nrow=N)
    newPar0[[feature]] <- matrix(rep(NA,(ncats-1)*N), nrow=N)
    for(i in 1:N){
      newPar0[[feature]][pairs[i,1],] <- oldPar0feature[pairs[i,2],]
    }
  }
}

if(model %in% c("fixed")){
  newPar0$label <- oldPar0$Par$label
} else if (model %in% c("random")) {
  oldPar0label <- matrix(oldPar0$Par$label,nrow=N)
  newPar0$label <- matrix(rep(NA,3*N),nrow=N)
  for(i in 1:N){
    newPar0$label[pairs[i,1],] <- oldPar0label[pairs[i,2],] + 1e-8
  }
}

# fix delta
newDelta <- oldPar0$delta[pairs[,2]]

# fix beta
oldBeta <- oldPar0$beta
newBeta <- c()
for(i in 1:(N-1)){
  newBeta <- c(newBeta,0,oldBeta[(N*(i-1)+1):(N*i)])
}
newBeta <- c(newBeta,0)
newBeta <- matrix(newBeta,nrow=N,byrow=T)
newBeta <- newBeta[pairs[,2],pairs[,2]]
newBeta0 <- c()
for(i in 1:N){
  newBeta0 <- c(newBeta0,newBeta[i,-i])
}
newBeta <- matrix(newBeta0,nrow=1)

# refit the hmm
hmm0 <- fitHMM(data=Data0,
               nbStates=N,
               dist=dist,
               DM=DM,
               beta0=newBeta,
               delta0=(1-eps)*newDelta + eps*rep(1/N,N),
               Par0=newPar0,
               fixPar=fixPar,
               knownStates=knownStates,
               lambda=lamb,
               userBounds=userBounds,
               workBounds=workBounds,
               nlmPar = list('iterlim'=iterlim,
                             'print.level'=0))

# Save hmm
saveRDS(hmm0,
        make_title(paste0(directory,"/params/"),
                   paste0(model,"-",
                          lamb,"-",
                          holdout_whale,"-",
                          rand_seed,"-",
                          "hmm.rds")))

}
