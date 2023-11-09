library(momentuHMM)
library(dplyr)
library(mclust)
library(data.table)
library(RcppHungarian)
library(mvtnorm)

#setwd("/Users/evsi8432/Documents/Research/PHMM/src/bash")

# get command-line arguments
args <- commandArgs(trailingOnly=TRUE)
#args <- c("hier_logMDDD_logWTotal_2-2-2_dd-02_2023-09-06.R",0)
#args <- c("hier_logMDDD_2-2-2_dd-02_2023-08-30.R",0)

opt_file <- args[1]
args <- as.integer(args[2])

# get options
source(paste0('../opt/',opt_file))

# find N for working scale
if(share_fine){
  N_working <- N0
} else {
  N_working <- N
}

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
  args_list <- sind:(120*n_retries-1)
} else {
  args_list <- c(args)
}

for(args in args_list){

# Set Seed
rand_seed <- (args[1] %% n_retries) + 1
set.seed(rand_seed)

# select holdout whale
whale_ind <- (floor(args[1] / n_retries) %% 20) + 1
whales <- c("none","A100a","A100b","A113a","A113b",
            "D21a","D21b","D26a","D26b",
            "I107a","I107b","I129","I145a","I145b",
            "L87","L88","R48a","R48b","R58a","R58b")
holdout_whale <- whales[whale_ind]

# Select Model
model_ind <- floor(args[1] / (20*n_retries)) + 1
models <- c("no","fixed_1","fixed_2","half_random","random","random_2")
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
  print(paste(file,"already exists. continuing..."))
  next
}

# get data
source("../HMM/load_data.R")

# get rid of held out whale
Data <- Data[!(Data$ID %in% holdout_whale),]

### Set distributions ###
if(model != "no"){
  dist$label = paste0("cat",nlevels(Data$label))
}

hierDist <- data.tree::Node$new("Killer Whale HHMM dist")
hierDist$AddChild(name="level1")
hierDist$AddChild(name="level2")
for(feature in names(dist)){
  hierDist$level2$AddChild(name=feature, dist=dist[feature])
}

### Define States ###
bhavs <- c("rest","trav","forg")
behaviours <- c("resting","travelling","foraging")

inds <- c(0,cumsum(statesPerBehaviour))
b_inds <- list()
b_inds$resting <- (inds[1]+1):inds[2]
b_inds$travelling <- (inds[2]+1):inds[3]
b_inds$foraging <- (inds[3]+1):inds[4]

hierStates <- data.tree::Node$new("Killer Whale HHMM states")
statenum <- 1
for(i in 1:3){
  hierStates$AddChild(name=bhavs[i])
  for(j in 1:statesPerBehaviour[i]){
    hierStates[[bhavs[i]]]$AddChild(name=paste(bhavs[i],j),
                                    state=statenum)
    statenum <- statenum + 1
  }
}

### set up DM matrix ###
DM <- list()

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
  if(dist[feature] == "mvnorm2"){
    if(share_fine){
      DM[[feature]] <- matrix(cbind(kronecker(c(1,1,1, 0,0,0, 0,0,0, 0,0,0,       0,0,0),diag(N0)), # mean x
                                    kronecker(c(0,0,0, 1,1,1, 0,0,0, 0,0,0,       0,0,0),diag(N0)), # mean y
                                    kronecker(c(0,0,0, 0,0,0, 1,1,1, 0.5,0.5,0.5, 0,0,0),diag(N0)), # sig x
                                    kronecker(c(0,0,0, 0,0,0, 0,0,0, 1,1,1,       0,0,0),diag(N0)), # sig xy
                                    kronecker(c(0,0,0, 0,0,0, 0,0,0, 0.5,0.5,0.5, 1,1,1),diag(N0))), # sig y
                              nrow=N*5,ncol=N0*5) 
    } else {
      DM[[feature]] <- matrix(cbind(kronecker(c(1, 0, 0, 0,   0),diag(N)), # mean x
                                    kronecker(c(0, 1, 0, 0,   0),diag(N)), # mean y
                                    kronecker(c(0, 0, 1, 0.5, 0),diag(N)), # sig x
                                    kronecker(c(0, 0, 0, 1,   0),diag(N)), # sig xy
                                    kronecker(c(0, 0, 0, 0.5, 1),diag(N))), # sig y
                              nrow=N*5,ncol=N*5)       
    }
  } else if(dist[feature] == "norm"){
    if(share_fine){
      DM[[feature]] <- matrix(cbind(kronecker(c(1,1,1,0,0,0),diag(N_working)),
                                    kronecker(c(0,0,0,1,1,1),diag(N_working))),
                              nrow=N*2,ncol=N_working*2)
    } else {
      DM[[feature]] <- matrix(cbind(kronecker(c(1,0),diag(N_working)),
                                    kronecker(c(0,1),diag(N_working))),
                              nrow=N*2,ncol=N_working*2)
    }

  } else if (substring(dist[[feature]], 1,3) == "cat"){
    ncats <- as.integer(substring(dist[[feature]], 4))
    if(share_fine){
      DM[[feature]] <- kronecker(kronecker(diag(ncats-1),c(1,1,1)),diag(N0))
    } else {
      DM[[feature]] <- matrix(diag(N*(ncats-1)))
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
    Par0[[feature]] <- matrix(rep(NA,2*N_working), nrow = N_working)
  } else if(dist[[feature]] == "mvnorm2"){
    Par0[[feature]] <- matrix(rep(NA,5*N_working), nrow = N_working)
  } else if (substring(dist[[feature]], 1,3) == "cat") {
    ncats <- as.integer(substring(dist[[feature]], 4))
    Par0[[feature]] <- matrix(rep(NA,(ncats-1)*N_working), nrow = N_working)
  } else {
    print("feature distribtuion not recognized")
  }
}

for(i in 1:N_working){
  for(feature in features1){
    if(dist[[feature]] == "norm"){
      Par0[[feature]][i,1] <- mean(Data[,feature],na.rm=T) + 0.1*rnorm(1)*sd(Data[,feature],na.rm=T) # mean
      Par0[[feature]][i,2] <- 0.1*rnorm(1) + log(sd(Data[,feature],na.rm=T)) - 0.1*log(N0)
    } else if (dist[[feature]] == "mvnorm2") {
      featurex <- paste0(feature,".x")
      featurey <- paste0(feature,".y")
      Par0[[feature]][i,1] <- mean(Data[,featurex],na.rm=T) + 0.1*rnorm(1)*sd(Data[,featurex],na.rm=T) # mean
      Par0[[feature]][i,2] <- mean(Data[,featurey],na.rm=T) + 0.1*rnorm(1)*sd(Data[,featurey],na.rm=T) # mean
      Par0[[feature]][i,3] <- 0.1*rnorm(1) + log(sd(Data[,featurex],na.rm=T)) - 0.1*log(N0) # sd x
      Par0[[feature]][i,4] <- 0.0
      Par0[[feature]][i,5] <- 0.1*rnorm(1) + log(sd(Data[,featurey],na.rm=T)) - 0.1*log(N0) # sd y
    } else if (substring(dist[[feature]], 1,3) == "cat") {
      ncats <- as.integer(substring(dist[[feature]], 4))
      for(j in 1:(ncats-1)){
        Par0[[feature]][i,j] <- log(mean(Data[,feature] == j,na.rm=T)) - log(mean(Data[,feature] == ncats,na.rm=T))
      }
      Par0[[feature]][i,] = Par0[[feature]][i,] + rnorm(ncats-1)
    }
  }
}

# get accuracy of mixture model for each behaviour and state
if (!share_fine){
  accs <- matrix(rep(0,N^2),nrow=N,ncol=N)
  for(i in 1:N){
    
    # get data associated with the behaviour in question
    behaviour_ind <- which(i <= cumsum(statesPerBehaviour))[1]
    behaviour_data <- Data[Data$knownState %in% behaviour_ind & !(Data$pseudolabel),]
    
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
}


if(model %in% c("fixed_1","fixed_2")){
  Par0$label = fixPar$label
} else if(model == "half_random"){
  Par0$label <- c(qlogis(0.1),qlogis(0.1),qlogis(0.001),-100)
} else if(model == "random"){
  Par0$label <- c(qlogis(0.1),qlogis(0.1),qlogis(0.001),qlogis(0.001))
}

### get initial TPMs from previous fits ###
hierBeta <- data.tree::Node$new("HHMM beta")
hierBeta$AddChild(name="level1",
                  beta=matrix(rnorm(3*(3-1),mean=-2,sd=1),nrow=1))

hierBeta$AddChild(name="level2")
hierBeta$level2$AddChild(name="rest",beta=matrix(rnorm(N0*(N0-1),mean=0,sd=1),nrow=1))
hierBeta$level2$AddChild(name="trav",beta=matrix(rnorm(N0*(N0-1),mean=0,sd=1),nrow=1))
hierBeta$level2$AddChild(name="forg",beta=matrix(rnorm(N0*(N0-1),mean=0,sd=1),nrow=1))

### get initial TPMs from previous fits ###
hierDelta <- data.tree::Node$new("HHMM delta")
hierDelta$AddChild(name="level1",
                   beta=matrix(rep(1/3,3),nrow=1))

hierDelta$AddChild(name="level2")
hierDelta$level2$AddChild(name="rest",beta=matrix(rep(1/N0,N0),nrow=1))
hierDelta$level2$AddChild(name="trav",beta=matrix(rep(1/N0,N0),nrow=1))
hierDelta$level2$AddChild(name="forg",beta=matrix(rep(1/N0,N0),nrow=1))

# update user and work bounds
for(feature in features1){
  if(dist[feature] == "mvnorm2"){
    userBounds[[feature]] <- matrix(c(rep(c(-Inf,-Inf,0.01,0.01,0.01),each=N),
                                      rep(c(Inf,Inf,Inf,Inf,Inf),each=N)),
                                    nrow=5*N,ncol=2)
    workBounds[[feature]] <- matrix(c(rep(c(-Inf,-Inf,-Inf,-Inf,-Inf),each=N_working),
                                      rep(c(Inf,Inf,Inf,-0.01,Inf),each=N_working)),
                                    nrow=5*N_working,ncol=2)
  } else if (dist[feature] == "norm"){
    userBounds[[feature]] <- matrix(c(rep(c(-Inf,0.01),each=N),
                                      rep(c(Inf,Inf),each=N)),
                                    nrow=2*N,ncol=2)
    workBounds[[feature]] <- matrix(c(rep(c(-Inf,-Inf),each=N_working),
                                      rep(c(Inf,Inf),each=N_working)),
                                    nrow=2*N_working,ncol=2)
  } else if (substring(dist[[feature]], 1,3) == "cat") {
    ncats <- as.integer(substring(dist[[feature]], 4))
    
    userBounds[[feature]] <- matrix(c(rep(0.00,(ncats-1)*N),
                                      rep(1.00,each=(ncats-1)*N)),
                                    nrow=(ncats-1)*N,ncol=2)
    workBounds[[feature]] <- matrix(c(rep(-Inf,each=(ncats-1)*N_working),
                                      rep(Inf,each=(ncats-1)*N_working)),
                                    nrow=2*N_working,ncol=2)
  }
}


# prep data
Data <- prepData(Data,
                 coordNames=NULL,
                 hierLevels = c("1", "2i", "2"))

Data[Data$level != 2,features2] <- NA

### fit model ###
hmm <- fitHMM(data=Data,
               hierStates=hierStates,
               hierDist=hierDist,
               hierBeta=hierBeta,
               hierDelta=hierDelta,
               Par0=Par0,
               fixPar=fixPar,
               userBounds=userBounds,
               workBounds=workBounds,
               DM=DM,
               nlmPar = list('iterlim'=iterlim))

# find best pairings with old HMM
probs <- stateProbs(hmm)
hmm$data <- cbind(hmm$data,probs)

# make accuracy matrix
accs <- matrix(nrow = 3, ncol = 3)
behaviours <- c("resting","travelling","foraging")
states <- c(paste("rest",1:statesPerBehaviour[1]),
            paste("trav",1:statesPerBehaviour[2]),
            paste("forg",1:statesPerBehaviour[3]))
inds <- c(0,cumsum(statesPerBehaviour))

for(i in 1:3){
  for(j in 1:3){
    behaviour_data <- hmm$data[hmm$data$knownState %in% i,]
    behaviour_data <- behaviour_data[!is.na(behaviour_data$boutnum),] %>% 
      distinct(boutnum,.keep_all=TRUE)
    accs[i,j] <- sum(behaviour_data[,states[(inds[j]+1):inds[j+1]]])
  }
}

# get pairs
hs <- HungarianSolver(-accs)
small_pairs <- hs$pairs

# expand pairs to all states
pairs <- matrix(nrow = N, ncol = 2)
pairs[,1] <- seq(1,N)
for (i in 1:3){
  pairs[(N0*(i-1)+1):(N0*i),2] <- (N0*(small_pairs[i,2]-1)+1):(N0*small_pairs[i,2])
}

# fix Par0
Par0 <- getPar0(hmm)
newPar0 <- list()
for(feature in features1){
  if(dist[feature] == "norm"){
    oldPar0 <- matrix(Par0$Par[[feature]],nrow=N)
    newPar0[[feature]] <- matrix(rep(NA,2*N), nrow = N)
    for(i in 1:N){
      newPar0[[feature]][pairs[i,1],] <- oldPar0[pairs[i,2],]
    }
  } else if (dist[feature] == "mvnorm2") {
    oldPar0 <- matrix(Par0$Par[[feature]],nrow=N)
    newPar0[[feature]] <- matrix(rep(NA,5*N), nrow=N)
    for(i in 1:N){
      newPar0[[feature]][pairs[i,1],] <- oldPar0[pairs[i,2],]
    }    
  } else if (substring(dist[[feature]], 1,3) == "cat") {
    ncats <- as.integer(substring(dist[[feature]], 4))
    oldPar0 <- matrix(Par0$Par[[feature]],nrow=N)
    newPar0[[feature]] <- matrix(rep(NA,(ncats-1)*N), nrow=N)
    for(i in 1:N){
      newPar0[[feature]][pairs[i,1],] <- oldPar0[pairs[i,2],]
    } 
  }
}

if(model != "no"){
  newPar0$label <- Par0$Par$label
}

# fix delta
oldDelta <- c(0,Par0$hierDelta$level1$delta)
newDelta <- oldDelta[small_pairs[,2]]
newDelta <- newDelta - newDelta[1]
newDelta <- matrix(newDelta[-1],nrow=1)
hierDelta <- data.tree::Node$new("HHMM delta")
hierDelta$AddChild(name="level1",delta=newDelta)

hierDelta$AddChild(name="level2")
hierDelta$level2$AddChild(name="rest",delta=Par0$hierDelta$level2[[bhavs[small_pairs[1,2]]]]$delta)
hierDelta$level2$AddChild(name="trav",delta=Par0$hierDelta$level2[[bhavs[small_pairs[2,2]]]]$delta)
hierDelta$level2$AddChild(name="forg",delta=Par0$hierDelta$level2[[bhavs[small_pairs[3,2]]]]$delta)

# fix beta
oldBeta <- Par0$hierBeta$level1$beta
newBeta <- c()
for(i in 1:(3-1)){
  newBeta <- c(newBeta,0,oldBeta[(3*(i-1)+1):(3*i)])
}
newBeta <- c(newBeta,0)
newBeta <- matrix(newBeta,nrow=3,byrow=T)
newBeta <- newBeta[small_pairs[,2],small_pairs[,2]]
newBeta0 <- c()
for(i in 1:3){
  newBeta0 <- c(newBeta0,newBeta[i,-i])
}
newBeta <- matrix(newBeta0,nrow=1)

hierBeta <- data.tree::Node$new("HHMM beta")
hierBeta$AddChild(name="level1",beta=newBeta)

hierBeta$AddChild(name="level2")
hierBeta$level2$AddChild(name="rest",beta=Par0$hierBeta$level2[[bhavs[small_pairs[1,2]]]]$beta)
hierBeta$level2$AddChild(name="trav",beta=Par0$hierBeta$level2[[bhavs[small_pairs[2,2]]]]$beta)
hierBeta$level2$AddChild(name="forg",beta=Par0$hierBeta$level2[[bhavs[small_pairs[3,2]]]]$beta)

# refit the model
hmm0 <- fitHMM(data=Data,
              hierStates=hierStates,
              hierDist=hierDist,
              hierBeta=hierBeta,
              hierDelta=hierDelta,
              Par0=newPar0,
              fixPar=fixPar,
              userBounds=userBounds,
              workBounds=workBounds,
              DM=DM,
              nlmPar = list('iterlim'=iterlim))

# Save hmm
saveRDS(hmm,
        make_title(paste0(directory,"/params/"),
                   paste0(model,"-",
                          holdout_whale,"-",
                          rand_seed,"-",
                          "hmm.rds")))
}