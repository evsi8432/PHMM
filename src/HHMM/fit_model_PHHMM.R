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

### change to hierarchical ###
Data0 <- NULL
for(ID in unique(Data$ID)){

  boutnum <- 1
  whale_data <- Data[Data$ID %in% ID,]
  start <- whale_data$stime[1]
  end <- whale_data$etime[nrow(whale_data)]
  time <- start

  while(time < end){

    # add levels 1 and 2i
    tmp1 <- as.data.frame(matrix(nrow=2,ncol=ncol(Data)+1))
    colnames(tmp1) <- c("boutnum",colnames(Data))
    tmp1$ID <- ID
    tmp1$level <- c("1","2i")

    # add level 2
    tmp2 <- whale_data[(whale_data$stime >= time) & (whale_data$stime < (time + span*60)),]

    rest <- any(c(1) %in% tmp2$knownState)
    trav <- any(c(2) %in% tmp2$knownState)
    forg <- any(c(3) %in% tmp2$knownState)

    # check for repeated knownState
    if(forg){
      tmp2$knownState <- 3
    }
    else if (rest + trav > 1){
      tmp2$knownState <- 4
    }

    if(nrow(tmp2) > 0){

      # combine it all together
      tmp2$level <- "2"

      tmp1$boutnum <- boutnum
      tmp2$boutnum <- boutnum
      boutnum <- boutnum + 1

      Data0 <- rbind(Data0,tmp1,tmp2)
    }

    # move the current time
    if(nrow(tmp2) > 0){
      time <- tail(tmp2,n=1)[1,"etime"]
    } else {
      time <- time + span*60
    }
  }
}
Data <- Data0
rownames(Data) <- 1:nrow(Data)
Data$label <- factor(Data$label,levels = 1:7)

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
    DM[[feature]] <- matrix(cbind(kronecker(c(1,1,1,0,0,0),diag(N_working)),
                                  kronecker(c(0,0,0,1,1,1),diag(N_working))),
                            nrow=N*2,ncol=N_working*2)
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

if(model == "no"){
  for(i in 1:N_working){
    for(feature in features1){
      if(dist[[feature]] == "norm"){
        Par0[[feature]][i,1] <- mean(Data[,feature],na.rm=T) + rnorm(1)*sd(Data[,feature],na.rm=T) # mean
        Par0[[feature]][i,2] <- 0.5*rnorm(1) + log(sd(Data[,feature],na.rm=T)) - log(N0)
      } else if (dist[[feature]] == "mvnorm2") {
        featurex <- paste0(feature,".x")
        featurey <- paste0(feature,".y")
        Par0[[feature]][i,1] <- mean(Data[,featurex],na.rm=T) + rnorm(1)*sd(Data[,featurex],na.rm=T) # mean
        Par0[[feature]][i,2] <- mean(Data[,featurey],na.rm=T) + rnorm(1)*sd(Data[,featurey],na.rm=T) # mean
        Par0[[feature]][i,3] <- 0.5*rnorm(1) + log(sd(Data[,featurex],na.rm=T)) - log(N0) # sd x
        Par0[[feature]][i,4] <- 0.0
        Par0[[feature]][i,5] <- 0.5*rnorm(1) + log(sd(Data[,featurey],na.rm=T)) - log(N0) # sd y
      } else if (substring(dist[[feature]], 1,3) == "cat") {
        ncats <- as.integer(substring(dist[[feature]], 4))
        for(j in 1:(ncats-1)){
          Par0[[feature]][i,j] <- log(mean(Data[,feature] == j,na.rm=T)) - log(mean(Data[,feature] == ncats,na.rm=T))
        }
        Par0[[feature]][i,] = Par0[[feature]][i,] + rnorm(ncats-1)
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

### get initial TPMs from previous fits ###
hierBeta <- data.tree::Node$new("HHMM beta")

if(model == "no"){
  hierBeta$AddChild(name="level1",
                    beta=matrix(rnorm(3*(3-1),mean=-2,sd=1),nrow=1))

  hierBeta$AddChild(name="level2")
  hierBeta$level2$AddChild(name="rest",beta=matrix(rnorm(N0*(N0-1),mean=0,sd=1),nrow=1))
  hierBeta$level2$AddChild(name="trav",beta=matrix(rnorm(N0*(N0-1),mean=0,sd=1),nrow=1))
  hierBeta$level2$AddChild(name="forg",beta=matrix(rnorm(N0*(N0-1),mean=0,sd=1),nrow=1))
} else {
  hierBeta <- getPar0(old_hmm)$hierBeta
}

### get initial TPMs from previous fits ###
hierDelta <- data.tree::Node$new("HHMM delta")
if(model == "no"){
  hierDelta$AddChild(name="level1",
                     beta=matrix(rep(1/3,3),nrow=1))

  hierDelta$AddChild(name="level2")
  hierDelta$level2$AddChild(name="rest",beta=matrix(rep(1/N0,N0),nrow=1))
  hierDelta$level2$AddChild(name="trav",beta=matrix(rep(1/N0,N0),nrow=1))
  hierDelta$level2$AddChild(name="forg",beta=matrix(rep(1/N0,N0),nrow=1))
} else {
  hierDelta <- getPar0(old_hmm)$hierDelta
}

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

### Check Parameters ###
checkPar0(data=Data,
          hierStates=hierStates,
          hierDist=hierDist,
          hierBeta=hierBeta,
          hierDelta=hierDelta,
          Par0=Par0,
          fixPar=fixPar,
          userBounds=userBounds,
          workBounds=workBounds,
          DM=DM)

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
               nlmPar = list('iterlim'=1000))

print(hmm)

# get labels switching values
Data <- cbind(Data,stateProbs(hmm))

accs <- matrix(nrow = 3, ncol = 3)
behaviours <- c("resting","travelling","foraging")
states <- c(paste("rest",1:statesPerBehaviour[1]),
            paste("trav",1:statesPerBehaviour[2]),
            paste("forg",1:statesPerBehaviour[3]))

# make accuracy matrix
for(i in 1:3){
  for(j in 1:3){
    behaviour_data <- Data[Data$knownState %in% i,]
    accs[i,j] <- sum(behaviour_data[,states[(N0*(j-1)+1):(N0*j)]])
  }
}

# get pairs
if (model %in% "no"){
  hs <- HungarianSolver(-accs)
  small_pairs <- hs$pairs

  # expand pairs to all states
  hmm$pairs <- matrix(nrow = N, ncol = 2)
  hmm$pairs[,1] <- seq(1,N)
  for (i in 1:3){
    hmm$pairs[(N0*(i-1)+1):(N0*i),2] <- (N0*(small_pairs[i,2]-1)+1):(N0*small_pairs[i,2])
  }
} else {
  hmm$pairs <- old_hmm$pairs
}

print(hmm$pairs)

# Save hmm
saveRDS(hmm,
        make_title(paste0(directory,"/params/"),
                   paste0(model,"-",
                          holdout_whale,"-",
                          rand_seed,"-",
                          "hmm.rds")))
}