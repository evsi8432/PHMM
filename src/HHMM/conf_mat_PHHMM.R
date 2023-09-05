library(momentuHMM)
library(dplyr)
library(mclust)
library(data.table)

#setwd("/Users/evsi8432/Documents/Research/PHMM/src")

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
  args_list <- sind:54
} else {
  args_list <- c(args)
}

for(args in args_list){

# set seed
set.seed(1)

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

# get data
source("../HMM/load_data.R")

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
    print(file)
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

# prep data
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
