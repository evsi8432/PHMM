library(momentuHMM)
library(circular)
library(CircStats)
library(runner)
library(ggplot2)
library(tidyr)
library(dplyr)
library(mclust)
library(mixreg)
library(readxl)
library(diveMove)
library(boot)
library(data.table)
library(party)
library(numDeriv)
library(data.tree)

### Parse Command line arguments ###
args = commandArgs(trailingOnly=TRUE)
args <- c("Female",1,"None")
#args <- 5

### Select Sex ###
sex <- args[1]
#sex = "Male"
#sex = "Female"

### set random seed ###
rand_seed <- args[2] 
#rand_seed <- 1
set.seed(rand_seed)

### hold out any whales? ###
holdout_whale <- args[3]
#holdout_whale <- "None"
#holdout_whale <- "I145"

### use other model as beginning values? (cross-validation) ###
#reference_model <- args[4]
if(sex == "Male"){
  if(holdout_whale == "None"){
    reference_model <- "../params/2023-04-03/hierHMM_None_Female_1.rds"
  } else {
    reference_model <- "../params/2023-04-03/hierHMM_None_Male_1.rds"
  }
} else {
  if(holdout_whale == "None"){
    reference_model <- "../params/2023-01-25/hierHmm_2023-01-25_13-16-15.rds" 
  } else {
    reference_model <- "../params/2023-04-03/hierHMM_None_Female_1.rds"
  }
}

# select number of states
N_coarse <- 3
N_fine <- 4
N <- N_coarse*N_fine

### Load Data ###
setwd("~/Documents/Research/PHMM/src")

# load data
Data <- data.frame(fread('../../dat/Final_Data_Beth.csv'))
ethogram <- data.frame(fread('../../dat/Final_ethogram_Beth.csv'))

# get post dive int
#Data[is.infinite(Data$postDiveInt),"postDiveInt"] <- NA

# turn data positive for gamma dists
Data$postDiveInt <- exp(Data$postDiveInt) + runif(nrow(Data))/2
Data$maxDepth <- exp(Data$maxDepth)
Data$diveDuration <- exp(Data$diveDuration) + runif(nrow(Data))/2

# change post-diveinterval to a categorical
Data$postDiveIntCat <- 0
Data$postDiveIntCat <- Data$postDiveIntCat + (Data$postDiveInt < 10)
Data$postDiveIntCat <- Data$postDiveIntCat + (Data$postDiveInt < 5)
Data$postDiveIntCat <- Data$postDiveIntCat + (Data$postDiveInt < 2)
Data$postDiveIntCat <- factor(Data$postDiveIntCat)

# make wiggliness a log-scale
Data$avg_w_low <- log(Data$avg_w_low)
Data$avg_w_high <- log(Data$avg_w_high)

# turn dive Type into a factor
Data$diveType <- factor(Data$diveType,
                        levels = c("Resting", "Travelling", "Foraging", "Logging"))
Data$diveType <- addNA(Data$diveType)
Data[!is.na(Data$postDiveInt) & 
       Data$postDiveInt >= 10,"diveType"] <- "Logging"

### Label Deep vs Shallow Dives ###
shallow_thresh <- c(0,7.5)
medium_thresh <- c(10,30)
deep_thresh <- c(50,Inf)

Data$broadDiveType <- NA
Data$broadDiveType[(shallow_thresh[1] < Data$maxDepth) & (Data$maxDepth < shallow_thresh[2])] <- 1
Data$broadDiveType[(medium_thresh[1] < Data$maxDepth) & (Data$maxDepth < medium_thresh[2])] <- 2
Data$broadDiveType[(deep_thresh[1] < Data$maxDepth) & (Data$maxDepth < deep_thresh[2])] <- 3
Data$broadDiveType[Data$diveType %in% "Logging"] <- 4
Data$broadDiveType[is.na(Data$broadDiveType)] <- 5

### add label to next deep dive for foraging dives

# deep dives mean foraging
Data[Data$broadDiveType %in% 3,"diveType"] <- "Foraging"

# deep-ish dives after labels mean foraging too
for(stime in Data$stime[Data$diveType %in% "Foraging"]){
  future_dives <- Data[(Data$stime >= stime) & (Data$stime < stime + 120),]
  future_deep_dives <- future_dives[future_dives$maxDepth > 30,]
  if(nrow(future_deep_dives) >= 1){
    dive_ind <- rownames(future_deep_dives)[1]
    Data[dive_ind,"diveType"] <- "Foraging"
  }
}

### Label knownStates ###
Data$knownState <- NA
Data$knownState[Data$diveType %in% "Resting"] <- 1
Data$knownState[Data$diveType %in% "Travelling"] <- 2
Data$knownState[Data$diveType %in% "Foraging"] <- 3
Data$knownState[is.na(Data$knownState)] <- 4

span <- 10 # minutes

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
    time <- tail(tmp2,n=1)[1,"etime"]
  }
}
Data <- Data0

rownames(Data) <- 1:nrow(Data)
  
### add sex to Data ###
Data$Sex <- "Female"
Data[Data0$ID %in% c("I107","D21","L87","L88"),"Sex"]  <- "Male"
Data$Sex <- factor(Data$Sex)

### split data into test and train ###
Data_all_labeled <- Data
Data_all_unlabeled <- Data
Data_all_unlabeled[(Data_all_unlabeled$ID %in% holdout_whale) & (Data$level %in% 2),"knownState"] <- 4
Data <- Data[!(Data$ID %in% holdout_whale),]

### Select features ###
features <- c("maxDepth","diveDuration","postDiveInt","knownState")

### Set distributions ###
hierDist <- data.tree::Node$new("Killer Whale HHMM dist")
hierDist$AddChild(name="level1")
hierDist$AddChild(name="level2")
for(feature in features){
  hierDist$level2$AddChild(name="maxDepth", dist="gamma")
  hierDist$level2$AddChild(name="diveDuration", dist="gamma")
  hierDist$level2$AddChild(name="postDiveInt", dist="gamma")
  hierDist$level2$AddChild(name="broadDiveType", dist="cat5")
  hierDist$level2$AddChild(name="knownState", dist="cat4")
}
plot(hierDist)

### Define States ###
hierStates <- data.tree::Node$new("Killer Whale HHMM states")
hierStates$AddChild(name="resting")
hierStates$resting$AddChild(name="srest", state=1) # shallow resting
hierStates$resting$AddChild(name="mrest", state=2) # medium resting
hierStates$resting$AddChild(name="drest", state=3) # deep resting
hierStates$resting$AddChild(name="lrest", state=4) # logging
hierStates$AddChild(name="travelling")
hierStates$travelling$AddChild(name="strav", state=5) # shallow travel
hierStates$travelling$AddChild(name="mtrav", state=6) # medium travel
hierStates$travelling$AddChild(name="dtrav", state=7) # deep travel
hierStates$travelling$AddChild(name="ltrav", state=8) # logging
hierStates$AddChild(name="foraging")
hierStates$foraging$AddChild(name="sforg", state=9)  # shallow foraging
hierStates$foraging$AddChild(name="mforg", state=10) # medium foraging
hierStates$foraging$AddChild(name="dforg", state=11) # deep foraging
hierStates$foraging$AddChild(name="lforg", state=12) # logging
plot(hierStates)

### Set initial parameters ###
#reference_model <- "../params/2023-01-25/hierHmm_2023-01-25_00-24-15.rds"
prev_hmm <- readRDS(reference_model)

get_norm_par <- function(Data,feature){
  
  Par0 <- rep(0,2*N_fine)
  
  #Par0[1] <- log(mean(Data[Data$broadDiveType %in% 1,feature],na.rm=TRUE))
  Par0[1] <- log(sample(Data[(Data$broadDiveType %in% 1) & (!is.na(Data[,feature])),
                             feature],1))
  Par0[N_fine+1] <- log(sd(Data[Data$broadDiveType %in% 1,feature],
                               na.rm=TRUE) + 0.1*sd(Data[,feature],na.rm=TRUE))
  
  #Par0[2] <- log(mean(Data[Data$broadDiveType %in% 2,feature],na.rm=TRUE))
  Par0[2] <- log(sample(Data[(Data$broadDiveType %in% 2) & (!is.na(Data[,feature])),
                             feature],1))
  Par0[N_fine+2] <- log(sd(Data[Data$broadDiveType %in% 2,feature],
                           na.rm=TRUE) + 0.1*sd(Data[,feature],na.rm=TRUE))
  
  #Par0[3] <- log(mean(Data[Data$broadDiveType %in% 3,feature],na.rm=TRUE))
  Par0[3] <- log(sample(Data[(Data$broadDiveType %in% 3) & (!is.na(Data[,feature])),
                         feature],1))
  Par0[N_fine+3] <- log(sd(Data[Data$broadDiveType %in% 3,feature],
                           na.rm=TRUE) + 0.1*sd(Data[,feature],na.rm=TRUE))
  
  #Par0[4] <- log(mean(Data[Data$broadDiveType %in% 4,feature],na.rm=TRUE))
  Par0[4] <- log(sample(Data[(Data$broadDiveType %in% 4) & (!is.na(Data[,feature])),
                         feature],1))
  Par0[N_fine+4] <- log(sd(Data[Data$broadDiveType %in% 4,feature],
                           na.rm=TRUE) + 0.1*sd(Data[,feature],na.rm=TRUE))
  
  Par0 <- getPar0(prev_hmm)$Par[[feature]]
  
  # add random noise
  if (rand_seed != 1){
    Par0[1:N_fine] <- rnorm(N_fine,
                            mean = Par0[1:N_fine],
                            sd = 0.1)
    
    Par0[(N_fine+1):(2*N_fine)] <- rnorm(N_fine,
                                         mean = Par0[(N_fine+1):(2*N_fine)],
                                         sd = 0.25) 
  }
  
  return(Par0)
}

ParMaxDepth0 <- get_norm_par(Data,"maxDepth")
ParDiveDuration0 <- get_norm_par(Data,"diveDuration")
ParPostDiveInt0 <- get_norm_par(Data,"postDiveInt")

eps <- 1e-15
ParKnownState0 <- c(
  0.1,0.1,0.1,0.1,eps,eps,eps,eps,eps,eps,eps,eps, # probability of label 1 (Resting), states 1-12
  eps,eps,eps,eps,0.1,0.1,0.1,0.1,eps,eps,eps,eps, # probability of label 2 (Travelling), states 1-12
  eps,eps,eps,eps,eps,eps,eps,eps,0.1,0.1,0.1,0.1  # probability of label 3 (Foraging), states 1-12
)
ParBroadDiveType0 <- c(
  0.1,eps,eps,eps,0.1,eps,eps,eps,0.1,eps,eps,eps, # probability of label 1 (shallow), states 1-12
  eps,0.1,eps,eps,eps,0.1,eps,eps,eps,0.1,eps,eps, # probability of label 2 (medium), states 1-12
  eps,eps,0.1,eps,eps,eps,0.1,eps,eps,eps,0.1,eps, # probability of label 3 (deep), states 1-12
  eps,eps,eps,1.0-4*eps,eps,eps,eps,1.0-4*eps,eps,eps,eps,1.0-4*eps  # probability of label 4 (logging), states 1-12
)

Par0 = list(maxDepth = ParMaxDepth0,
            diveDuration = ParDiveDuration0,
            postDiveInt = ParPostDiveInt0,
            knownState = ParKnownState0,
            broadDiveType = ParBroadDiveType0)

### get initial PTMs from previous fits ###
hierBeta <- data.tree::Node$new("HHMM beta")

coarse_betas <- unique(prev_hmm$mle$beta[1,])[-1]
if(rand_seed != 1){
  coarse_betas <- rnorm(n = length(coarse_betas),
                        mean = coarse_betas,
                        sd = 1)
}

hierBeta$AddChild(name="level1",beta=matrix(coarse_betas,1))
hierBeta$AddChild(name="level2")

fine_betas <- unique(prev_hmm$mle$beta[3,])[-N_fine]
if(rand_seed != 1){
  fine_betas <- rnorm(n = length(coarse_betas),
                      mean = coarse_betas,
                      sd = 1)
}

fine_betas[c(2,5,8,11,12+2,12+5,12+8,12+11)] <- -100
hierBeta$level2$AddChild(name="resting",beta=matrix(fine_betas[1:12],1))
hierBeta$level2$AddChild(name="travelling",beta=matrix(fine_betas[13:24],1))
hierBeta$level2$AddChild(name="foraging",beta=matrix(fine_betas[25:36],1))

### Fix Known Parameters ###
fixPar = list()

fixHierBeta <- data.tree::Node$new("HHMM beta")
fixHierBeta$AddChild(name="level1",beta=matrix(rep(NA,N_coarse*(N_coarse-1)),1))
fixHierBeta$AddChild(name="level2")
no_deep <- matrix(c(   NA,-1e10,NA,
                       NA,-1e10,NA,
                       NA,-1e10,NA,
                       NA,-1e10,NA),
                         1)

fixHierBeta$level2$AddChild(name="resting",beta=no_deep)
fixHierBeta$level2$AddChild(name="travelling",beta=no_deep)
fixHierBeta$level2$AddChild(name="foraging",beta=matrix(rep(NA,N_fine*(N_fine-1)),1))

fixPar$knownState <- Par0$knownState
fixPar$broadDiveType <- Par0$broadDiveType
fixPar$beta <- fixHierBeta
  
### make coarse-scale shared ###
DM0 <- matrix(cbind(kronecker(c(1,1,1,0,0,0),diag(4)),
                    kronecker(c(0,0,0,1,1,1),diag(4))),
              nrow=12*2,
              ncol=8,
              dimnames=list(c(paste0("mean_",1:12),
                              paste0("sd_",1:12)),
                            paste0(rep(c("mean","sd"),each=4),
                                   c("_1-5-9:(Intercept)",
                                     "_2-6-10:(Intercept)",
                                     "_3-7-11:(Intercept)",
                                     "_4-8-12:(Intercept)"))))

DM <- list(maxDepth=DM0,
           diveDuration=DM0,
           postDiveInt=DM0)

# prep data
Data <- prepData(Data,
                 coordNames=NULL,
                 hierLevels = c("1", "2i", "2"))

Data_all_labeled <- prepData(Data_all_labeled,
                             coordNames=NULL,
                             hierLevels = c("1", "2i", "2"))

Data_all_unlabeled <- prepData(Data_all_unlabeled,
                               coordNames=NULL,
                               hierLevels = c("1", "2i", "2"))

### Check Parameters ###
if(FALSE){#(sex == "Male"){
  checkPar0(data=Data[Data$Sex == sex,],
            hierStates=hierStates,
            hierDist=hierDist,
            hierBeta=hierBeta,
            Par0=Par0,
            fixPar=fixPar,
            DM=DM)
  
  hhmm <- fitHMM(data=Data[Data$Sex == sex,],
                 hierStates=hierStates,
                 hierDist=hierDist,
                 hierBeta=hierBeta,
                 Par0=Par0,
                 DM=DM,
                 fixPar=fixPar,
                 nlmPar = list('print.level'=2,
                               'iterlim'=10000))
} else if(FALSE){#(sex == "Female") {
  checkPar0(data=Data[Data$Sex == sex,],
            hierStates=hierStates,
            hierDist=hierDist,
            hierBeta=hierBeta,
            Par0=Par0,
            fixPar=fixPar,
            DM=DM)
  
  hhmm <- fitHMM(data=Data[Data$Sex == sex,],
                 hierStates=hierStates,
                 hierDist=hierDist,
                 hierBeta=hierBeta,
                 Par0=Par0,
                 DM=DM,
                 fixPar=fixPar,
                 nlmPar = list('print.level'=2,
                               'iterlim'=10000))
}

#hhmm

# Save HHMM
#date = Sys.Date()
date = "2023-04-03"

time = Sys.time()
time <- gsub(" ","_",time)
time <- gsub(":","-",time)

dir.create(file.path("../params/",date))
dir.create(file.path("../plots/",date))

#saveRDS(hhmm,paste0("../params/",date,"/hierHmm_", 
#                    holdout_whale,"_",
#                    sex,"_",
#                    rand_seed,".rds"))

hhmm <- readRDS(paste0("../params/",date,"/hierHmm_", 
                       holdout_whale,"_",
                       sex,"_",
                       rand_seed,".rds"))

# Do some small EDA
if (sex == "Male"){
  #hhmm <- readRDS("../params/2023-01-25/hierHmm_2023-01-25_00-24-15.rds")
  Data <- Data[Data$Sex %in% "Male",]
  Data_all_unlabeled <- Data_all_unlabeled[Data_all_unlabeled$Sex %in% "Male",]
  Data_all_labeled <- Data_all_labeled[Data_all_labeled$Sex %in% "Male",]

} else {
  #hhmm <- readRDS("../params/2023-01-24/hierHmm_2023-01-24_23-45-39.rds")
  Data <- Data[Data$Sex %in% "Female",]
  Data_all_unlabeled <- Data_all_unlabeled[Data_all_unlabeled$Sex %in% "Female",]
  Data_all_labeled <- Data_all_labeled[Data_all_labeled$Sex %in% "Female",]
}

# define dive types
group.names <- c("1" = "Resting (Shallow)",
                 "2" = "Resting (Medium)",
                 "3" = "Resting (Deep)",
                 "4" = "Resting (Logging)",
                 "5" = "Travelling (Shallow)",
                 "6" = "Travelling (Medium)",
                 "7" = "Travelling (Deep)",
                 "8" = "Travelling (Logging)",
                 "9" = "Foraging (Shallow)",
                 "10" = "Foraging (Medium)",
                 "11" = "Foraging (Deep)",
                 "12" = "Foraging (Logging)")

group.names.coarse <- c("1" = "Resting",
                        "2" = "Travelling",
                        "3" = "Foraging",
                        "4" = "Logging")

group.names.fine <- c("1" = "Shallow",
                      "2" = "Medium",
                      "3" = "Deep",
                      "4" = "Logging")

# remake hhmm with labels removed
Par0 <- getPar0(hhmm)

hhmm_unlabeled <- fitHMM(data=Data_all_unlabeled,
                         hierStates=hierStates,
                         hierDist=hierDist,
                         hierBeta=Par0$hierBeta,
                         hierDelta=Par0$hierDelta,
                         Par0=Par0$Par,
                         DM=DM,
                         fixPar=fixPar,
                         nlmPar = list('print.level'=2,
                                       'stepmax'=1e-100,
                                       'iterlim'=1))

# decode vitirbi 
vstates <- viterbi(hhmm_unlabeled)

vstates_fine <- rep(NA,length(vstates))
vstates_fine[vstates %in% c(1,5,9)] <- 1
vstates_fine[vstates %in% c(2,6,10)] <- 2
vstates_fine[vstates %in% c(3,7,11)] <- 3
vstates_fine[vstates %in% c(4,8,12)] <- 4

vstates_coarse <- rep(NA,length(vstates))
vstates_coarse[vstates %in% c(1,2,3)] <- 1
vstates_coarse[vstates %in% c(5,6,7)] <- 2
vstates_coarse[vstates %in% c(9,10,11)] <- 3
vstates_coarse[vstates %in% c(4,8,12)] <- 4

fb <- stateProbs(hhmm_unlabeled)
colnames(fb) <- group.names

# add dive types to unlabelled and labelled Data
Data_all_unlabeled$vstate1 <- vstates
Data_all_unlabeled$vstate1_coarse <- vstates_coarse
Data_all_unlabeled$vstate1_fine <- vstates_fine

Data_all_unlabeled <- left_join(Data_all_unlabeled,
                                data.frame(vstate1=1:N,vdivetype=group.names),
                                by="vstate1")
Data_all_unlabeled <- cbind(Data_all_unlabeled,fb)

# add dive types to labelled Data
Data_all_labeled$vstate1 <- vstates
Data_all_labeled$vstate1_coarse <- vstates_coarse
Data_all_labeled$vstate1_fine <- vstates_fine

Data_all_labeled <- left_join(Data_all_labeled,
                              data.frame(vstate1=1:N,vdivetype=group.names),
                              by="vstate1")
Data_all_labeled <- cbind(Data_all_labeled,fb)

### run data on test whales ###
conf_matrix <- matrix(rep(0,N_coarse*N_coarse),
                      nrow=N_coarse,
                      ncol=N_coarse)

behaviours <- c("Resting","Travelling","Foraging")
for(i in 1:N_coarse){
  dives_true_behaviour <- Data_all_labeled[(Data_all_labeled$knownState %in% i) & 
                                           (Data_all_labeled$diveType %in% behaviours[i]) &
                                           (Data_all_labeled$ID %in% holdout_whale),]  
  
  dives_true_behaviour <- dives_true_behaviour %>% distinct(boutnum,.keep_all=TRUE)
  for(j in 1:N_coarse){
    conf_matrix[i,j] = sum(dives_true_behaviour$vstate1_coarse == j)
  }
}

rownames(conf_matrix) <- c("True Resting", "True Travelling", "True Foraging")
colnames(conf_matrix) <- c("Predicted Resting", "Predicted Travelling", "Predicted Foraging")

write.csv(conf_matrix,paste0("../params/",date,"/conf_matrix_",
                             holdout_whale,"_",
                             sex,"_",
                             rand_seed,".csv"))
