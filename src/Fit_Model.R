library(momentuHMM)
library(circular)
library(CircStats)
library(runner)
library(ggplot2)
library(tidyr)
library(dplyr)
set.seed(0)

### Parse Command line arguments ###
args = commandArgs(trailingOnly=TRUE)
#args <- 5

### Load Data ###
setwd("~/Documents/Research/PHMM/src")
whale = "I107" #"I145" "D26"

### start by saying how many states we will have ###
N_deep <- 2
N_shallow <- 2

shallow_thresh <- 5
deep_thresh <- 10

N_fine <- 2

N_coarse <- N_deep + N_shallow
N_total <- N_coarse * N_fine


Data <- read.csv(paste0('../Dat/diary/',whale,'/processed_Data_',whale,'.csv'),
                 colClasses=c("level"="character"))

### preprare data for HHMM ###
Data$stime <- as.POSIXct(Data$stime, origin = '1970-01-01')
Data$etime <- as.POSIXct(Data$etime, origin = '1970-01-01')

date0 = format(Data$Time[1], "%Y-%m-%d")

Data[is.na(Data$head.tm1),"head.tm1"] <- 0
Data[is.na(Data$pitch.tm1),"pitch.tm1"] <- 0
Data[is.na(Data$roll.tm1),"roll.tm1"] <- 0

Data[is.na(Data$aw1.tm1),"aw1.tm1"] <- 0
Data[is.na(Data$aw2.tm1),"aw2.tm1"] <- 0
Data[is.na(Data$aw3.tm1),"aw3.tm1"] <- 0

Data$w.cov <- Data$w
Data[is.na(Data$w.cov),"w.cov"] <- min(Data$w.cov,na.rm=TRUE)

Data$VeDBA.cov <- Data$VeDBA
Data[is.na(Data$VeDBA.cov),"VeDBA.cov"] <- min(Data$VeDBA.cov,na.rm=TRUE)

#summary(Data)

### Label Deep vs Shallow Dives ###
Data$broadDiveType[is.na(Data$maxDepth)] <- NA  # initialize column
Data$broadDiveType[!is.na(Data$maxDepth)] <- 3  # unknown
Data$broadDiveType[Data$maxDepth > deep_thresh] <- 2     # deep
Data$broadDiveType[Data$maxDepth < shallow_thresh] <- 1      # shallow
Data$broadDiveType <- factor(Data$broadDiveType,levels = 1:3)

### label the data ###
Data$labelnum <- 6
Data$labelnum[Data$label == "clicks"] <- 1
Data$labelnum[Data$label == "click train"] <- 2
Data$labelnum[Data$label == "buzz"] <- 3
Data$labelnum[Data$label == "punch"] <- 4
Data$labelnum[Data$label == "crunch"] <- 5
Data$labelnum <- factor(Data$labelnum,levels = 1:6)

# prep the data idk why
Data <- prepData(Data,
                 coordNames=NULL,
                 hierLevels=c("1","2i","2"))

Data$labelnum[Data$level == "1"] <- NA
Data$labelnum[Data$level == "2i"] <- NA

# Start with a Coarse-Scale model only



# get states
hierStates <- data.tree::Node$new("Orca HHMM states")
state_num <- 1
for (n in 1:N_coarse){
  name_coarse = paste0("D",n)
  hierStates$AddChild(name = name_coarse) 
  for (n_star in 1:N_fine){
    name_fine = paste0("S",n,n_star)
    hierStates[[name_coarse]]$AddChild(name_fine, state=state_num)
    state_num <- state_num + 1
  }
}


### coarse scale Distributions ###

hierDist <- data.tree::Node$new("orca HHMM dist")
hierDist$AddChild("level1")
hierDist$level1$AddChild("diveDuration", dist="gamma")
hierDist$level1$AddChild("maxDepth", dist="gamma")
hierDist$level1$AddChild("broadDiveType", dist="cat3")


### fine scale Distributions ###

hierDist$AddChild("level2")

#hierDist$level2$AddChild("hv", dist="norm")
#hierDist$level2$AddChild("jp", dist="norm")
#hierDist$level2$AddChild("rajp", dist="vm")

hierDist$level2$AddChild("w", dist="norm")
hierDist$level2$AddChild("aw1", dist="norm")
hierDist$level2$AddChild("aw2", dist="norm")
hierDist$level2$AddChild("aw3", dist="norm")
hierDist$level2$AddChild("labelnum", dist="cat6")

# get formula using depth as covariate

# hierFormula <- data.tree::Node$new("orca HHMM formula")
# hierFormula$AddChild(name="level1", formula=~1)
# hierFormula$AddChild(name="level2", formula=~ad)
# plot(hierFormula)


### Initalize Gamma ###

# initialize beta (initial ptm)
hierBeta <- data.tree::Node$new("orca beta")
hierBeta$AddChild("level1",
                  beta=matrix(rep(-1.00,N_coarse*(N_coarse-1)),
                              nrow=1,
                              ncol=N_coarse*(N_coarse-1)))
hierBeta$AddChild("level2")

for(n_coarse in 1:N_coarse){
  print(n_coarse)
  hierBeta$level2$AddChild(paste0("D",n_coarse),
                           beta=matrix(rep(-1.00,N_fine*(N_fine-1)),
                                       nrow=1,
                                       ncol=N_fine*(N_fine-1),
                                       byrow=TRUE)) 
}
#plot(hierBeta)

# initialize delta (initial distribution)
hierDelta <- data.tree::Node$new("orca delta")
hierDelta$AddChild("level1",delta=matrix(rep(0.0,N_coarse-1),nrow=1))
hierDelta$AddChild("level2")
for (n_coarse in 1:N_coarse){
  hierDelta$level2$AddChild(paste0("D",n_coarse),
                            delta=matrix(rep(0.0,N_fine-1),nrow=1))
}
#plot(hierDelta)

### make Design Matrix ###

DM <- list()
userBounds <- list()

dd_names <- list(paste0(rep(c("mean_","sd_"),each=N_total),1:N_total),
                 c(paste0(rep(c("mean_","sd_"),each=N_coarse),1:N_coarse,":(Intercept)")))
md_names <- list(paste0(rep(c("mean_","sd_"),each=N_total),1:N_total),
                 c(paste0(rep(c("mean_","sd_"),each=N_coarse),1:N_coarse,":(Intercept)")))
dt_names <- list(paste0(rep(c("prob1_","prob2_"),each=N_total),1:N_total),
                 c(paste0(rep(c("prob1_","prob2_"),each=N_coarse),1:N_coarse,":(Intercept)")))

DM$diveDuration <- matrix(kronecker(diag(2*N_coarse),rep(1,N_fine)),
                          nrow=2*N_total,
                          ncol=2*N_coarse,
                          dimnames=dd_names)

DM$maxDepth <- matrix(kronecker(diag(2*N_coarse),rep(1,N_fine)),
                      nrow=2*N_total,
                      ncol=2*N_coarse,
                      dimnames=md_names)

DM$broadDiveType <- matrix(kronecker(diag(2*N_coarse),rep(1,N_fine)),
                           nrow=2*N_total,
                           ncol=2*N_coarse,
                           dimnames=dt_names)

# userBounds$diveDuration <- matrix(c(rep(c(0.0,Inf),N_total),
#                                     rep(c(1.0,Inf),N_total)),
#                                   nrow=N_total*2,byrow=TRUE,
#                                   dimnames=list(rownames(DM$diveDuration),c("lower","upper")))
# 
# userBounds$maxDepth <- matrix(c(rep(c(0.0,Inf),N_total),
#                                 rep(c(0.1,Inf),N_total)),
#                               nrow=N_total*2,byrow=TRUE,
#                               dimnames=list(rownames(DM$maxDepth),c("lower","upper")))

# fine-scale design matrix
DM$w = list(mean = ~1, sd = ~1)
DM$aw1 <- list(mean = ~1 + aw1.tm1, sd = ~1)
DM$aw2 <- list(mean = ~1 + aw2.tm1, sd = ~1)
DM$aw3 <- list(mean = ~1 + aw3.tm1, sd = ~1)

DM$labelnum <- list(prob1 = ~ 1,
                    prob2 = ~ 1,
                    prob3 = ~ 1,
                    prob4 = ~ 1,
                    prob5 = ~ 1 + w.cov)


### Get the old HHMM ###
old_hhmm <- readRDS("../params/hhmm_Mar_1.rds")
Par <- getPar0(model = old_hhmm)$Par
hierBeta <- getPar0(model = old_hhmm)$hierBeta
hierDelta <- getPar0(model = old_hhmm)$hierDelta

### Set initial parameters ###

# Coarse-scale 
#Par$maxDepth[3] <- log(15) # mean of maxDepth, Dive type 3
#Par$maxDepth[7] <- log(7) # sd of maxDepth, Dive type 3

#Par$diveDuration[3] <- log(93) # mean of diveDuration, Dive type 3
#Par$diveDuration[7] <- log(51) # sd of diveDuration, Dive type 3

#Par$maxDepth[4] <- log(100) # mean of maxDepth, Dive type 4
#Par$maxDepth[8] <- log(40) # sd of maxDepth, Dive type 4

#Par$diveDuration[4] <- log(280) # mean of diveDuration, Dive type 4
#Par$diveDuration[8] <- log(116) # sd of diveDuration, Dive type 4

  
ParBroadDiveType0 <- c(
  rep(0.00, N_shallow), rep(-100,N_deep), # prob shallow, each state (intercept * 8)
  rep(-100, N_shallow), rep(0.00,N_deep)  # prob deep, each state (intercept * 8)
)

# fine-scale
Parlabelnum0 <- c(
  rep(-4.0, 8),                     # prob clicks, each state (intercept * 8)
  rep(-4.0, 8),                     # prob click train, each state (intercept * 8)
  rep(-10,6),rep(-4.0,2),           # prob buzz, each state (intercept * 8)
  rep(-10,6),rep(-4.0,2),           # prob punch, each state (intercept * 8)
  -10.0, 0.0,                       # prob crunch, state 1 (S11) (intercept,w)
  -10.0, 0.0,                       # prob crunch, state 2 (S12) (intercept,w)
  -10.0, 0.0,                       # prob crunch, state 3 (S21) (intercept,w)
  -10.0, 0.0,                       # prob crunch, state 4 (S22) (intercept,w)
  -10.0, 0.0,                       # prob crunch, state 5 (S31) (intercept,w)
  -10.0, 0.0,                       # prob crunch, state 6 (S32) (intercept,w)
  -10.0, 1.0,                       # prob crunch, state 7 (S41) (intercept,w)
  -10.0, 1.0)                       # prob crunch, state 8 (S42) (intercept,w)

#Par$labelnum = Parlabelnum0
#Par$broadDiveType = ParBroadDiveType0

### Fix Known Parameters ###

fixPar = list()

fixParlabelnum0 <- c(
  rep(NA, 8),           # prob clicks, each state (intercept * 8)
  rep(NA, 8),           # prob click train, each state (intercept * 8)
  rep(-10,6),rep(NA,2), # prob buzz, each state
  rep(-10,6),rep(NA,2), # prob punch, each state
  -10, 0,               # prob crunch, state 1 (S11) (intercept,w)
  -10, 0,               # prob crunch, state 2 (S12) (intercept,w)
  -10, 0,               # prob crunch, state 3 (S21) (intercept,w)
  -10, 0,               # prob crunch, state 4 (S22) (intercept,w)
  -10, 0,               # prob crunch, state 5 (S31) (intercept,w)
  -10, 0,               # prob crunch, state 6 (S32) (intercept,w)
  NA, NA,               # prob crunch, state 7 (S41) (intercept,w)
  NA, NA)               # prob crunch, state 8 (S42) (intercept,w)

fixParBroadDiveType0 <- c(
  rep(0.00, N_shallow), rep(-100,N_deep), # prob shallow, each state (intercept * 8)
  rep(-100, N_shallow), rep(0.00,N_deep)  # prob deep, each state (intercept * 8)
)

#fixPar$labelnum = fixParlabelnum0
fixPar$broadDiveType <- fixParBroadDiveType0
fixPar$labelnum <- fixParlabelnum0

### Check Parameters ###
checkPar0(Data[1:10,],
         hierStates = hierStates,
         hierDist = hierDist,
         hierBeta = hierBeta,
         hierDelta = hierDelta,
         fixPar=fixPar,
         DM=DM,
         Par0=Par)


### Fit Model ###

maxiters <- as.numeric(args[1])

ptm <- proc.time()

if (maxiters == 0 | is.na(maxiters)){
  hhmm <- fitHMM(Data,
                 hierStates = hierStates,
                 hierDist = hierDist,
                 hierBeta = hierBeta,
                 hierDelta = hierDelta,
                 DM=DM,
                 Par0=Par,
                 fixPar = fixPar,
                 nlmPar = list('print.level'=2))
} else {
  print(paste("max iters:",maxiters))
  hhmm <- fitHMM(Data,
                 hierStates = hierStates,
                 hierDist = hierDist,
                 hierBeta = hierBeta,
                 hierDelta = hierDelta,
                 DM=DM,
                 Par0=Par,
                 fixPar = fixPar,
                 nlmPar = list('print.level'=2,
                               'iterlim'=maxiters))
}

print(proc.time()-ptm)

hhmm

saveRDS(hhmm,"../params/hhmm_Mar_1.rds")
