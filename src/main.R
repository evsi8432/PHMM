library(momentuHMM)
library(circular)
library(CircStats)
library(runner)
library(ggplot2)
library(tidyr)
library(dplyr)

setwd("~/Documents/Research/PHMM/src")
whale = "I107" # I145
date0 = "2020-08-25" # "2020-08-30

# start by saying how many states we will have
N_coarse <- 3
N_fine <- 3

# This script loads in the data, makes a "RawData" df, and adds the important
# stuff to a data frame called "Data"
#source("Prep_Data.R")
Data <- read.csv(paste0('../Dat/diary/processed_Data_',whale,'.csv'),
                 colClasses=c("level"="character"))

# preprare data for HHMM 
Data$stime <- as.POSIXct(Data$stime, origin = '1970-01-01')
Data$etime <- as.POSIXct(Data$etime, origin = '1970-01-01')
Data$label_num[!is.na(Data$label)] <- N_coarse*N_fine
Data$label_num <- as.numeric(Data$label_num)
Data <- prepData(Data,coordNames=NULL,
                 hierLevels=c("1","2i","2"))

summary(Data)

rawData <- read.csv(paste0('../Dat/diary/processed_rawData_',whale,'.csv'),
                    colClasses=c("Time"="double"))
rawData$Time <- as.POSIXct(rawData$Time, origin = '1970-01-01')
rawData <- rawData[,!(names(rawData) %in% c("vstate1","vstate2","label"))]

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
plot(hierStates)

# get distributions
hierDist <- data.tree::Node$new("orca HHMM dist")
hierDist$AddChild("level1")
hierDist$level1$AddChild("diveDuration", dist="gamma")
hierDist$level1$AddChild("maxDepth", dist="gamma")
hierDist$AddChild("level2")
hierDist$level2$AddChild("hv", dist="norm")
hierDist$level2$AddChild("jp", dist="norm")
hierDist$level2$AddChild("rajp", dist="vm")
#hierDist$level2$AddChild("w", dist="gamma")
plot(hierDist)

# get formula using depth as covariate in 
#hierFormula <- data.tree::Node$new("orca HHMM formula")
#hierFormula$AddChild(name="level1", formula=~1)
#hierFormula$AddChild(name="level2", formula=~ad)
#plot(hierFormula)

# initialize beta (initial ptm)
hierBeta <- data.tree::Node$new("orca beta")
hierBeta$AddChild("level1",
                  beta=matrix(rep(-1.00,N_coarse*(N_coarse-1)),
                              nrow=1,
                              ncol=N_coarse*(N_coarse-1)))
hierBeta$AddChild("level2")

for (n_coarse in 1:N_coarse){
  print(n_coarse)
  hierBeta$level2$AddChild(paste0("D",n_coarse),
                           beta=matrix(rep(-1.00,N_fine*(N_fine-1)),
                                       nrow=1,
                                       ncol=N_fine*(N_fine-1),
                                       byrow=TRUE)) 
}

# initialize delta (initial distribution)
hierDelta <- data.tree::Node$new("orca delta")
hierDelta$AddChild("level1",delta=matrix(rep(0.0,N_coarse-1),nrow=1))
hierDelta$AddChild("level2")
for (n_coarse in 1:N_coarse){
  hierDelta$level2$AddChild(paste0("D",n_coarse),
                            delta=matrix(rep(0.0,N_fine-1),nrow=1))
}

# constrain coarse scale states by making the coarse-scale 
# parameters the same for all fine-scale states
N_total <- N_coarse*N_fine

DM <- list()
userBounds <- list()

A <- diag(2*N_coarse)
for (row in 1:N_coarse){
  for (col in row:N_coarse){
    A[row,col] <- 1
  }
}

DM$diveDuration <- matrix(kronecker(diag(2*N_coarse),rep(1,N_fine)),
                          nrow=2*N_total,
                          ncol=2*N_coarse,
                          dimnames=list(paste0(rep(c("mean_","sd_"),each=N_total),
                                               1:N_total),
                                        c(paste0(rep(c("mean_","sd_"),each=N_coarse),
                                                 1:N_coarse,
                                                 ":(Intercept)"))))

DM$maxDepth <- matrix(kronecker(diag(2*N_coarse),rep(1,N_fine)),
                      nrow=2*N_total,
                      ncol=2*N_coarse,
                      dimnames=list(paste0(rep(c("mean_","sd_"),each=N_total),
                                           1:N_total),
                                    c(paste0(rep(c("mean_","sd_"),each=N_coarse),
                                             1:N_coarse,
                                             ":(Intercept)"))))

# userBounds$diveDuration <- matrix(c(rep(c(0.0,Inf),N_total),
#                                     rep(c(1.0,Inf),N_total)),
#                                   nrow=N_total*2,byrow=TRUE,
#                                   dimnames=list(rownames(DM$diveDuration),c("lower","upper")))
# 
# userBounds$maxDepth <- matrix(c(rep(c(0.0,Inf),N_total),
#                                 rep(c(0.1,Inf),N_total)),
#                               nrow=N_total*2,byrow=TRUE,
#                               dimnames=list(rownames(DM$maxDepth),c("lower","upper")))
# 
# # constrain fine-scale states to be the same
# A <- diag(N_fine)
# for (row in 1:N_fine){
#   for (col in row:N_fine){
#     A[row,col] <- 1
#   }
# }
# 
# DM$jp <- matrix(cbind(kronecker(c(rep(1,N_coarse),rep(0,N_coarse)),diag(N_fine)),
#                       kronecker(c(rep(0,N_coarse),rep(1,N_coarse)),diag(N_fine))),
#                 nrow=2*N_total,ncol=2*N_fine,
#                 dimnames=list(c(paste0("mean_",1:N_total),
#                                 paste0("sd_",1:N_total)),
#                               paste0(rep(c("mean","sd"),each=N_fine),
#                               paste0("_",1:N_fine,":(Intercept)"))))
# 
# userBounds$jp <- matrix(c(rep(c(-Inf,Inf),N_total),
#                           rep(c(0.1,Inf),N_total)),
#                         nrow=N_total*2,byrow=TRUE,
#                         dimnames=list(rownames(DM$jp),c("lower","upper")))
# 
# DM$hv <- matrix(cbind(kronecker(c(rep(1,N_coarse),rep(0,N_coarse)),diag(N_fine)),
#                       kronecker(c(rep(0,N_coarse),rep(1,N_coarse)),diag(N_fine))),
#                 nrow=2*N_total,ncol=2*N_fine,
#                 dimnames=list(c(paste0("mean_",1:N_total),
#                                 paste0("sd_",1:N_total)),
#                               paste0(rep(c("mean","sd"),each=N_fine),
#                                      paste0("_",1:N_fine,":(Intercept)"))))
# 
# userBounds$hv <- matrix(c(rep(c(-Inf,Inf),N_total),
#                           rep(c(0.1,Inf),N_total)),
#                         nrow=N_total*2,byrow=TRUE,
#                         dimnames=list(rownames(DM$hv),c("lower","upper")))
# 
# DM$rajp <- matrix(cbind(kronecker(c(rep(1,N_coarse),rep(0,N_coarse)),diag(N_fine)),
#                 kronecker(c(rep(0,N_coarse),rep(1,N_coarse)),diag(N_fine))),
#                 nrow=2*N_total,ncol=2*N_fine,
#                 dimnames=list(c(paste0("mean_",1:N_total),
#                                 paste0("concentration_",1:N_total)),
#                               paste0(rep(c("mean","sd"),each=N_fine),
#                                      paste0("_",1:N_fine,":(Intercept)"))))
# 
# userBounds$rajp <- matrix(c(rep(c(0.5,Inf),N_total)),
#                           nrow=N_total,byrow=TRUE,
#                           dimnames=list(rownames(DM$rajp[(N_total+1):(2*N_total),]),
#                                         c("lower","upper")))

get_par_coarse0 <- function(feature_data){
  mus <- quantile(feature_data,
                  probs = seq(0,1,1/(N_coarse+1)), 
                  na.rm=T)[2:(N_coarse+1)]
  
  sigs <- rep(sd(feature_data,na.rm=TRUE),N_coarse)
  
  pars <- c(rep(mus,each=N_fine),
            rep(sigs,each=N_fine))
  
  return(pars)
}

ParDiveDuration0 <- get_par_coarse0(Data$diveDuration)
ParMaxDepth0 <- get_par_coarse0(Data$maxDepth)

get_par_fine0 <- function(feature_data){
  mus <- seq(min(feature_data,na.rm=TRUE),
             max(feature_data,na.rm=TRUE),
             length.out = N_fine+2)[1:N_fine+1]
  
  sigs <- rep(sd(feature_data,na.rm=TRUE),N_fine)
  
  pars <- c(rep(mus,times=N_coarse),
            rep(sigs,times=N_coarse))
  
  return(pars)
}

Parhv0 <- get_par_fine0(Data$hv)
Parjp0 <- get_par_fine0(Data$jp)
Parrajp0 <- get_par_fine0(Data$rajp)

# get initial parameters
Par0 = list(diveDuration = ParDiveDuration0,
            maxDepth = ParMaxDepth0,
            hv = Parhv0,
            jp = Parjp0,
            rajp = Parrajp0)

Par <- getParDM(Data,
                hierStates=hierStates,
                hierDist=hierDist,
                Par=Par0,
                DM=DM,
                estAngleMean = list(rajp=TRUE))

# check hierarchical model specifications
checkPar0(Data,
          hierStates = hierStates,
          hierDist = hierDist,
          hierBeta = hierBeta, 
          hierDelta = hierDelta,
          DM=DM,
          Par0=Par,
          estAngleMean = list(rajp=TRUE))

# fit the HHMM
hhmm <- fitHMM(Data,
               hierStates = hierStates,
               hierDist = hierDist,
               hierBeta = hierBeta, 
               hierDelta = hierDelta,
               DM=DM,
               Par0=Par,
               estAngleMean = list(rajp=TRUE))#,
               #knownStates = Data$label)


labels <- NULL
labels$level1 <- Data$label[Data$level == 1]
labels$level2 <- Data$label[Data$level != 1]
vstates <- viterbi(hhmm, hierarchical = TRUE)


dive_types <- data.frame(vstate1 = vstates$level1,
                         label1 = labels$level1,
                         divenum = seq(1,length(vstates$level1)))
if (!("vstate1" %in% colnames(rawData))){
  rawData <- left_join(rawData,dive_types,by="divenum")
}

subdive_states <- data.frame(vstate2 = vstates$level2,
                             label2 = labels$level2,
                             divesegnum = seq(1,length(vstates$level2)))
rawData$divesegnum = rawData$segnum + rawData$divenum
if (!("vstate2" %in% colnames(rawData))){
  rawData <- left_join(rawData,subdive_states,by="divesegnum")
}
rawData$vstate2[rawData$segnum == 0] <- NA
rawData$label2[rawData$segnum == 0] <- NA

# now do the same for Data
Data$vstate1 <- NA
Data$vstate2 <- NA
Data$vstate1[Data$level == 1] <- vstates$level1
Data$vstate2[Data$level != 1] <- vstates$level2

Data$vstate1[Data$vstate2 %in% c("S11","S12","S13")] <- "D1"
Data$vstate1[Data$vstate2 %in% c("S21","S22","S23")] <- "D2"
Data$vstate1[Data$vstate2 %in% c("S31","S32","S33")] <- "D3"


# add elevation
rawData$Elevation <- -rawData$p

# prepare the data for plotting
rawDataDownsampled <- rawData[seq(1,nrow(rawData),50),]
rawDataDownsampled$jp <- log(rawDataDownsampled$jp)
rawDataDownsampled$hv <- logit(rawDataDownsampled$hv)

rawDataDownLong <- rawDataDownsampled %>% 
  pivot_longer(cols = c("Elevation","jp","rajp","hv","Aw_1","Aw_2","Aw_3"), 
               names_to = "feature")

# plot the results
ggplot(rawDataDownLong[rawDataDownLong$divenum %in% 150:200,], 
       aes(x=Time, y=value, group=divenum)) +
  geom_line(aes(color=interaction(vstate1,is.na(label1)))) + 
  facet_wrap(~feature,ncol=1,scales='free_y')

ggplot(rawDataDownLong[(rawDataDownLong$divenum %in% 325:350) &
                       (rawDataDownLong$segnum != 0),], 
       aes(x=Time, y=value, group=segnum)) +
  geom_line(aes(color=vstate2),size=1) + 
  facet_wrap(~feature,ncol=1,scales='free_y')

# plot each feature individually
plot_coarse_feature <- function(feature){
  shape1 = hhmm$mle[[feature]]['mean','S11']^2 / hhmm$mle[[feature]]['sd','S11']^2
  rate1 = hhmm$mle[[feature]]['mean','S11'] / hhmm$mle[[feature]]['sd','S11']^2
  shape2 = hhmm$mle[[feature]]['mean','S21']^2 / hhmm$mle[[feature]]['sd','S21']^2
  rate2 = hhmm$mle[[feature]]['mean','S21'] / hhmm$mle[[feature]]['sd','S21']^2
  shape3 = hhmm$mle[[feature]]['mean','S31']^2 / hhmm$mle[[feature]]['sd','S31']^2
  rate3 = hhmm$mle[[feature]]['mean','S31'] / hhmm$mle[[feature]]['sd','S31']^2
  
  plot <- ggplot(Data,aes_string(x=feature)) +
    geom_histogram(aes(fill=vstate1,y=stat(density)),alpha=0.5,position = "identity") + 
    stat_function(fun = function(x) {dgamma(x,shape=shape1,rate=rate1)},
                  aes(color = 'D1')) +
    stat_function(fun = function(x) {dgamma(x,shape=shape2,rate=rate2)},
                  aes(color = 'D2')) +
    stat_function(fun = function(x) {dgamma(x,shape=shape3,rate=rate3)},
                  aes(color = 'D3'))
  
  return(plot)
}

print(plot_coarse_feature('diveDuration'))
print(plot_coarse_feature('maxDepth'))

# plot each feature individually
plot_fine_feature_norm <- function(feature){
  
  plt <- ggplot(Data,aes_string(x=feature)) +
    geom_histogram(aes(fill=vstate2,y=stat(density)),alpha=0.5,position = "identity")
  
  states <- list()
  
  #for (state in c('S11','S12','S21','S22','S31','S32')){
  #  mu = hhmm$mle[[feature]]['mean',state]
  #  sig = hhmm$mle[[feature]]['sd',state]
  #  
  #  plt <- plt + stat_function(fun = dnorm,
  #                             args = list(mean=mu,sd=sig),
  #                             aes(color = state))
  #}
  
  plt <- plt + facet_wrap(~vstate1,ncol=1,scales='free_y')
  
  return(plt)
}

print(plot_fine_feature_norm('hv'))
print(plot_fine_feature_norm('jp'))

plot_fine_feature_vm <- function(feature){

  plt <- ggplot(Data,aes_string(x=feature)) +
    geom_histogram(aes(fill=vstate2,y=stat(density)),alpha=0.5,position = "identity")

  states <- list()

  for (state in c('S11','S12','S21','S22','S31','S32')){
    mu = hhmm$mle[[feature]]['mean',state]
    kappa = hhmm$mle[[feature]]['concentration',state]

    plt <- plt + stat_function(fun = dvonmises,
                               args = list(mu=mu,kappa=kappa),
                               aes(color = state))
  }

  plt <- plt + facet_wrap(~vstate1,ncol=1,scales='free_y')

  return(plt)
}

print(plot_fine_feature_vm('rajp'))

