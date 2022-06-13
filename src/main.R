library(momentuHMM)
library(circular)
library(CircStats)
library(runner)
library(ggplot2)
library(tidyr)
library(dplyr)
library(mclust)
library(mixreg)
set.seed(0)

setwd("~/Documents/Research/PHMM/src")
whale = "I107" # I145
date0 = "2020-08-25" # "2020-08-30

# start by saying how many states we will have
N_coarse <- 4
N_fine <- 2
N_total <- N_coarse*N_fine

# This script loads in the data, makes a "RawData" df, and adds the important
# stuff to a data frame called "Data"
#source("Prep_Data.R")
Data <- read.csv(paste0('../Dat/diary/',whale,'/processed_Data_',whale,'.csv'),
                 colClasses=c("level"="character"))

# prepare data for HHMM 
Data$stime <- as.POSIXct(Data$stime, origin = '1970-01-01')
Data$etime <- as.POSIXct(Data$etime, origin = '1970-01-01')
Data[is.na(Data$aw1.tm1),"aw1.tm1"] <- 0
Data[is.na(Data$aw2.tm1),"aw2.tm1"] <- 0
Data[is.na(Data$aw3.tm1),"aw3.tm1"] <- 0

Data$diveDuration <- log(Data$diveDuration)
Data$maxDepth <- log(Data$maxDepth)
Data$label_num <- NA


# label the data
Data$label_num[Data$label == 'crunch'] <- N_coarse*N_fine
Data$label_num <- as.numeric(Data$label_num)
Data <- prepData(Data,coordNames=NULL,
                 hierLevels=c("1","2i","2"))

summary(Data)

rawData <- read.csv(paste0('../Dat/diary/',whale,'/processed_rawData_',whale,'.csv'),
                    colClasses=c("Time"="double"))
rawData$Time <- as.POSIXct(rawData$Time, origin = '1970-01-01')
rawData <- rawData[,!(names(rawData) %in% c("vstate1","vstate2","label"))]

# do some exploratory data analysis
Data$d_aw1 <- c(0,diff(Data$aw1))
Data$d_aw2 <- c(0,diff(Data$aw2))
Data$d_aw3 <- c(0,diff(Data$aw3))

ggplot(Data) +
  geom_histogram(aes(x=jp,y=..density..)) +
  facet_wrap(~label,ncol=1,scales='free_y')

ggplot(Data) +
  geom_histogram(aes(x=rajp,y=..density..)) +
  facet_wrap(~label,ncol=1,scales='free_y')

ggplot(Data) +
  geom_histogram(aes(x=hv,y=..density..)) +
  facet_wrap(~label,ncol=1,scales='free_y')

ggplot(Data) +
  geom_histogram(aes(x=log(htv),y=..density..)) +
  facet_wrap(~label,ncol=1,scales='free_y')

ggplot(Data) +
  geom_histogram(aes(x=aw1,y=..density..)) +
  facet_wrap(~label,ncol=1,scales='free_y')

ggplot(Data) +
  geom_histogram(aes(x=aw2,y=..density..)) +
  facet_wrap(~label,ncol=1,scales='free_y')

ggplot(Data) +
  geom_histogram(aes(x=aw3,y=..density..)) +
  facet_wrap(~label,ncol=1,scales='free_y')

ggplot(Data) +
  geom_histogram(aes(x=d_aw1,y=..density..)) +
  facet_wrap(~label,ncol=1,scales='free_y')

ggplot(Data) +
  geom_histogram(aes(x=d_aw2,y=..density..)) +
  facet_wrap(~label,ncol=1,scales='free_y')

ggplot(Data) +
  geom_histogram(aes(x=d_aw3,y=..density..)) +
  facet_wrap(~label,ncol=1,scales='free_y')

ggplot(Data) +
  geom_histogram(aes(x=w,y=..density..)) +
  facet_wrap(~label,ncol=1,scales='free_y')

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
hierDist$level1$AddChild("diveDuration", dist="norm")
hierDist$level1$AddChild("maxDepth", dist="norm")
hierDist$AddChild("level2")
#hierDist$level2$AddChild("hv", dist="norm")
#hierDist$level2$AddChild("jp", dist="norm")
#hierDist$level2$AddChild("rajp", dist="vm")
#hierDist$level2$AddChild("w", dist="norm")
#hierDist$level2$AddChild("aw1", dist="norm")
#hierDist$level2$AddChild("aw2", dist="norm")
#hierDist$level2$AddChild("aw3", dist="norm")
plot(hierDist)

# get formula using depth as covariate in 
#hierFormula <- data.tree::Node$new("orca HHMM formula")
#hierFormula$AddChild(name="level1", formula=~1)
#hierFormula$AddChild(name="level2", formula=~ad)
#plot(hierFormula)

# initialize beta (initial ptm)
hierBeta <- data.tree::Node$new("orca beta")
hierBeta$AddChild("level1",
                  beta=matrix(rep(0.00,N_coarse*(N_coarse-1)),
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

DM <- list()
userBounds <- list()

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

# DM$w = list(mean = ~1, sd = ~1)

get_par_coarse0 <- function(feature_data,feature){
  
  # do kmeans clustering
  clust <- kmeans(feature_data,centers = matrix(c(1,3,4.5,3,4,5.5),nrow=3,ncol=2))
  
  mus <- rep(NA,N_coarse)
  sigs <- rep(NA,N_coarse)
  
  for (i in 1:N_coarse){
    mus[i]  <- mean(feature_data[clust$cluster == i,feature])
    sigs[i] <- log(sd  (feature_data[clust$cluster == i,feature]))
  }
  pars <- c(mus,sigs)
  return(pars)
}

get_par_fine0 <- function(feature_data,include_corrs=T){
  
  mus <- quantile(feature_data,
                  probs = seq(0,1,1/(N_coarse+1)),
                  na.rm=TRUE)[2:(N_fine+1)]
  
  sigs <- log(rep(sd(feature_data,na.rm=TRUE),N_fine))
  
  if (include_corrs){
    mus <- 0.5*mus
    corrs <- rep(0.5,N_fine)
    
    pars <- c()
    for (i in 1:N_fine){
      pars <- c(pars,c(mus[i],corrs[i]))
    }
    pars <- rep(pars,N_coarse)
    pars <- c(pars,rep(sigs,times=N_coarse))
      
  } else {

    pars <- c(rep(mus,times=N_coarse),
              rep(sigs,times=N_coarse))    
  }
  
  return(pars)
}

coarse_Data <- Data[Data$level == 1,c("diveDuration","maxDepth")]
ParDiveDuration0 <- get_par_coarse0(coarse_Data,"diveDuration")
ParMaxDepth0 <- get_par_coarse0(coarse_Data,"maxDepth")

#ParDiveDuration0 <- rep(20,4)
#ParMaxDepth0  <- rep(20,4)

#Parhv0 <- get_par_fine0(Data$hv)
#Parjp0 <- get_par_fine0(Data$jp)
#Parrajp0 <- get_par_fine0(Data$rajp)
Paraw10 = rep(0,3*N_total)
#Paraw10[seq(1,2*N_total,2)] = sample(Data$aw1[!is.na(Data$aw1)], N_total)
Paraw10[seq(2,2*N_total,2)] = 1.0
Paraw10[(2*N_total+1):(3*N_total)] = log(sd(Data$d_aw1,na.rm=TRUE))

Paraw20 = rep(0,3*N_total)
#Paraw20[seq(1,2*N_total,2)] = sample(Data$aw2[!is.na(Data$aw2)], N_total)
Paraw20[seq(2,2*N_total,2)] = 1.0
Paraw20[(2*N_total+1):(3*N_total)] = log(sd(Data$d_aw2,na.rm=TRUE))

Paraw30 = rep(0,3*N_total)
#Paraw30[seq(1,2*N_total,2)] = sample(Data$aw3[!is.na(Data$aw3)], N_total)
Paraw30[seq(2,2*N_total,2)] = 1.0
Paraw30[(2*N_total+1):(3*N_total)] = log(sd(Data$d_aw3,na.rm=TRUE))

Parw0 = c(sample(Data$w[!is.na(Data$w)], N_total),
          log(rep(sd(Data$w,na.rm=TRUE),N_total)))

# get initial parameters
Par = list(diveDuration = rep(ParDiveDuration0,each=1),
           maxDepth = rep(ParMaxDepth0,each=1))

#Par <- getParDM(Data,
#                hierStates=hierStates,
#                hierDist=hierDist,
#                Par=Par,
#                DM=DM)

# include the other stuff after the fact
# hierDist$level2$AddChild("aw1", dist="norm")
# hierDist$level2$AddChild("aw2", dist="norm")
# hierDist$level2$AddChild("aw3", dist="norm")

# DM$aw1 <- list(mean = ~1 + aw1.tm1, sd = ~1)
# DM$aw2 <- list(mean = ~1 + aw2.tm1, sd = ~1)
# DM$aw3 <- list(mean = ~1 + aw3.tm1, sd = ~1)

# Par$aw1 = Paraw10
# Par$aw2 = Paraw20
# Par$aw3 = Paraw30

# check hierarchical model specifications
checkPar0(Data[1:1000,],
          hierStates = hierStates,
          hierDist = hierDist,
          hierBeta = hierBeta, 
          hierDelta = hierDelta,
          DM=DM,
          Par0=Par)

# # fit the HHMM
hhmm <- fitHMM(Data,
               hierStates = hierStates,
               hierDist = hierDist,
               hierBeta = hierBeta,
               hierDelta = hierDelta,
               DM=DM,
               Par0=Par,
               nlmPar = list('print.level'=2))
hhmm
saveRDS(hhmm,"hhmm_test.rds")
hhmm <- readRDS("../params/hhmm_02-25-2022.rds")

labels <- NULL
labels$level1 <- Data$label[Data$level == 1]
labels$level2 <- Data$label[Data$level != 1]
vstates <- viterbi(hhmm, hierarchical = TRUE)

# add the dive types and subdive states to the RawData
dive_types <- data.frame(vstate1 = vstates$level1,
                         label1 = labels$level1,
                         divenum = seq(1,length(vstates$level1)))

rawData <- rawData[,!(names(rawData) %in% c("vstate1","label1"))]
rawData <- left_join(rawData,dive_types,by="divenum")

subdive_states <- data.frame(vstate2 = vstates$level2,
                             label2 = labels$level2,
                             segnum = seq(1,length(vstates$level2)))
rawData$divesegnum = rawData$segnum + rawData$divenum
rawData <- rawData[,!(names(rawData) %in% c("vstate2","label2"))]
rawData <- left_join(rawData,subdive_states,by="segnum")

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
Data$vstate1[Data$vstate2 %in% c("S41","S42","S43")] <- "D4"

# add elevation
rawData$Elevation <- -rawData$p

# prepare the data for plotting
rawDataDownsampled <- rawData[seq(1,nrow(rawData),5),]
rawDataDownsampled$jp <- log(rawDataDownsampled$jp)
rawDataDownsampled$hv <- logit(rawDataDownsampled$hv)
rawDataDownsampled$roll <- rawDataDownsampled$roll*180/pi
rawDataDownsampled$head <- rawDataDownsampled$head*180/pi

# columns to plot

cols_to_plot = c("Elevation","w","Aw_1","Aw_2","Aw_3")

labs <- c(Elevation = "Depth (m)", 
          head = "Heading (degrees)",
          roll = "Roll (degrees)",
          pitch = "Pitch (degrees)",
          jp = "Jerk Peak (m/s^3)",
          htv = "heading total variation",
          w = "wiggliness",
          Aw_1 = "x-acc (m/s^2)",
          Aw_2 = "y-acc (m/s^2)",
          Aw_3 = "z-acc (m/s^2)")

rawDataDownLong <- rawDataDownsampled %>% 
  pivot_longer(cols = cols_to_plot, 
               names_to = "feature")

# plot the results
A <- rawDataDownLong[rawDataDownLong$divenum %in% 110:115,]
inds = rawDataDownLong[rawDataDownLong$divenum %in% 110:115,] %>% pull("divenum") %>% diff
inds = which(inds %in% c(1))
times = A[inds,] %>% pull("Time")

times <- Data[Data$label == "crunch","stime"]

# ggplot(rawDataDownLong[rawDataDownLong$divenum %in% 110:115,], 
#        aes(x=Time, y=value, group=divenum)) +
#   geom_line(aes(color=vstate1)) + 
#   geom_vline(xintercept = times) +
#   facet_wrap(~feature,ncol=1,scales='free_y')

ggplot(rawDataDownLong[(rawDataDownLong$divenum %in% 180:200) &
                       (rawDataDownLong$segnum != 0),], 
       aes(x=Time, y=value)) +
  geom_line(aes(color=vstate2, group=segnum)) + 
  geom_vline(xintercept = times) +
  facet_wrap(~feature,
             ncol=1,
             scales='free_y',
             strip.position = "left",
             labeller = as_labeller(labs)) +
  labs(color="Behavioural State",y=NULL) + ggtitle("Shallow Behaviour of Whale I107 on August 25, 2020") +
  theme(strip.background = element_blank(),
        strip.placement = "outside")
  
ggplot(rawDataDownLong[(rawDataDownLong$divenum %in% 180:200) &
                         (rawDataDownLong$segnum != 0),], 
       aes(x=Time, y=value)) +
  geom_line(aes(color=vstate2, group=segnum)) + 
  geom_vline(xintercept = times) +
  facet_wrap(~feature,
             ncol=1,
             scales='free_y',
             strip.position = "left",
             labeller = as_labeller(labs)) +
  labs(color="Subdive State",y=NULL) + ggtitle("Within Dive Behaviour of Whale I107 on August 25, 2020") +
  theme(strip.background = element_blank(),
        strip.placement = "outside")

ggplot(rawDataDownLong[(rawDataDownLong$divenum %in% 180:200) &
                         (rawDataDownLong$segnum != 0),], 
       aes(x=Time, y=value)) +
  geom_line(aes(color=vstate1, group=divenum)) + 
  geom_vline(xintercept = times) +
  facet_wrap(~feature,
             ncol=1,
             scales='free_y',
             strip.position = "left",
             labeller = as_labeller(labs)) +
  labs(color="Dive Type",y=NULL) + ggtitle("Dive Types of Whale I107 on August 25, 2020") +
  theme(strip.background = element_blank(),
        strip.placement = "outside")

# plot each feature individually
plot_coarse_feature <- function(feature){
  shape1 = hhmm$mle[[feature]]['mean','S11']^2 / hhmm$mle[[feature]]['sd','S11']^2
  rate1 = hhmm$mle[[feature]]['mean','S11'] / hhmm$mle[[feature]]['sd','S11']^2
  shape2 = hhmm$mle[[feature]]['mean','S21']^2 / hhmm$mle[[feature]]['sd','S21']^2
  rate2 = hhmm$mle[[feature]]['mean','S21'] / hhmm$mle[[feature]]['sd','S21']^2
  shape3 = hhmm$mle[[feature]]['mean','S31']^2 / hhmm$mle[[feature]]['sd','S31']^2
  rate3 = hhmm$mle[[feature]]['mean','S31'] / hhmm$mle[[feature]]['sd','S31']^2
  shape4 = hhmm$mle[[feature]]['mean','S41']^2 / hhmm$mle[[feature]]['sd','S41']^2
  rate4 = hhmm$mle[[feature]]['mean','S41'] / hhmm$mle[[feature]]['sd','S41']^2
  
  plot <- ggplot(Data,aes_string(x=feature)) +
    geom_histogram(aes(fill=vstate1,y=stat(density)),alpha=0.5,position = "identity") + 
    stat_function(fun = function(x) {dgamma(x,shape=shape1,rate=rate1)},
                  aes(color = 'D1'),) +
    stat_function(fun = function(x) {dgamma(x,shape=shape2,rate=rate2)},
                  aes(color = 'D2')) +
    stat_function(fun = function(x) {dgamma(x,shape=shape3,rate=rate3)},
                  aes(color = 'D3')) +
    stat_function(fun = function(x) {dgamma(x,shape=shape4,rate=rate4)},
                  aes(color = 'D4'))
  
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

print(plot_fine_feature_norm('w'))
print(plot_fine_feature_norm('aw1'))
print(plot_fine_feature_norm('aw2'))
print(plot_fine_feature_norm('aw3'))

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

