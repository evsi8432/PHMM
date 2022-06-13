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
set.seed(0)

### Parse Command line arguments ###
args = commandArgs(trailingOnly=TRUE)
#args <- 5

### Load Data ###
setwd("~/Documents/Research/PHMM/src")

### start by saying how many states we will have ###
N_deep <- 2
N_shallow <- 2

N_coarse <- N_deep + N_shallow

# load data
Data <- data.frame(fread('../dat/diary/Final_Data.csv'))
rawData <- data.frame(fread('../dat/diary/Final_rawData.csv'))
ethogram <- data.frame(fread('../dat/diary/Final_ethogram.csv'))

### Label Deep vs Shallow Dives ###
shallow_thresh <- log(5)
deep_thresh <- log(20)

Data$broadDiveType <- NA #3  # unknown
Data$broadDiveType[Data$maxDepth > deep_thresh] <- TRUE  # deep
Data$broadDiveType[Data$maxDepth < shallow_thresh] <- FALSE  # shallow

### Set initial parameters ###
features <- c("maxDepth","diveDuration",
              "avg_bot_htv","avg_bot_abs_roll","max_bot_jp")

shallowData <- Data[Data$maxDepth < log(shallow_thresh) & complete.cases(Data[,features]),]
shallowFit <- Mclust(shallowData[,features],
                     G=N_shallow,modelNames= "VVI",na.rm=TRUE)

deepData <- Data[Data$maxDepth > log(deep_thresh) & complete.cases(Data[,features]),]
deepFit <- Mclust(deepData[,features], 
                   G=N_deep,modelNames = "VVI", na.rm=TRUE)

get_norm_par <- function(shallowData,shallowFit,
                         deepData,deepFit,
                         feature){
  
  Par0 <- rep(0,2*N_coarse)
  
  for (i in 1:N_shallow){
    Par0[i] <- mean(shallowData[shallowFit$classification == i,feature],
                    na.rm=TRUE)
    Par0[N_coarse + i] <- sd(shallowData[shallowFit$classification == i,feature],
                           na.rm=TRUE)
  }
  for (i in 1:N_deep){
    Par0[N_shallow + i] <- mean(deepData[deepFit$classification == i,feature],
                                na.rm=TRUE)
    Par0[N_coarse + N_shallow + i] <- sd(deepData[deepFit$classification == i,feature],
                                         na.rm=TRUE)
  }
  
  Par0 <- replace_na(Par0,sd(deepData$feature,na.rm=TRUE))
  
  return(Par0)
}

get_beta_par <- function(shallowData,shallowFit,
                         deepData,deepFit,
                         feature){
  
  Par0 <- rep(0,2*N_coarse)
  
  for (i in 1:N_shallow){
    Par0[0 * N_coarse + i] <- 1
    Par0[1 * N_coarse + i] <- 1
    #Par0[2 * N_coarse + i] <- max(min(sum(shallowData[shallowFit$classification == i,][[feature]] == 0,
    #                              na.rm=TRUE),0.99),0.01)
    #Par0[3 * N_coarse + i] <- max(min(sum(shallowData[shallowFit$classification == i,][[feature]] == 1,
    #                              na.rm=TRUE),0.99),0.01)
  }
  for (i in 1:N_deep){
    Par0[0 * N_coarse + N_shallow + i] <- 1
    Par0[1 * N_coarse + N_shallow + i] <- 1
    #Par0[2 * N_coarse + N_shallow + i] <- max(min(sum(deepData[deepFit$classification == i,][[feature]] == 0,
    #                                          na.rm=TRUE),0.99),0.01) 
    #Par0[3 * N_coarse + N_shallow + i] <- max(min(sum(deepData[deepFit$classification == i,][[feature]] == 1,
    #                                          na.rm=TRUE),0.99),0.01)
  }
  
  Par0 <- replace_na(Par0,0.5)
  
  return(Par0)
}

ParMaxDepth0 <- get_norm_par(shallowData,shallowFit,
                             deepData,deepFit,
                             "maxDepth")
ParDiveDuration0 <- get_norm_par(shallowData,shallowFit,
                                 deepData,deepFit,
                                 "diveDuration")
ParAscentSpeed0 <- get_norm_par(shallowData,shallowFit,
                                deepData,deepFit,
                                "ascentSpeed")
ParBottomProp0 <- get_norm_par(shallowData,shallowFit,
                               deepData,deepFit,
                               "bottomProp")
ParHTV0 <- get_norm_par(shallowData,shallowFit,
                        deepData,deepFit,
                        "avg_bot_htv")
ParRoll0 <- get_norm_par(shallowData,shallowFit,
                         deepData,deepFit,
                         "avg_bot_abs_roll")
ParJerkPeak0 <- get_norm_par(shallowData,shallowFit,
                             deepData,deepFit,
                             "max_bot_jp")
ParBroadDiveType0 <- c(
  rep(1e-10, N_shallow), rep(1.0-1e-10,N_deep)
)

get_bern_par <- function(shallowData,shallowFit,
                         deepData,deepFit,
                         feature){
  
  Par0 <- rep(0,N_coarse)
  
  for (i in 1:N_shallow){
    Par0[i] <- mean(shallowData[shallowFit$classification==i,feature],
                    na.rm=TRUE)
  }
  for (i in 1:N_deep){
    Par0[N_shallow + i] <- mean(deepData[deepFit$classification==i,feature],
                                na.rm=TRUE)
  }
  
  Par0 <- pmax(Par0,0.01)
  Par0 <- pmin(Par0,0.99)
  Par0 <- replace_na(Par0,0.5)
  
  return(Par0)
}

ParForaging0 <- get_bern_par(shallowData,shallowFit,
                             deepData,deepFit,
                             "foraging")

ParScales0 <- get_bern_par(shallowData,shallowFit,
                             deepData,deepFit,
                             "Scales")

ParClicks0 <- get_bern_par(shallowData,shallowFit,
                           deepData,deepFit,
                           "clicks")

ParClickTrain0 <- get_bern_par(shallowData,shallowFit,
                           deepData,deepFit,
                           "click.train")

ParBuzz0 <- get_bern_par(shallowData,shallowFit,
                         deepData,deepFit,
                         "buzz")

ParWhales0 <- get_bern_par(shallowData,shallowFit,
                           deepData,deepFit,
                           "whales")

ParTouch0 <- get_bern_par(shallowData,shallowFit,
                           deepData,deepFit,
                           "touch")

ParKelp0 <- get_bern_par(shallowData,shallowFit,
                          deepData,deepFit,
                          "kelp")

Par0 = list(maxDepth = ParMaxDepth0,
            diveDuration = ParDiveDuration0,
            broadDiveType = ParBroadDiveType0,
            #bottomProp = ParBottomProp0,
            #ascentSpeed = ParAscentSpeed0,
            avg_bot_htv=ParHTV0,
            avg_bot_abs_roll=ParRoll0,
            max_bot_jp=ParJerkPeak0,
            foraging = ParForaging0,
            scales = ParScales0,
            #clicks = ParClicks0,
            click.train = ParClickTrain0,
            #buzz = ParBuzz0,
            whales = ParWhales0,
            touch = ParTouch0,
            kelp = ParKelp0)

### coarse scale Distributions ###

dist <- list(diveDuration="norm",
             maxDepth="norm",
             broadDiveType="bern",
             #bottomProp="norm",
             #ascentSpeed="norm",
             avg_bot_htv="norm",
             avg_bot_abs_roll="norm",
             max_bot_jp="norm",
             foraging="bern",
             scales="bern",
             #clicks="bern",
             click.train="bern",
             #buzz="bern",
             whales="bern",
             touch="bern",
             kelp="bern")


### Fix Known Parameters ###

fixPar = list()

fixParBroadDiveType0 <- c(
  #rep(0.5, N_shallow), rep(1e-10,N_deep), # prob shallow, each state (intercept * 4)
  rep(1e-10, N_shallow), rep(1.0-1e-10,N_deep)  # prob deep, each state (intercept * 4)
)

fixPar$broadDiveType <- fixParBroadDiveType0

# prep data idk why
Data <- prepData(Data,
                 coordNames=NULL)

### Check Parameters ###
checkPar0(Data,
          nbStates=N_coarse,
          dist=dist,
          fixPar=fixPar,
          Par0=Par0)


### Fit Model ###

maxiters <- as.numeric(args[1])
#maxiters <- 130

ptm <- proc.time()

if (maxiters == 0 | is.na(maxiters)){
  hhmm <- fitHMM(Data,
                 nbStates = N_coarse,
                 dist = dist,
                 fixPar=fixPar,
                 Par0=Par0,
                 #retryFits = 5,
                 nlmPar = list('print.level'=2))
} else {
  print(paste("max iters:",maxiters))
  hhmm <- fitHMM(Data,
                 nbStates = N_coarse,
                 dist = dist,
                 fixPar=fixPar,
                 Par0=Par0,
                 #retryFits = 5,
                 nlmPar = list('print.level'=2,
                               'iterlim'=maxiters))
}

print(proc.time()-ptm)

hhmm

saveRDS(hhmm,"../params/hhmm_June_7.rds")





# Do some small EDA

setwd("~/Documents/Research/PHMM/src")
hhmm <- readRDS("../params/hhmm_June_7.rds")

vstates <- factor(viterbi(hhmm))

# add dive types to Data
Data$vstate1 <- vstates

# add the dive types to the RawData
dive_types <- data.frame(vstate1 = vstates,
                         noBottom = Data$noBottom,
                         divenum = seq(1,length(vstates)))

rawData <- rawData[,!(names(rawData) %in% c("vstate1","label1"))]
rawData <- left_join(rawData,dive_types,by="divenum")

rawData$Elevation <- -rawData$p

# prepare the data for plotting
rawDataDownsampled <- rawData[seq(1,nrow(rawData),5),]
rawDataDownsampled$jp <- log(rawDataDownsampled$jp)
rawDataDownsampled$htv <- log(rawDataDownsampled$htv)
rawDataDownsampled$roll <- rawDataDownsampled$roll*180/pi
rawDataDownsampled$head <- rawDataDownsampled$head*180/pi

# columns to plot

cols_to_plot = c("Elevation","head","pitch","roll")

labs <- c(Elevation = "Depth (m)", 
          head = "Heading (degrees)",
          roll = "Roll (degrees)",
          pitch = "Pitch (degrees)",
          jp = "Jerk Peak (m/s^3)",
          htv = "Heading Total Variation",
          w = "wiggliness",
          Aw_1 = "x-acc (m/s^2)",
          Aw_2 = "y-acc (m/s^2)",
          Aw_3 = "z-acc (m/s^2)",
          maxDepth = "log(Maximum Depth (m))",
          diveDuration = "log(Dive Duration (s))",
          max_bot_jp = "log(Bottom Peak Jerk (m/s^3))",
          avg_bot_htv = "log(Average Heading Variation Rate at Bottom (rad/s))",
          avg_bot_abs_roll = "log(Average Absolute Roll at Bottom (rad))")

rawDataDownLong <- rawDataDownsampled %>% 
  pivot_longer(cols = cols_to_plot, 
               names_to = "feature")

# plot the results
ac_ethogram <- read.csv("../dat/ethograms/echolocation_ethogram.csv")

ac_ethogram$date0 <- NA
ac_ethogram[ac_ethogram$Animal == "I107","date0"] <- "2020-08-25"
ac_ethogram[ac_ethogram$Animal == "I145","date0"] <- "2020-08-30"
ac_ethogram[ac_ethogram$Animal == "D26" ,"date0"] <- "2020-08-31"

ac_ethogram$stime <- as.POSIXct(paste(ac_ethogram$date0,ac_ethogram$time),format="%Y-%m-%d %H:%M:%S")
ac_ethogram$etime <- ac_ethogram$stime + ac_ethogram$duration.sec

ac_times <- ac_ethogram[ac_ethogram$broad.type %in% c('clicks','click train','buzz'),
                        c('stime','etime','broad.type')]

names(ac_times) <- c('stime','etime','label')

times <- ethogram[ethogram$Behavior %in% c('Kill','Confirmed Fish Capture','Scales','whale','touch'),
                  c('stime','etime','Behavior')]
names(times) <- c('stime','etime','label')

times <- rbind(ac_times,times)

#sclicks <- ac_ethogram$stime[!is.na(ac_ethogram$broad.type) & 
#                             (ac_ethogram$broad.type == "clicks")]

#eclicks <- ac_ethogram$etime[!is.na(ac_ethogram$broad.type) & 
#                               (ac_ethogram$broad.type == "clicks")]
#nclicks <- length(stimes)

#times <- data.frame(xmin=stimes,
#                    xmax=etimes, 
#                    ymin = rep(-Inf,ntimes),
#                    ymax = rep(Inf,ntimes))

times <- transform(times, 
                   id = 1:nrow(times),
                   ymin = -Inf,
                   ymax = Inf)

# ggplot(rawDataDownLong[rawDataDownLong$divenum %in% 110:115,], 
#        aes(x=Time, y=value, group=divenum)) +
#   geom_line(aes(color=vstate1)) + 
#   geom_vline(xintercept = times) +
#   facet_wrap(~feature,ncol=1,scales='free_y')

print(head(Data[Data$ID == "R48",]))
print(tail(Data[Data$ID == "R48",]))

sdive = 7927 #3800
edive = 8235 #3950

ggplot(rawDataDownLong[(rawDataDownLong$divenum %in% sdive:edive) &
                       (rawDataDownLong$segnum != 0),], 
       aes(x=Time, y=value)) +
  geom_rect(data=times, inherit.aes=FALSE,
            aes(xmin=stime,xmax=etime,ymin=ymin,ymax=ymax,
                group=id), alpha=0.5) +
  #geom_vline(xintercept = times$stime) +
  geom_line(aes(color=factor(vstate1), group=divenum)) + 
  facet_wrap(~feature,
             ncol=1,
             scales='free_y',
             strip.position = "left",
             labeller = as_labeller(labs)) +
  labs(color="Dive Type",y=NULL) + ggtitle("Behaviour of Whale R48 on August 25, 2020") +
  theme(strip.background = element_blank(),
        strip.placement = "outside") + 
  xlim(Data$stime[sdive],Data$etime[edive])

for(feature in c("maxDepth","diveDuration","avg_bot_htv","avg_bot_abs_roll","max_bot_jp")){
  mu1 = hhmm$mle[[feature]]['mean','state 1']
  sig1 = hhmm$mle[[feature]]['sd','state 1']
  mu2 = hhmm$mle[[feature]]['mean','state 2']
  sig2 = hhmm$mle[[feature]]['sd','state 2']
  mu3 = hhmm$mle[[feature]]['mean','state 3']
  sig3 = hhmm$mle[[feature]]['sd','state 3']
  mu4 = hhmm$mle[[feature]]['mean','state 4']
  sig4 = hhmm$mle[[feature]]['sd','state 4']
  mu5 = hhmm$mle[[feature]]['mean','state 5']
  sig5 = hhmm$mle[[feature]]['sd','state 5']
  #mu6 = hhmm$mle[[feature]]['mean','state 6']
  #sig6 = hhmm$mle[[feature]]['sd','state 6']
  
  plot <- ggplot(Data,aes_string(x=feature)) +
    geom_histogram(aes(fill=vstate1,y=stat(density)),alpha=0.5,position = "identity") + 
    stat_function(fun = function(x) {dnorm(x,mean=mu1,sd=sig1)},
                  aes(color = 'state 1'),) +
    stat_function(fun = function(x) {dnorm(x,mean=mu2,sd=sig2)},
                  aes(color = 'state 2')) +
    stat_function(fun = function(x) {dnorm(x,mean=mu3,sd=sig3)},
                  aes(color = 'state 3')) +
    stat_function(fun = function(x) {dnorm(x,mean=mu4,sd=sig4)},
                  aes(color = 'state 4')) +
    stat_function(fun = function(x) {dnorm(x,mean=mu5,sd=sig5)},
                  aes(color = 'state 5')) + 
    labs(color = "Dive Type", x = labs[feature]) + 
    ggtitle(paste("Histogram for Killer Whale Dives from 2020 Field Season"))
    #+
    #stat_function(fun = function(x) {dnorm(x,mean=mu6,sd=sig6)},
    #              aes(color = 'state 6'))
  print(plot)
}

for(feature in c("bottomProp")){
  a1 = hhmm$mle[[feature]]['shape1','state 1']
  b1 = hhmm$mle[[feature]]['shape2','state 1']
  a2 = hhmm$mle[[feature]]['shape1','state 2']
  b2 = hhmm$mle[[feature]]['shape2','state 2']
  a3 = hhmm$mle[[feature]]['shape1','state 3']
  b3 = hhmm$mle[[feature]]['shape2','state 3']
  a4 = hhmm$mle[[feature]]['shape1','state 4']
  b4 = hhmm$mle[[feature]]['shape2','state 4']
  a5 = hhmm$mle[[feature]]['shape1','state 5']
  b5 = hhmm$mle[[feature]]['shape2','state 5']
  
  plot <- ggplot(Data,aes_string(x=feature)) +
    geom_histogram(aes(fill=vstate1,y=stat(density)),alpha=0.5,position = "identity") + 
    stat_function(fun = function(x) {dbeta(x,shape1=a1,shape2=b1)},
                  aes(color = 'state 1'),) +
    stat_function(fun = function(x) {dbeta(x,shape1=a2,shape2=b2)},
                  aes(color = 'state 2')) +
    stat_function(fun = function(x) {dbeta(x,shape1=a3,shape2=b3)},
                  aes(color = 'state 3')) +
    stat_function(fun = function(x) {dbeta(x,shape1=a4,shape2=b4)},
                  aes(color = 'state 4')) +
    stat_function(fun = function(x) {dbeta(x,shape1=a5,shape2=b5)},
                  aes(color = 'state 5'))
  print(plot)
}

ggplot(data = Data[Data$vstate1 %in% c(1,2),], aes_string(x="max_bot_jp",
                                                          y="avg_bot_htv", 
                                                          color="vstate1")) +
  geom_point() + 
  labs(color="Dive Type", 
       x = "log(Bottom Peak Jerk (m/s^3))", 
       y = "log(Average Heading Variation Rate at Bottom (rad/s))") + 
  ggtitle("Kinematics of Shallow Killer Whale Dives from 2020 Field Season")

ggplot(data = Data[Data$vstate1 %in% c(1,2),], aes_string(x="avg_bot_abs_roll",
                                                            y="avg_bot_htv", 
                                                            color="vstate1")) +
  geom_point() + 
  labs(color="Dive Type", 
       x = "log(Average Absolute Roll at Bottom (rad))", 
       y = "log(Average Heading Variation Rate at Bottom (rad/s))") + 
  ggtitle("Kinematics of Shallow Killer Whale Dives from 2020 Field Season")

ggplot(data = Data[Data$vstate1 %in% c(3,4,5),], aes_string(x="max_bot_jp",
                                                          y="avg_bot_htv", 
                                                          color="vstate1")) +
  geom_point() + 
  labs(color="Dive Type", 
       x = "log(Bottom Peak Jerk (m/s^3))", 
       y = "log(Average Heading Variation Rate at Bottom (rad/s))") + 
  ggtitle("Kinematics of Deep Killer Whale Dives from 2020 Field Season")

ggplot(data = Data[Data$vstate1 %in% c(3,4,5),], aes_string(x="avg_bot_abs_roll",
                                                            y="avg_bot_htv", 
                                                            color="vstate1")) +
  geom_point() + 
  labs(color="Dive Type", 
       x = "log(Average Absolute Roll at Bottom (rad))", 
       y = "log(Average Heading Variation Rate at Bottom (rad/s))") + 
  ggtitle("Kinematics of Deep Killer Whale Dives from 2020 Field Season")
  

ggplot(data = Data, aes_string(x="maxDepth",
                               y="diveDuration", 
                               color="vstate1")) +
  geom_point() + 
  labs(color="Dive Type", x = "log(Maximum Depth (m))", y = "log(Dive Duration (s))") + 
  ggtitle("Shapes of all Killer Whale Dives from 2020 Field Season")

probs <- data.frame(rbind(hhmm$mle$kelp,
                          hhmm$mle$touch,
                          hhmm$mle$whales,
                          hhmm$mle$click.train,
                          hhmm$mle$scales,
                          hhmm$mle$foraging))

probs$event <- c("kelp interaction","whale interaction","whales in view",
                 "click train","visible scales","confirmed foraging")
colnames(probs) <- c(1:5,"event")

probs <- probs %>% pivot_longer(cols = 1:5,
                                names_to = "Dive Type")

ggplot(data=probs[probs$event %in% c("click train","visible scales","confirmed foraging"),], 
       aes(x=event,y=value,fill=`Dive Type`)) + 
  geom_bar(stat="identity",position="dodge") + 
  labs(y = "Probability of Observation") + 
  ggtitle("Probability of Foraging Events by Dive Type")

ggplot(data=probs[probs$event %in% c("kelp interaction","whale interaction","whales in view"),], 
       aes(x=event,y=value,fill=`Dive Type`)) + 
  geom_bar(stat="identity",position="dodge") + 
  labs(y = "Probability of Observation") + 
  ggtitle("Probability of Social Events by Dive Type")

A <- rep(1/N_coarse,N_coarse)

for(i in 1:100){
  A <- A %*% hhmm$mle$gamma
}
A <- A * exp(hhmm$mle$diveDuration["mean",])
A <- A / sum(A)
A <- t(A)
colnames(A) <- "prop"
A <- data.frame(A)
A$dive.type <- factor(1:5)

ggplot(data=data.frame(A), 
       aes(x=dive.type, y=prop, fill=dive.type)) + 
  geom_bar(stat="identity") + 
  labs(y = "Propotion of Time") + 
  ggtitle("Proportion of Time in Each Dive Type")

