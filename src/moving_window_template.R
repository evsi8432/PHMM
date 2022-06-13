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
set.seed(0)

setwd("~/Documents/Research/PHMM/src")

# define summary stats function
SumStats <- function(x){
  
  # get rid of NA rows
  if(sum(complete.cases(x[,c("head","pitch","roll")])) > (nrow(x)/2)){
    x <- x[complete.cases(x[,c("head","pitch","roll")]),]
  }
  
  # smooth out heading, pitch, roll, Aw_1-3
  smooth_head <- ksmooth(x$Time,abs(x$head),kernel="normal",bandwidth=0.25)$y
  smooth_pitch <- ksmooth(x$Time,x$pitch,kernel="normal",bandwidth=0.25)$y
  smooth_roll <- ksmooth(x$Time,abs(x$roll),kernel="normal",bandwidth=0.25)$y
  
  smooth_Aw_1 <- ksmooth(x$Time,x$dyn.aw1,kernel="normal",bandwidth=0.25)$y
  smooth_Aw_2 <- ksmooth(x$Time,x$dyn.aw2,kernel="normal",bandwidth=0.25)$y
  smooth_Aw_3 <- ksmooth(x$Time,x$dyn.aw3,kernel="normal",bandwidth=0.25)$y
  
  # get Jerk Peak 
  jerk <- 50*sqrt(diff(x$Aw_1)^2 + 
                  diff(x$Aw_2)^2 + 
                  diff(x$Aw_3)^2)
  
  jp <- max(jerk,na.rm=TRUE)
  
  # find roll at Jerk Peak
  if (jp > -Inf){
    jpind <- which(jerk == jp) 
    rajp <- x$roll[jpind]
  } 
  else {
    rajp <- NaN
  }
  
  # find total variation of heading, pitch, and roll
  htv <- sum(abs((diff(smooth_head)+pi) %% (2*pi) - pi)) / 2
  ptv <- sum(abs((diff(smooth_pitch)+pi) %% (2*pi) - pi)) / 2
  rtv <- sum(abs((diff(smooth_roll)+pi) %% (2*pi) - pi)) / 2
  
  # get average depth
  ad <- mean(x$p)
  
  # get average acceleration
  aw1 <- mean(smooth_Aw_1)
  aw2 <- mean(smooth_Aw_2)
  aw3 <- mean(smooth_Aw_3)
  
  # get average heading, pitch, and roll
  head <- mean(x$head)
  pitch <- mean(x$pitch)
  roll <- mean(x$roll)
  
  # get wiggliness
  w <- Mod(fft(x$Aw_1))^2 + Mod(fft(x$Aw_2))^2 + Mod(fft(x$Aw_3))^2
  w <- sum(w[2:10])
  
  # get time
  stime <- min(x$Time)
  etime <- max(x$Time)
  
  # get VeDBA after getting dynamic acc
  VeDBA <- mean(x$VeDBA)
  
  # return everything
  return(c(jp,rajp,ad,htv,ptv,rtv,aw1,aw2,aw3,w,head,pitch,roll,VeDBA,stime,etime))
}

Prep_Data <- function(whale){
  
  # load in data for whale of interest
  print(whale)
  print("reading data")
  rawData <- read.csv(paste0("../dat/diary/",whale,"/kinematic_data_calibrated_",
                             whale,"_50Hz.csv"))
  
  # convert to a date time
  rawData$Time <- as.POSIXct(rawData$Time,format="%d-%b-%Y %H:%M:%OS")
  date0 = format(rawData$Time[1], "%Y-%m-%d")  # prone to error for overnight data
  
  # get dynamic Acceleration and VeDBA(subtract gravity)
  g = 1 #g's
  
  gw1 <- g*sin(rawData$pitch)
  gw2 <- -g*cos(rawData$pitch)*sin(rawData$roll)
  gw3 <- -g*cos(rawData$pitch)*cos(rawData$roll)
  
  rawData$dyn.aw1 <- rawData$Aw_1 - gw1
  rawData$dyn.aw2 <- rawData$Aw_2 - gw2
  rawData$dyn.aw3 <- rawData$Aw_3 - gw3
  
  rawData$VeDBA <- sqrt(rawData$dyn.aw1^2+rawData$dyn.aw2^2+rawData$dyn.aw3^2)
  
  # label dives and get phases using diveMove
  rawDataDown <- rawData[seq(1,nrow(rawData),25),]
  TDRData <- createTDR(rawDataDown$Time, rawDataDown$p, file="temp.csv")
  
  
  # use 3 m threshold
  TDRData <- calibrateDepth(TDRData, zoc.method = "filter",
                            dive.thr=1, k=c(5,30*60) ,probs=c(0.5,0.02))
  diveSums <- diveStats(TDRData)
  
  divestarts <- diveSums$begdesc
  bottstarts <- diveSums$enddesc
  bottends <- diveSums$begasc
  diveends <- divestarts + diveSums$divetim
  
  h <- 100 # number of observations in a moving window
  
  # initialize things in raw data
  rawData$divenum <- 0
  rawData$segnum <- 0
  rawData$bottom <- 0
  segnum <- 0
  
  # define sum stats
  rawData$jp <- NA
  rawData$rajp <- NA
  rawData$htv <- NA
  rawData$ptv <- NA
  rawData$rtv <- NA
  rawData$w <- NA
  
  # define a full data matrix
  Data <- NULL
  
  print("Processing Dives")
  for (i in 1:(length(divestarts))){
    
    print(i/length(divestarts))
    
    # get dive start / end
    startind <- min(which(rawData$Time >= divestarts[i]))
    endind <- max(which(rawData$Time <= diveends[i]))
      
    # get dataframe with dives
    diveData <- rawData[startind:endind,]
    
    # get bottom phase
    startBotind <- min(which((rawData$Time >= divestarts[i]) & 
                               rawData$p >= 0.8*max(diveData$p)))
    endBotind <- max(which((rawData$Time <= diveends[i]) & 
                             rawData$p >= 0.8*max(diveData$p)))
    
    # label the rawData
    rawData$divenum[startind:endind] <- i
    rawData$bottom[startBotind:endBotind] <- 1
    
    # get coarse-scale data
    dd <- as.double(difftime(max(diveData$Time),min(diveData$Time)),units='secs')
    md <- max(diveData$p)
    bd <- diveSums$bottdist[i]/ diveSums$botttim[i] # average vertical dist at bottom (m/s)
    pb <- diveSums$botttim[i]/ diveSums$divetim[i] # proportion of time at bottom
    as <- diveSums$ascdist[i]/ diveSums$asctim[i] # speed of ascent
    
    coarseScale <- data.frame(ID = whale,
                              level="1",
                              diveDuration = dd,
                              maxDepth = md,
                              bottomDist = bd,
                              bottomProp = pb,
                              ascentSpeed = as,
                              jp = NA,
                              rajp = NA,
                              htv = NA,
                              ptv = NA,
                              rtv = NA,
                              ad = NA,
                              aw1 = NA,
                              aw2 = NA,
                              aw3 = NA,
                              aw1.tm1 = NA,
                              aw2.tm1 = NA,
                              aw3.tm1 = NA,
                              w = NA,
                              head = NA,
                              pitch = NA,
                              roll = NA,
                              head.tm1 = NA,
                              pitch.tm1 = NA,
                              roll.tm1 = NA,
                              VeDBA = NA,
                              bottom = NA,
                              stime = as.numeric(rawData$Time[startind]),
                              etime = as.numeric(rawData$Time[endind]))
    
    print(h)
    print(nrow(diveData))
    print("")
    
    # get moving window values
    if(nrow(diveData) > 2*h){
      RawFineScale <- runner(
        x = diveData,
        k = h,
        at = seq(h+1,nrow(diveData),h),
        f = SumStats
      )
    } else {
      RawFineScale <- matrix(,nrow = 16,ncol = 0)
    }
    
    nsegs <- dim(RawFineScale)[2]
    print(nsegs)
    if (typeof(RawFineScale) == "list"){
      RawFineScale <- t(do.call(rbind, RawFineScale))
    }
    
    # label the segment numbers in the rawData
    for (j in seq_along(RawFineScale[15,])){
      segnum = segnum + 1
      segstartind <- which(rawData$Time == RawFineScale[15,j])
      segendind <- which(rawData$Time == RawFineScale[16,j])
      rawData$segnum[segstartind:segendind] <- segnum
      rawData$jp[segstartind:segendind] <- RawFineScale[1,j]
      rawData$rajp[segstartind:segendind] <- RawFineScale[2,j]
      rawData$htv[segstartind:segendind] <- RawFineScale[4,j]
      rawData$ptv[segstartind:segendind] <- RawFineScale[5,j]
      rawData$rtv[segstartind:segendind] <- RawFineScale[6,j]
      rawData$w[segstartind:segendind] <- RawFineScale[10,j]
    }
    
    if(nrow(diveData) > 2*h){
      fineScale <- data.frame(ID = whale,
                              level="2",
                              diveDuration = NA,
                              maxDepth = NA,
                              bottomDist = NA,
                              bottomProp = NA,
                              ascentSpeed = NA,
                              jp = RawFineScale[1,-1],
                              rajp = RawFineScale[2,-1],
                              ad = RawFineScale[3,-1],
                              htv = RawFineScale[4,-1],
                              ptv = RawFineScale[5,-1],
                              rtv = RawFineScale[6,-1], 
                              aw1 = RawFineScale[7,-1],
                              aw2 = RawFineScale[8,-1],
                              aw3 = RawFineScale[9,-1],
                              aw1.tm1 = RawFineScale[7,-nsegs],
                              aw2.tm1 = RawFineScale[8,-nsegs],
                              aw3.tm1 = RawFineScale[9,-nsegs],
                              w = RawFineScale[10,-1],
                              head = RawFineScale[11,-1],
                              pitch = RawFineScale[12,-1],
                              roll = RawFineScale[13,-1],
                              head.tm1 = RawFineScale[11,-nsegs],
                              pitch.tm1 = RawFineScale[12,-nsegs],
                              roll.tm1 = RawFineScale[13,-nsegs],
                              VeDBA = RawFineScale[14,-1],
                              bottom = (RawFineScale[16,-1] > bottstarts[i]) &
                                       (RawFineScale[15,-1] < bottends[i]),
                              stime = RawFineScale[15,-1],
                              etime = RawFineScale[16,-1])
    } else {
      fineScale <- data.frame(ID = whale,
                              level="2",
                              diveDuration = NA,
                              maxDepth = NA,
                              bottomDist = NA,
                              bottomProp = NA,
                              ascentSpeed = NA,
                              jp = NA,
                              rajp = NA,
                              ad = NA,
                              htv = NA,
                              ptv = NA,
                              rtv = NA, 
                              aw1 = NA,
                              aw2 = NA,
                              aw3 = NA,
                              aw1.tm1 = NA,
                              aw2.tm1 = NA,
                              aw3.tm1 = NA,
                              w = NA,
                              head = NA,
                              pitch = NA,
                              roll = NA,
                              head.tm1 = NA,
                              pitch.tm1 = NA,
                              roll.tm1 = NA,
                              VeDBA = NA,
                              bottom = NA,
                              stime = NA,
                              etime = NA)
    }
      
      fineScale <- rbind(data.frame(ID = whale,
                                    level="2i",
                                    diveDuration = NA,
                                    maxDepth = NA,
                                    bottomDist = NA,
                                    bottomProp = NA,
                                    ascentSpeed = NA,
                                    jp = NA,
                                    rajp = NA,
                                    ad = NA,
                                    htv = NA,
                                    ptv = NA,
                                    rtv = NA,
                                    aw1 = NA,
                                    aw2 = NA,
                                    aw3 = NA,
                                    aw1.tm1 = NA,
                                    aw2.tm1 = NA,
                                    aw3.tm1 = NA,
                                    w = NA,
                                    head = NA,
                                    pitch = NA,
                                    roll = NA,
                                    head.tm1 = NA,
                                    pitch.tm1 = NA,
                                    roll.tm1 = NA,
                                    VeDBA = NA,
                                    bottom = NA,
                                    stime = NA,
                                    etime = NA),
                         fineScale)
    print("")
    Data <- rbind(Data,coarseScale,fineScale)
  }
  
  summary(Data)
  
  # add labels
  ethogram <- read.csv("../dat/ethograms/echolocation_ethogram.csv")
  ethogram <- ethogram[ethogram$Animal == whale,]
  
  if (nrow(ethogram) > 0) {
    ethogram$stime <- as.POSIXct(paste(date0, ethogram$time))
    ethogram$stime <- as.numeric(ethogram$stime)
    ethogram$etime <- ethogram$stime + ethogram$duration.sec
  }
  
  Data$label <- 0
  
  for (i in 1:length(ethogram$broad.type)){
    inds <- which((Data$stime <= ethogram$etime[i]) & 
                    (Data$etime >= ethogram$stime[i]) & 
                    (Data$level != 1))
    print(ethogram$broad.type[i])
    print(ethogram$stime[i])
    print(ethogram$etime[i])
    print(inds)
    Data$label[inds] <- ethogram$broad.type[i]
  }
  
  # preprare data somehow idk why this is needed
  Data <- prepData(Data,coordNames=NULL,
                   hierLevels=c("1","2i","2"))
  
  summary(Data)
  
  # write Data to a csv file
  options(digits = 15)
  
  write.csv(Data, row.names=F, file=paste0('../dat/diary/',whale,'/processed_data_',whale,'.csv'))
  
  rawData$Time <- as.numeric(rawData$Time)
  write.csv(rawData, row.names=F, file=paste0('../dat/diary/',whale,'/processed_rawData_',whale,'.csv'))
}

Prep_Data("A100")
Prep_Data("A113")
Prep_Data("D21")
Prep_Data("D26")
Prep_Data("I107")
Prep_Data("I129")
Prep_Data("I145")
Prep_Data("L87")
Prep_Data("L88")
Prep_Data("R48")
Prep_Data("R58")