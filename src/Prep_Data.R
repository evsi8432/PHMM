library(momentuHMM)
library(circular)
library(CircStats)
library(runner)
library(ggplot2)
library(tidyr)

setwd("~/Documents/Research/PHMM/src")

# define summary stats function
SumStats <- function(x){
  
  # get Jerk Peak 
  jerk <- sqrt(diff(x$Aw_1)^2 + diff(x$Aw_2)^2 + diff(x$Aw_3)^2)
  jp <- max(jerk,na.rm=TRUE)
  
  # find roll at Jerk Peak
  jpind <- which(jerk == jp) 
  rajp <- x$roll[jpind]
  
  # find heading variance
  hv <- circ.disp(x$roll)$var
  
  # get average depth
  ad <- mean(x$p)
  
  # get average acceleration
  aw1 <- mean(x$Aw_1)
  aw2 <- mean(x$Aw_2)
  aw3 <- mean(x$Aw_3)
  
  # get wiggliness
  w <- Mod(fft(x$Aw_1))^2 + Mod(fft(x$Aw_2))^2 + Mod(fft(x$Aw_3))^2
  w <- sum(w[2:10])
  
  # get time
  stime <- min(x$Time)
  etime <- max(x$Time)
  
  # return everything
  return(c(jp,rajp,hv,ad,aw1,aw2,aw3,w,stime,etime))
}

# load in data for whale of interest
whale = "I107" # I145
date0 = "2020-08-25" # "2020-08-30
rawData <- read.csv(paste0("../dat/diary/kinematic_data_calibrated_",
                           whale,"_2020_50Hz.csv"))

# convert to a date time
rawData$Time <- as.POSIXct(rawData$Time,format="%d-%b-%Y %H:%M:%OS")

# label the dives by dive number
surfaceThresh <- 0.25 # threshold to begin dive in meters
divedurationThresh <- 2 # minimum duration for a dive in seconds
h <- 100 # number of observations in a moving window

rawData$divenum <- 0
segnum <- 0
rawData$segnum <- segnum
rawDataAtSurface <- rawData[rawData$p < surfaceThresh,]
inds <- which(as.double(diff(rawDataAtSurface$Time),units='secs') > 2)
divestarts <- rawDataAtSurface$Time[inds]

# define sum stats
rawData$jp <- NA
rawData$rajp <- NA
rawData$hv <- NA
rawData$w <- NA

# define a full data matrix
Data <- NULL

for (i in 1:(length(divestarts))){
  
  print(i/length(divestarts))
  
  # get dive start / end
  startind <- which(rawData$Time == divestarts[i])
  if (i == length(divestarts)){
    endind <- nrow(rawData)
  } else {
    endind <- which(rawData$Time == divestarts[i+1])
  }
  
  # label the rawData
  rawData$divenum[startind:endind] <- i
  
  # get dataframe with dives
  diveData <- rawData[startind:endind,]
  
  # get coarse-scale data
  dd <- as.double(difftime(max(diveData$Time),min(diveData$Time)),units='secs')
  md <- max(diveData$p)
  
  coarseScale <- data.frame(ID = whale,
                            level="1",
                            diveDuration = dd,
                            maxDepth = md,
                            jp = NA,
                            rajp = NA,
                            hv = NA,
                            ad = NA,
                            aw1 = NA,
                            aw2 = NA,
                            aw3 = NA,
                            w = NA,
                            stime = as.numeric(rawData$Time[startind]),
                            etime = as.numeric(rawData$Time[endind]))
  
  # get moving window values
  RawFineScale <- runner(
    x = diveData,
    k = h,
    at = seq(h+1,nrow(diveData),h),
    f = SumStats
  )
  
  print(dim(RawFineScale)[2])
  
  # label the segment numbers in the rawData
  for (j in 1:length(RawFineScale[9,])){
    segnum = segnum + 1
    segstartind <- which(rawData$Time == RawFineScale[9,j])
    segendind <- which(rawData$Time == RawFineScale[10,j])
    rawData$segnum[segstartind:segendind] <- segnum
    rawData$jp[segstartind:segendind] <- RawFineScale[1,j]
    rawData$rajp[segstartind:segendind] <- RawFineScale[2,j]
    rawData$hv[segstartind:segendind] <- RawFineScale[3,j]
    rawData$w[segstartind:segendind] <- RawFineScale[8,j]
  }
  
  
  fineScale <- data.frame(ID = whale,
                          level="2",
                          diveDuration = NA,
                          maxDepth = NA,
                          jp = RawFineScale[1,],
                          rajp = RawFineScale[2,],
                          hv = RawFineScale[3,],
                          ad = RawFineScale[4,],
                          aw1 = RawFineScale[5,],
                          aw2 = RawFineScale[6,],
                          aw3 = RawFineScale[7,],
                          w = RawFineScale[8,],
                          stime = RawFineScale[9,],
                          etime = RawFineScale[10,])
  
  fineScale <- rbind(data.frame(ID = whale,
                                level="2i",
                                diveDuration = NA,
                                maxDepth = NA,
                                jp = NA,
                                rajp = NA,
                                hv = NA,
                                ad = NA,
                                aw1 = NA,
                                aw2 = NA,
                                aw3 = NA,
                                w = NA,
                                stime = NA,
                                etime = NA),
                     fineScale)
  
  Data <- rbind(Data,coarseScale,fineScale)
}

# change etime and stime to POSIXct

# log-transform the variables
Data$jp <- log(Data$jp)
Data$hv <- logit(Data$hv)

summary(Data)

# add labels
ethogram <- read.csv("../dat/ethograms/echolocation_ethogram.csv")
ethogram <- ethogram[ethogram$Animal == whale,]
ethogram$stime <- as.POSIXct(paste(date0, ethogram$time))
ethogram$stime <- as.numeric(ethogram$stime)
ethogram$etime <- ethogram$stime + ethogram$duration.sec

Data$label <- NA

for (i in 1:length(ethogram$broad.type)){
  inds <- which((Data$stime < ethogram$etime[i]) & 
                (Data$etime > ethogram$stime[i]) & 
                (Data$level != 1))
  print(ethogram$broad.type[i])
  print(inds)
  #print(as.POSIXct(Data$etime[inds],origin = '1970-01-01'))
  Data$label[inds] <- ethogram$broad.type[i]
}

# preprare data somehow idk why this is needed
Data <- prepData(Data,coordNames=NULL,
                 hierLevels=c("1","2i","2"))

summary(Data)

# write Data to a csv file
write.csv(Data, row.names=F, file=paste0('../Dat/diary/processed_data_',whale,'.csv'))

rawData$Time <- as.numeric(rawData$Time)
write.csv(rawData, row.names=F, file=paste0('../Dat/diary/processed_rawData_',whale,'.csv'))
