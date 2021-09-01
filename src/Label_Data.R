library(momentuHMM)
library(CircStats)
library(runner)

setwd("~/Documents/Research/PHMM/src")

ethogram <- read.csv("../dat/ethograms/echolocation_ethogram.csv")
ethogram <- ethogram[ethogram$Animal == 'I145',]
ethogram$stime <- as.POSIXct(paste("2020-08-30", ethogram$time))
ethogram$etime <- ethogram$stime + seconds(ethogram$duration.sec)



Data$label <- NA

for (i in 1:length(ethogram)){
  inds <- which((Data$stime < ethogram$etime[i]) & 
                (Data$etime > ethogram$stime[i]) & 
                (Data$level != 1))
  Data$label[inds] <- ethogram$broad.type[i]
}
# convert to a date time
rawData$Time <- as.POSIXct(rawData$Time,format='%d-%b-%Y %H:%M:%OS')

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
  
  coarseScale <- data.frame(ID = "I145",
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
                            stime = rawData$Time[startind],
                            etime = rawData$Time[endind])
  
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
  
  
  fineScale <- data.frame(ID = "I145",
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
                          stime = as.POSIXct(RawFineScale[9,],
                                             format='%d-%b-%Y %H:%M:%OS',
                                             origin = origin),
                          etime = as.POSIXct(RawFineScale[10,],
                                             format='%d-%b-%Y %H:%M:%OS',
                                             origin = origin))
  
  fineScale <- rbind(data.frame(ID = "I145",
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
                                stime = as.POSIXct(NA),
                                etime = as.POSIXct(NA)),
                     fineScale)
  
  Data <- rbind(Data,coarseScale,fineScale)
}

# log-transform the variables
Data$jp <- log(Data$jp)
Data$hv <- log(Data$hv)

summary(Data)

# preprare data somehow idk why this is needed
Data <- prepData(Data,coordNames=NULL,
                 hierLevels=c("1","2i","2"))

summary(Data)

