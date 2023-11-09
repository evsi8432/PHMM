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

### Load Data ###
setwd("~/Documents/Research/PHMM/src")
whales = c("A100","A113",
           "D21","D26",
           "I107","I129","I145",
           "L87","L88",
           "R48","R58")

options(digits.secs = 3) 

all_Data = NULL
#all_rawData = NULL
#all_ethogram = NULL

for (whale in whales){
  
  print(whale)
  
  Data <- data.frame(fread(paste0('../../dat/diary/',whale,'/processed_Data_',whale,'.csv'),
                           colClasses=c("level"="character")))
  
  ### prepare data for HHMM ###
  Data$stime <- as.POSIXct(Data$stime, origin = '1970-01-01', tz = "America/Vancouver")
  Data$etime <- as.POSIXct(Data$etime, origin = '1970-01-01', tz = "America/Vancouver")
  
  # we have to do this maybe?
  Data <- prepData(Data, coordNames = NULL)
  
  # get the raw data
  rawData <- data.frame(fread(paste0('../../dat/diary/',whale,'/processed_rawData_',whale,'.csv'),
                              colClasses=c("Time"="double")))
  rawData$Time <- as.POSIXct(rawData$Time, origin = '1970-01-01', tz = "America/Vancouver")
  rawData$ID <- whale
  rawData <- rawData[,!(names(rawData) %in% c("vstate1","vstate2","label"))]
  
  date0 = format(rawData$Time[1], "%Y-%m-%d")  # prone to error for overnight data
  
  # add dive number 
  Data$divenum <- NA
  Data$divenum[Data$level == 1] <- (1:sum(Data$level == 1))
  Data <- Data %>% fill(divenum)
  Data$divenum <- Data$divenum + max(c(0,all_Data$divenum),na.rm=T)
  rawData$divenum <- rawData$divenum + max(c(0,all_rawData$divenum),na.rm=T)
  
  # add seg number
  Data$segnum <- NA
  Data$segnum[Data$level == 2] <- (1:sum(Data$level == 2))
  
  # convert roll
  Data$abs_roll <- abs(Data$roll)
  
  # add labels
  ethogram <- read_excel(paste0('../../dat/diary/',whale,'/',whale,'.xlsx'))
  
  media.file = names(ethogram)[4]
  temp <- strsplit(ethogram[[media.file]],"-")
  start_inds <- lapply(lapply(temp,function(x){grepl("20200",x)}),which)
  temp <- mapply(function(X,Y){paste0(X[Y],"-",X[Y+1])},temp,start_inds)
  
  ethogram$time <- as.POSIXct(temp,format="%Y%m%d-%H%M%S")
  ethogram$`Start (s)` <- as.numeric(ethogram$`Start (s)`)
  ethogram$`Stop (s)` <- as.numeric(ethogram$`Stop (s)`)
  ethogram$`Duration (s)` <- as.numeric(ethogram$`Duration (s)`)
  
  ethogram$stime <- ethogram$time + ethogram$`Start (s)`
  ethogram$etime <- ethogram$time + ethogram$`Stop (s)`
  
  # get segments when cam is off
  camchange_inds <- c(1,which(diff(rawData$camon) != 0) + 1)
  camchange_df <- rawData[camchange_inds,c("Time","camon")]
  Data$camon <- NA
  for (i in 1:nrow(camchange_df)) {
    inds <- (Data$stime >= camchange_df[i,"Time"]) & !is.na(Data$stime)
    Data$camon[inds] <- camchange_df[i,"camon"]
  }
  Data$camon[Data$level != 2] <- NA
  
  ### initialize other labels ###
  Data$foraging <- FALSE # confirmed prey, kill
  Data$chase <- FALSE
  Data$scales <- FALSE
  Data$echo.steady <- FALSE
  Data$echo.rapid <- FALSE
  Data$crunch <- FALSE
  
  ### Label Foraging Success indicators ###
  for (i in seq_along(ethogram$stime)){
    
    # foraging labels
    if (ethogram$Behavior[i] %in% c('Kill','Confirmed Fish Capture')){
      inds <- which((Data$stime <= ethogram$etime[i]) & 
                    (Data$etime >= ethogram$stime[i]))
      print(inds)
      print(ethogram$Behavior[i])
      print(ethogram$stime[i])
      Data$foraging[inds] <- TRUE
    }
    
    # chase labels
    if (ethogram$Behavior[i] %in% c("chase")){
      inds <- which((Data$stime <= ethogram$etime[i]) & 
                    (Data$etime >= ethogram$stime[i]))
      print(inds)
      print(ethogram$Behavior[i])
      Data$chase[inds] <- TRUE
    }
    
    # scales labels
    if (ethogram$Behavior[i] %in% c("Scales")){
      inds <- which((Data$stime <= ethogram$etime[i]) & 
                    (Data$etime >= ethogram$stime[i]))
      print(inds)
      print(ethogram$Behavior[i])
      Data$scales[inds] <- TRUE
    }
    
    # echolocation labels
    if (ethogram$Behavior[i] %in% c("Echolocation")){
      inds <- which((Data$stime <= ethogram$etime[i]) & 
                    (Data$etime >= ethogram$stime[i]))
      print(inds)
      print(ethogram$Behavior[i])
      print(ethogram$Modifiers[i])
      
      if (ethogram$Modifiers[i] == "Steady") {
        Data$echo.steady[inds] <- TRUE
      } else if (ethogram$Modifiers[i] == "Rapid"){
        Data$echo.rapid[inds] <- TRUE
      }
    }
    
    # crunch labels
    if (ethogram$Behavior[i] %in% c('crunch')){
      inds <- which((Data$stime <= ethogram$etime[i]) & 
                    (Data$etime >= ethogram$stime[i]))
      Data$crunch[inds] <- TRUE
    }
  }
  
  # make dives with no camera on "NA". This may label some dives as FALSE when
  # they should be TRUE if a behaviour occurred when the cam was off.
  Data[!(Data$camon %in% T),"foraging"] <- NA
  Data[!(Data$camon %in% T),"chase"] <- NA
  Data[!(Data$camon %in% T),"scales"] <- NA
  Data[!(Data$camon %in% T),"echo.steady"] <- NA
  Data[!(Data$camon %in% T),"echo.rapid"] <- NA
  Data[!(Data$camon %in% T),"crunch"] <- NA
  
  # change ethograms to be consistent
  if ("Observation id...1" %in% colnames(ethogram)){
    ethogram$`Observation id` <- ethogram$`Observation id...1`
  }
  if ("Observation date...2" %in% colnames(ethogram)){
    ethogram$`Observation date` <- ethogram$`Observation date...2`
  }
  if ("Description...3" %in% colnames(ethogram)){
    ethogram$Description <- ethogram$`Description...3`
  }
  if ("Media file...4" %in% colnames(ethogram)){
    ethogram$`Media file` <- ethogram$`Media file...4`
  }
  if ("Total length...5" %in% colnames(ethogram)){
    ethogram$`Total length` <- ethogram$`Total length...5`
  }
  
  if (!is.null(all_ethogram)){
    cols_to_keep <- intersect(colnames(ethogram),colnames(all_ethogram))
    ethogram <- ethogram[,c(cols_to_keep)]
    all_ethogram <- all_ethogram[,c(cols_to_keep)]
  }
  
  # concatenate Data
  all_Data <- rbind(all_Data,Data)
  #all_rawData <- rbind(all_rawData,rawData)
  #all_ethogram <- rbind(all_ethogram,ethogram)
  
  # display results
  #print(max(c(0,all_rawData$divenum),na.rm=T))
}

# add sex to Data
all_Data$Sex <- "Female"
all_Data[all_Data$ID %in% c("I107","D21","L87","L88"),"Sex"]  <- "Male"
all_Data$Sex <- factor(all_Data$Sex)

# adjust/add columns
all_Data$postDiveInt <- all_Data$postDiveInt + runif(nrow(all_Data))/2
all_Data$maxDepth[all_Data$ID == "D21"] <- all_Data$maxDepth[all_Data$ID == "D21"] + 1
all_Data$diveDuration <- all_Data$diveDuration + runif(nrow(all_Data))/2
all_Data$logMaxDepth <- log(all_Data$maxDepth)
all_Data$logDiveDuration <- log(all_Data$diveDuration)
all_Data$logWLow <- log(all_Data$w_low)
all_Data$logWHigh <- log(all_Data$w_high)
all_Data$logWTotal <- log(all_Data$w_low + all_Data$w_high)
all_Data$logMDDD.x <- all_Data$logMaxDepth
all_Data$logMDDD.y <- all_Data$logDiveDuration
all_Data$logW.x <- all_Data$logWLow
all_Data$logW.y <- all_Data$logWHigh

# add labels
all_Data$knownState <- NA
all_Data$knownState[all_Data$diveType %in% "Resting"] <- 1
all_Data$knownState[all_Data$diveType %in% "Travelling"] <- 2
all_Data$knownState[all_Data$diveType %in% "Foraging"] <- 3
all_Data$knownState[is.na(all_Data$knownState)] <- 4
all_Data$knownState <- as.factor(all_Data$knownState)

# save results
write.csv(all_Data, row.names=F, file=paste0('../../dat/Final_Data_fine.csv'))
#write.csv(all_rawData, row.names=F, file=paste0('../../dat/Final_rawData_fine.csv'))
