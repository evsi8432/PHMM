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

#whales = c("A100","A113",
#           "D21","D26",
#           "I107","I145",
#           "R48","R58")

options(digits.secs = 3) 

all_Data = NULL
all_rawData = NULL
all_ethogram = NULL

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
  rawData <- data.frame(fread(paste0('../dat/diary/',whale,'/processed_rawData_',whale,'.csv'),
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
  
  # convert roll
  Data$abs_roll <- abs(Data$roll)
  
  # pick columns to take averages for
  cols <- c("head","pitch","abs_roll",
            "htv","ptv","rtv","VeDBA",
            "w_low","w_high")
  
  # add averages for pitch, head, roll, VeDBA, ptv, htv, and rtv
  avg_df <- aggregate(Data[,cols], list(divenum = Data$divenum), 
                      mean, na.rm = TRUE)
  names(avg_df)[2:ncol(avg_df)] <- paste0("avg_",names(avg_df)[2:ncol(avg_df)])
  
  # add averages for the bottom
  Data[,paste0("bot_",cols)] <- Data[,cols] * Data[,"bottom"]
  
  bot_df <- aggregate(Data[,c(paste0("bot_",cols),"bottom")], 
                      list(divenum = Data$divenum), 
                      sum, na.rm = TRUE)
  
  bot_df[,paste0("bot_",cols)] <- bot_df[,paste0("bot_",cols)] / bot_df$bottom
  names(bot_df)[2:ncol(bot_df)] <- paste0("avg_",names(bot_df)[2:ncol(bot_df)])
  
  # now take jerk peak at bottom and over dive
  Data[,"bot_jp"] <- Data[,"jp"] * Data[,"bottom"]
  max_df <- aggregate(Data[,c("jp","bot_jp")], list(divenum = Data$divenum), 
                      max, na.rm = TRUE)
  names(max_df)[2:ncol(max_df)] <- paste0("max_",names(max_df)[2:ncol(max_df)])
  
  max_df[max_df$max_jp == -Inf,"max_jp"] <- NA
  max_df[max_df$max_bot_jp == -Inf,"max_bot_jp"] <- NA
  
  ### reduce Data to only dives ###
  
  Data <- Data[Data$level == 1,]
  
  Data <- merge(x = Data, y = avg_df, by = "divenum", all.x = TRUE)
  Data <- merge(x = Data, y = bot_df, by = "divenum", all.x = TRUE)
  Data <- merge(x = Data, y = max_df, by = "divenum", all.x = TRUE)
  
  Data <- Data[,colSums(is.na(Data))<nrow(Data)]
  
  # add labels
  ac_ethogram <- read.csv("../dat/ethograms/echolocation_ethogram.csv")
  ac_ethogram <- ac_ethogram[ac_ethogram$Animal == whale,]
  
  if (nrow(ac_ethogram > 0)){
    ac_ethogram$stime <- as.POSIXct(paste(date0,ac_ethogram$time),format="%Y-%m-%d %H:%M:%OS")
    ac_ethogram$etime <- ac_ethogram$stime + ac_ethogram$duration.sec
  }
  
  ethogram <- read_excel(paste0('../dat/diary/',whale,'/',whale,'.xlsx'))
  
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
  
  # get dives when cam is off
  camDives <- unique(rawData[rawData$camon==1,"divenum"])
  NAdives <- setdiff(1:nrow(Data),camDives)
  
  ### initialize other labels ###
  
  # WARNING: we may be preferentially sampling "TRUE" here, but IDK since
  # most trues have longer dives which are more likely to have some cam missing
  # so we can't just throw out dives without full cam footage
  
  Data$foraging <- FALSE
  Data$scales <- FALSE
  Data$clicks <- FALSE
  Data$click.train <- FALSE
  Data$buzz <- FALSE
  Data$crunch <- FALSE
  Data$whales <- FALSE
  Data$touch <- FALSE
  Data$kelp <- FALSE
  
  ### Label Foraging Success indicators ###
  
  forage_grace_period = 0*60 # seconds
  scales_grace_period = 5*60 # seconds
  
  # go through regular ethogram
  for (i in seq_along(ethogram$stime)){
    
    # start with foraging labels with grace period
    if (ethogram$Behavior[i] %in% c('Kill','Confirmed Fish Capture')){
      inds <- max(which((Data$stime <= ethogram$etime[i]) & 
                          ((Data$etime+forage_grace_period) >= ethogram$stime[i]) &
                          (Data$broadDiveType %in% TRUE)))
      print(inds)
      print(ethogram$Behavior[i])
      Data$foraging[inds] <- TRUE
    }
    
    # start with foraging labels with grace period
    if (ethogram$Behavior[i] %in% c("Scales")){
      inds <- max(which((Data$stime <= ethogram$etime[i]) & 
                          ((Data$etime+scales_grace_period) >= ethogram$stime[i]) &
                          (Data$broadDiveType %in% TRUE)))
      print(inds)
      print(ethogram$Behavior[i])
      Data$scales[inds] <- TRUE
    }
    
    # then do whales
    if (ethogram$Behavior[i] %in% c('whale')){
      inds <- which((Data$stime <= ethogram$etime[i]) & 
                      (Data$etime >= ethogram$stime[i]))
      Data$whales[inds] <- TRUE
    }
    
    # then do touch
    if (ethogram$Behavior[i] %in% c('touch','rub')){
      inds <- which((Data$stime <= ethogram$etime[i]) & 
                      (Data$etime >= ethogram$stime[i]))
      Data$touch[inds] <- TRUE
    }
    
    # finally do kelp
    if (ethogram$Behavior[i] %in% c('kelp')){
      inds <- which((Data$stime <= ethogram$etime[i]) & 
                      (Data$etime >= ethogram$stime[i]))
      Data$kelp[inds] <- TRUE
    }  
  }
  
  
  # go through acoustic ethogram
  for (i in seq_along(ac_ethogram$stime)){
    
    # first do clicks
    if (ac_ethogram$broad.type[i] %in% c('clicks')){
      inds <- which((Data$stime <= ac_ethogram$etime[i]) & 
                      (Data$etime >= ac_ethogram$stime[i]))
      Data$clicks[inds] <- TRUE
    }
    
    # now do click trains
    if (ac_ethogram$broad.type[i] %in% c('click train')){
      inds <- which((Data$stime <= ac_ethogram$etime[i]) & 
                      (Data$etime >= ac_ethogram$stime[i]))
      print(inds)
      Data$click.train[inds] <- TRUE
    }
    
    # now do buzzes
    if (ac_ethogram$broad.type[i] %in% c('buzz')){
      inds <- which((Data$stime <= ac_ethogram$etime[i]) & 
                      (Data$etime >= ac_ethogram$stime[i]))
      Data$buzz[inds] <- TRUE
    }
    
    # Finally, do crunches with grace period
    if (ac_ethogram$broad.type[i] %in% c('crunch')){
      inds <- max(which((Data$stime <= ac_ethogram$etime[i]) & 
                          ((Data$etime+forage_grace_period) >= ac_ethogram$stime[i]) &
                          Data$broadDiveType %in% TRUE))
      print(inds)
      print(ac_ethogram$broad.type[i])
      Data$foraging[inds] <- TRUE
    }  
  }
  
  # make dives with no camera on "NA". This may label some dives as FALSE when
  # they should be TRUE if a behaviour occured when the cam was off.
  
  Data[NAdives,"whales"] <- NA
  Data[NAdives,"touch"] <- NA
  Data[NAdives,"scales"] <- NA
  Data[NAdives,"kelp"] <- NA
  
  # get rid of click, click.train, buzz, and crunch if there is no echolocation
  if (nrow(ac_ethogram) == 0){
    Data[,"clicks"] <- NA
    Data[,"click.train"] <- NA
    Data[,"buzz"] <- NA
    Data[,"crunch"] <- NA
    Data[NAdives,"foraging"] <- NA
  }
  
  # preprare data somehow idk why this is needed
  Data <- prepData(Data,coordNames=NULL)
  
  ### convert columns ###
  Data$maxDepth <- log(Data$maxDepth)
  Data$diveDuration <- log(Data$diveDuration)
  Data$postDiveInt <- log(Data$postDiveInt)
  Data$avg_bot_abs_roll <- log(Data$avg_bot_abs_roll)
  Data$avg_bot_htv <- log(Data$avg_bot_htv)
  Data$max_bot_jp <- log(Data$max_bot_jp)
  
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
  all_rawData <- rbind(all_rawData,rawData)
  all_ethogram <- rbind(all_ethogram,ethogram)
  
  # display results
  print(max(c(0,all_rawData$divenum),na.rm=T))
}

# add sex to Data
Data$Sex <- "Female"
Data[Data$ID %in% c("I107","D21","L87","L88"),"Sex"]  <- "Male"
Data$Sex <- factor(Data$Sex)

# adjust/add columns
Data$postDiveInt <- exp(Data$postDiveInt) + runif(nrow(Data))/2
Data$maxDepth <- exp(Data$maxDepth)
Data$maxDepth[Data$ID == "D21"] <- Data$maxDepth[Data$ID == "D21"] + 0.5
Data$diveDuration <- exp(Data$diveDuration) + runif(nrow(Data))/2
Data$logMaxDepth <- log(Data$maxDepth)
Data$logDiveDuration <- log(Data$diveDuration)
Data$logWLow <- log(Data$avg_w_low)
Data$logWHigh <- log(Data$avg_w_high)
Data$logWTotal <- log(Data$avg_w_low + Data$avg_w_high)
Data$logMDDD.x <- Data$logMaxDepth
Data$logMDDD.y <- Data$logDiveDuration
Data$logW.x <- Data$logWLow
Data$logW.y <- Data$logWHigh

# add labels
Data$knownState <- NA
Data$knownState[Data$diveType %in% "Resting"] <- 1
Data$knownState[Data$diveType %in% "Travelling"] <- 2
Data$knownState[Data$diveType %in% "Foraging"] <- 3
Data$knownState[is.na(Data$knownState)] <- 4
Data$knownState <- as.factor(Data$knownState)

# save results
write.csv(all_Data, row.names=F, file=paste0('../../dat/Final_Data_Beth.csv'))
write.csv(all_rawData, row.names=F, file=paste0('../../dat/Final_rawData_Beth.csv'))
write.csv(all_ethogram, row.names=F, file=paste0('../../dat/Final_ethogram_Beth.csv'))
