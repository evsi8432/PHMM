library(momentuHMM)
library(dplyr)
library(mclust)
library(data.table)

# Coarse Data options
sex <- c("Male","Female")
dd_thresh <- 2
md_thresh <- 0.5
md_threshs <- c(5,10,30,50)
hier = F

data_file <- '../../dat/Final_Data_Beth.csv'

# load data
Data <- data.frame(fread(data_file))

# mark the data we are going to drop
Data$keep <- TRUE
Data$keep <- Data$keep & (Data$Sex %in% sex)
Data$keep <- Data$keep & (Data$diveDuration >= dd_thresh)
Data$keep <- Data$keep & (Data$maxDepth >= md_thresh)
Data$keep <- Data$keep & !is.na(Data$ID)

# assign labels of dives that have been deleted to nearest dive
for(i in c()){#1:nrow(Data)){
  if((Data[i,"knownState"] != 4) & (Data[i,"diveDuration"] < dd_thresh)){
    
    close_dives <- Data[(Data$stime < Data[i,"etime"]+300) & # 5 minutes
                        (Data$etime > Data[i,"stime"]-300) & # 5 minutes
                        (Data$ID == Data[i,"ID"]) &
                        (Data$keep),]
    
    if(nrow(close_dives) > 0){
      close_dives$dist <- pmin(abs(close_dives$etime - Data[i,"stime"]),
                               abs(close_dives$stime - Data[i,"etime"]))
      
      dive_ind <- as.numeric(rownames(close_dives)[which.min(close_dives$dist)])
      
      if(Data[dive_ind,"knownState"] == 4){
        Data[dive_ind,"knownState"] <- Data[i,"knownState"]
      }
    }
  }
}

Data$label <- factor(Data$knownState, levels = 1:4)

# get values for D21 (max Depth is off for some reason)
D21_data <- read.csv("../../dat/dive_summaries/D21_diveMove_summary_0.5m.csv")
Data$maxDepth[Data$ID %in% "D21"] <- D21_data$max.dive.depth.m
Data$logMaxDepth[Data$ID %in% "D21"] <- log(D21_data$max.dive.depth.m)
Data$logMDDD.x[Data$ID %in% "D21"] <- log(D21_data$max.dive.depth.m)

# Assign categorical variables
Data$maxDepthCat <- 1
for(thresh in md_threshs){
  Data$maxDepthCat[Data$maxDepth > thresh] <- Data$maxDepthCat[Data$maxDepth > thresh] + 1
}

# drop the rows we said we were going to drop
Data <- Data[Data$keep,!(names(Data) %in% "keep")]

# divide Data ID into two for each whale
for(ID in unique(Data$ID)){
  nlabels <- cumsum(Data$knownState %in% c(1,2,3) & Data$ID == ID)
  split_ind <- which(nlabels > tail(nlabels,1) / 2)[1]
  Data$ID[(1:nrow(Data) < split_ind) & (Data$ID == ID)] <- paste0(ID,"a")
  Data$ID[(1:nrow(Data) >= split_ind) & (Data$ID == ID)] <- paste0(ID,"b")
}

### change to hierarchical ###
if(hier){
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
      
      # turn dive-level label to bout level label
      if(forg){
        tmp2$knownState[1] <- 3
        tmp2$knownState[-1] <- 4
      }
      else if (rest + trav > 1){
        tmp2$knownState <- 4
      }
      else if (rest){
        tmp2$knownState[1] <- 1
        tmp2$knownState[-1] <- 4
      }
      else if (trav){
        tmp2$knownState[1] <- 2
        tmp2$knownState[-1] <- 4
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
      if(nrow(tmp2) > 0){
        time <- tail(tmp2,n=1)[1,"etime"]
      } else {
        time <- time + span*60
      }
    }
  }
  Data <- Data0
  rownames(Data) <- 1:nrow(Data)
  Data$label <- factor(Data$label,levels = 1:7)
}
