library(momentuHMM)
library(dplyr)
library(mclust)
library(data.table)

# load data
Data <- data.frame(fread('../../dat/Final_Data_Beth.csv'))

# assign labels of dives to dive before
Data$pseudolabel <- FALSE
for(i in 1:nrow(Data)){
  if((Data[i,"knownState"] != 4) & (Data[i,"diveDuration"] < dd_thresh)){
    past_dives <- Data[(Data$stime < Data[i,"stime"]) &
                         (Data$etime > Data[i,"stime"]-300) &
                         (Data$ID == Data[i,"ID"]) &
                         (Data$diveDuration >= dd_thresh),]
    if(nrow(past_dives) >= 1){
      dive_ind <- tail(rownames(past_dives),1)
      if(Data[dive_ind,"knownState"] == 4){
        Data[dive_ind,"knownState"] <- Data[i,"knownState"]
        Data[dive_ind,"pseudolabel"] <- TRUE
      }
    }
  }
}

Data$label <- (2*Data$knownState - 1) + Data$pseudolabel
Data$label <- as.factor(Data$label)

# only keep a subset of the entire data set
Data <- Data[Data$Sex %in% sex,]
Data <- Data[Data$diveDuration >= dd_thresh,]
Data <- Data[Data$maxDepth >= md_thresh,]
Data <- Data[!is.na(Data$ID),]