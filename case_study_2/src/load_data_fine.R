library(momentuHMM)
library(dplyr)
library(mclust)
library(data.table)

#data_file <- '../../../dat/Final_Data_Beth.csv'
data_file <- '../../dat/Final_Data_fine.csv'

# load data
Data <- data.frame(fread(data_file))

# mark the data we are going to drop
Data$keep <- TRUE
#Data$keep <- Data$keep & (Data$Sex %in% sex)

# drop the rows we said we were going to drop
Data <- Data[Data$keep,!(names(Data) %in% "keep")]

# add labels
Data$label <- factor(Data$knownState, levels = 1:4)

# get normalized jerk peak readings
rawData <- data.frame(fread('../../dat/Final_rawData_Beth.csv'))

# calculate jerk
rawData$elev <- -rawData$p
rawData$jerk <- c(NA,50*sqrt(diff(rawData$Aw_1)^2 + 
                             diff(rawData$Aw_2)^2 + 
                             diff(rawData$Aw_3)^2))

# add max depth of each dive
max_depths <- rawData[,c("divenum","p")] %>% 
  group_by(divenum) %>% 
  slice_max(p,with_ties = F)
names(max_depths) <- c("divenum","maxDepth")
rawData <- left_join(x = rawData, y = max_depths)

# add bottom times of each dive
rawData0 <- rawData[rawData$p > 0.7*rawData$maxDepth,c("divenum","Time")]
bottom_times <- rawData0 %>% 
  group_by(divenum) %>%
  summarize(sbot = min(Time), 
            ebot = max(Time))
rawData <- left_join(x = rawData, y = bottom_times)

# add median jerk of each dive
rawDataBottom <- rawData[rawData$Time > rawData$sbot & rawData$Time < rawData$ebot,
                         c("divenum","jerk")]

med_jerk <- rawDataBottom %>% 
  group_by(divenum) %>%
  summarize(med_jerk = median(jerk, na.rm = T))

rawData <- left_join(x = rawData, y = med_jerk)
Data <- left_join(x = Data, y = med_jerk)

# divide jerk by median jerk of each dive
Data$jp_normed <- Data$jp / Data$med_jerk
rawData$jp_normed <- rawData$jp / rawData$med_jerk

# store as Data_fine
Data_fine <- Data
Data_fine <- Data_fine[Data_fine$level == 2,]

# add features
Data_fine$delt_d <- c(NA,diff(Data_fine$ad))
Data_fine$logHtv <- log(Data_fine$htv)
Data_fine$logJpNorm <- log(Data_fine$jp_normed)

# save the transformed dataset
write.csv(Data_fine,"../../dat/Final_Data_fine1.csv")
