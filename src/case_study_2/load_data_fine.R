library(momentuHMM)
library(dplyr)
library(mclust)
library(data.table)

#data_file <- '../../../dat/Final_Data_Beth.csv'
data_file <- '../../../dat/Final_Data_fine.csv'

# load data
Data <- data.frame(fread(data_file))

# mark the data we are going to drop
Data$keep <- TRUE
#Data$keep <- Data$keep & (Data$Sex %in% sex)

# drop the rows we said we were going to drop
Data <- Data[Data$keep,!(names(Data) %in% "keep")]

# add labels
Data$label <- factor(Data$knownState, levels = 1:4)
