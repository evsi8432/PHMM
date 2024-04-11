library(MASS)

Data_fine_LDA <- Data_fine[Data_fine$divenum %in% train_dives,]
rawDataTrain <- rawData[rawData$segnum %in% unique(Data_fine_LDA$segnum),]

# get ceptrum and time
cs <- matrix(nrow=0,ncol=200)
segnums <- c()
IDs <- c()
labs <- c()

for(rownum in 1:nrow(Data_fine_LDA)){
  
  dive_df <- rawDataTrain[rawDataTrain$segnum %in% Data_fine_LDA$segnum[rownum] & 
                          rawDataTrain$ID %in% Data_fine_LDA$ID[rownum],]
  
  print(rownum/nrow(Data_fine_LDA))
  
  # get segnum and ID
  segnums <- c(segnums,Data_fine_LDA$segnum[rownum])
  IDs <- c(IDs,Data_fine_LDA$ID[rownum])
  
  # get cepstrum
  cept <- Mod(fft(dive_df$dyn.aw1)^2 + fft(dive_df$dyn.aw2)^2 + fft(dive_df$dyn.aw3)^2)
  cs <- rbind(cs,cept)
}

# make a dataframe
cept_df <- data.frame(segnum = Data_fine_LDA$segnum,
                      ID = IDs,
                      cept = cs)

write.csv(cept_df,"cept.csv")
