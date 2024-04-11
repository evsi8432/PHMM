# add label to rawData
rawData$label <- NA
Data_fine0 <- Data_fine[Data_fine$knownState %in% 4,
                        c("stime","etime","foraging_signs")]

for(seg_ind in 1:nrow(Data_fine0)){
  stime <- Data_fine0$stime[seg_ind]
  etime <- Data_fine0$etime[seg_ind]
  foraging_sign <- Data_fine0$foraging_signs[seg_ind]
  rawData$label[rawData$Time >= stime & rawData$Time <= etime] <- T
  rawData$label[rawData$Time >= stime & rawData$Time <= etime] <- T
}

# plot features
ggplot(Data_fine[Data_fine$segnum != 0,]) +
  geom_histogram(aes(x = logWLow, y = ..density.., 
                     fill = knownState %in% 4),
                 alpha = 0.5, position = "identity")

ggplot(Data_fine[Data_fine$segnum != 0,]) +
  geom_histogram(aes(x = logJpNorm, y = ..density.., 
                     fill = knownState %in% 4),
                 alpha = 0.5, position = "identity")

ggplot(Data_fine[Data_fine$segnum != 0,]) +
  geom_histogram(aes(x = logWHigh, y = ..density.., 
                     fill = knownState %in% 4),
                 alpha = 0.5, position = "identity")

ggplot(Data_fine[Data_fine$segnum != 0,]) +
  geom_histogram(aes(x = logHtv, y = ..density.., 
                     fill = knownState %in% 4),
                 alpha = 0.5, position = "identity")

ggplot(Data_fine[Data_fine$segnum != 0,]) +
  geom_histogram(aes(x = delt_d, y = ..density.., 
                     fill = knownState %in% 4),
                 alpha = 0.5, position = "identity")

ggplot(Data_fine[Data_fine$segnum != 0,]) +
  geom_histogram(aes(x = rajp, y = ..density.., 
                     fill = knownState %in% 4),
                 alpha = 0.5, position = "identity")

ggplot(Data_fine[Data_fine$knownState %in% c(4,5),]) +
  geom_point(aes(x = logWLow, y = logHtv, color = knownState))

# plot dives
for(dive in pos_dives){
  print(dive)
  dive_df <- rawData[rawData$divenum %in% dive,]
  start_ascent <- max(dive_df$Time[dive_df$elev < 0.7*min(dive_df$elev)])
  dive_df <- dive_df[(dive_df$Time > start_ascent - 60) & (dive_df$Time < start_ascent+60),]
  dive_df$label[is.na(dive_df$label)] <- "NA"
  dive_df$abs_rajp <- abs(dive_df$rajp)
  
  dive_df0 <- dive_df %>%
    pivot_longer(cols = c("elev","jp_normed","w_high","w_low"),
                 names_to = "feature")
  
  print(ggplot(dive_df0[dive_df0$segnum != 0,]) +
          geom_line(aes(x = Time, y = value, 
                        color = label, group = segnum)) +
          facet_wrap(~feature, ncol = 1, scales = "free_y") +
          ggtitle(toString(dive)))
  
  dive_df1 <- dive_df %>%
    pivot_longer(cols = c("elev","dyn.aw1","dyn.aw2","dyn.aw3"),
                 names_to = "feature")
  
  print(ggplot(dive_df1[dive_df1$segnum != 0,]) +
          geom_line(aes(x = Time, y = value, 
                        color = label, group = segnum)) +
          facet_wrap(~feature, ncol = 1, scales = "free_y") +
          ggtitle(toString(dive)))
}
