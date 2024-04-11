rawDataPos <- rawData[rawData$divenum %in% train_dives,]
rawDataPos$elev <- -rawDataPos$p
rawDataPos$logWLow <- log(rawDataPos$w_low)
rawDataPos$dLogWlow <- c(0,diff(rawDataPos$logWLow))
rawDataPos$d2LogWlow <- c(0,diff(rawDataPos$dWlow))

# add labels to raw data frame
rawDataPos$label <- NA
Data_fine0 <- Data_fine[Data_fine$foraging_signs %in% "crunch",
                        c("stime","etime","foraging_signs")]

for(seg_ind in 1:nrow(Data_fine0)){
  stime <- Data_fine0$stime[seg_ind]
  etime <- Data_fine0$etime[seg_ind]
  foraging_sign <- Data_fine0$foraging_signs[seg_ind]
  rawDataPos$label[rawDataPos$Time >= stime & rawDataPos$Time <= etime] <- T
  rawDataPos$label[rawDataPos$Time >= stime & rawDataPos$Time <= etime] <- T
}

# window size
w0 <- 200/2 # half-window for data
w1 <- 500/2 # half-window for normalization
w2 <- 200/2 # half-window for labels 

# get ceptrum and time
ts <- c(rawDataPos$Time[1])
cs <- matrix(nrow=0,ncol=2*w0)
peak_jerk_norm <- c()
labs <- c()
dives0 <- c()

for(dive in dives){
  
  dive_df <- rawDataPos[rawDataPos$divenum == dive,]
    
  print(dive)
  print(min(dive_df$Time))
  print(max(dive_df$Time))
  print("")
  
  for(i in 1:(floor(nrow(dive_df)/w0))){
    
    ind <- w0 * (2*i - 1)
    
    # extract data
    x0 <- dive_df[max(1,ind-w0+1):min(nrow(dive_df),ind+w0),]
    x1 <- dive_df[max(1,ind-w1+1):min(nrow(dive_df),ind+w1),]
    x2 <- dive_df[max(1,ind-w2+1):min(nrow(dive_df),ind+w2),]
    
    # get normalization constants
    mux <- mean(diff(x1$dyn.aw1))
    sigx <- sd(diff(x1$dyn.aw1))
    muy <- mean(diff(x1$dyn.aw2))
    sigy <- sd(diff(x1$dyn.aw2))
    muz <- mean(diff(x1$dyn.aw3))
    sigz <- sd(diff(x1$dyn.aw3))
    
    # get cepstrum
    cept <- Mod(fft((x0$dyn.aw1 - mean(x1$dyn.aw1))/sd(x1$dyn.aw1)))
    
    #cept <- Mod(fft(x$dyn.aw1))
    cs <- rbind(cs,cept)
    
    # get time
    ts <- c(ts,x0$Time[1])
    
    # get dive
    dives0 <- c(dives0,dive)
    
    # get label 
    if(any(format(x2$Time, "%Y-%m-%d %H:%M:%S") %in% format(Data_fine0$stime, "%Y-%m-%d %H:%M:%S"))){
      labs <- c(labs,T)
      col = "blue"
    } else if(any(!x2$camon,na.rm=T)) {
      labs <- c(labs,NA)
      col = "black"
    } else {
      labs <- c(labs,F)
      col = "red"
    }
    
    # get peak jerk normalized
    diff_aw1 <- diff(x0$dyn.aw1)
    diff_aw2 <- diff(x0$dyn.aw2)
    diff_aw3 <- diff(x0$dyn.aw3)
    
    diff_aw1 <- (diff_aw1 - mux) / sigx
    diff_aw2 <- (diff_aw2 - muy) / sigy
    diff_aw3 <- (diff_aw3 - muz) / sigz
    
    jerk <- sqrt(diff_aw1^2 + diff_aw2^2 + diff_aw3^2)
    peak_jerk_norm <- c(peak_jerk_norm,max(jerk))
    
    #plot(x$dyn.aw1,
    #     ylab = "aw1",
    #     type = "l",
    #     main = paste(x$ID[1],x$divenum[1],x$Time[1]),
    #     col=col)
    #plot(x$dyn.aw1,
    #     type = "l",
    #     main = paste(x$ID[1],x$divenum[1],x$Time[1]),
    #     col=col)
    #plot(jerk,
    #     type = "l",
    #     main = paste(x$ID[1],x$divenum[1],x$Time[1]),
    #     col=col)
    #hist(diff(aw1_norm),
    #     breaks = 25,
    #     main = paste(x$ID[1],x$divenum[1],x$Time[1]),
    #     col=col)
  }
}

# plot raw acceleration at the end of the dive
for(dive in dives){
  print(dive)
  dive_df <- rawDataPos[rawDataPos$divenum %in% dive,]
  start_ascent <- max(dive_df$Time[dive_df$elev < 0.8*min(dive_df$elev)])
  dive_df <- dive_df[(dive_df$Time > start_ascent - 60) & (dive_df$Time < start_ascent),]
  dive_df$label[is.na(dive_df$label)] <- "NA"
  dive_df$abs_rajp <- abs(dive_df$rajp)
  
  dive_df0 <- dive_df %>%
    pivot_longer(cols = c("elev","jp","abs_rajp","htv"),
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


# plot raw acceleration data
ggplot(rawDataPos[rawDataPos$divenum == 3818 & rawDataPos$segnum != 0,]) +
  geom_line(aes(x = Time, y = dyn.aw1, color = label)) + 
  facet_wrap(~divenum, scales = "free")

ggplot(rawDataPos[rawDataPos$segnum != 0,]) +
  geom_line(aes(x = Time, y = dyn.aw1, color = label, group = segnum)) + 
  facet_wrap(~divenum, scales = "free")

ggplot(rawDataPos[rawDataPos$label %in% T,]) +
  geom_line(aes(x = Time, y = dyn.aw1, color = label, group = segnum)) + 
  facet_wrap(~divenum, scales = "free")

ggplot(rawDataPos[rawDataPos$segnum != 0,]) +
  geom_line(aes(x = Time, y = w_low, color = label, group = segnum)) + 
  facet_wrap(~divenum, scales = "free")

ggplot(rawDataPos[rawDataPos$segnum != 0,]) +
  geom_line(aes(x = Time, y = elev, color = label, group = segnum)) + 
  facet_wrap(~divenum, scales = "free")

# make a dataframe
cept_df <- data.frame(time = ts[-1],
                      dive = dives0,
                      cept = cs[,1:10],
                      peak_jerk_norm = peak_jerk_norm,
                      label = labs)

#plot some data
ggplot(cept_df) +
  geom_point(aes_string(x = "cept.1", y = "cept.2", color = "label"))
ggplot(cept_df) +
  geom_point(aes_string(x = "cept.1", y = "cept.3", color = "label"))
ggplot(cept_df) +
  geom_point(aes_string(x = "cept.1", y = "cept.4", color = "label"))
ggplot(cept_df) +
  geom_point(aes_string(x = "cept.2", y = "cept.3", color = "label"))
ggplot(cept_df) +
  geom_point(aes_string(x = "cept.2", y = "cept.4", color = "label"))
ggplot(cept_df) +
  geom_point(aes_string(x = "cept.3", y = "cept.4", color = "label"))

for(i in 1:10){
  print(ggplot(cept_df) +
          geom_point(aes_string(x = "time", 
                                y = paste0("cept.",i), 
                                color = "label")) + 
          facet_wrap(~dive, scales = "free"))
  
  print(ggplot(cept_df) +
          geom_point(aes_string(x = "time", 
                                y = paste0("c(0,diff(cept.",i,"))"), 
                                color = "label")) + 
          facet_wrap(~dive, scales = "free"))
}

cept_df <- cept_df %>% pivot_longer(
  cols = starts_with("cept"),
  names_to = "cept",
  names_prefix = "cept.",
  values_to = "value"
)

cept_df$cept <- as.numeric(cept_df$cept)

library(viridis)

ggplot(cept_df[cept_df$label,]) +
  geom_tile(aes_string(x = "time", 
                       y = "cept", 
                       fill = "value")) + 
  scale_fill_viridis(limits = c(0,100)) +
  labs(title = "positive segments") +
  facet_wrap(~dive, scales = "free")

ggplot(cept_df[!cept_df$label,]) +
  geom_tile(aes_string(x = "time", 
                       y = "cept", 
                       fill = "value")) + 
  scale_fill_viridis(limits = c(0,100)) +
  labs(title = "negative segments") +
  facet_wrap(~dive, scales = "free")

ggplot(rawDataPos[rawDataPos$segnum != 0,]) +
  geom_line(aes(x = Time, y = c(0,diff(w_low)), color = label, group = segnum)) + 
  facet_wrap(~divenum, scales = "free")

ggplot(rawDataPos[rawDataPos$segnum != 0,]) +
  geom_line(aes(x = Time, y = w_high, color = label, group = segnum)) + 
  facet_wrap(~divenum, scales = "free")

ggplot(rawDataPos[rawDataPos$segnum != 0,]) +
  geom_line(aes(x = Time, y = elev, color = label, group = segnum)) + 
  facet_wrap(~divenum, scales = "free")

#ggplot(rawDataPos[rawDataPos$segnum != 0,]) +
#  geom_point(aes(x = Time, y = dLogWlow, 
#                 color = label, group = segnum)) + 
#  facet_wrap(~divenum, scales = "free")

#ggplot(rawDataPos[rawDataPos$segnum != 0,]) +
#  geom_point(aes(x = Time, y = d2LogWlow, 
#                 color = label, group = segnum)) + 
#  facet_wrap(~divenum, scales = "free")

#ggplot(rawDataPos[rawDataPos$segnum != 0,]) +
#  geom_line(aes(x = Time, y = dyn.aw1, color = label)) + 
#  facet_wrap(~divenum, scales = "free")

#ggplot(rawDataPos[rawDataPos$segnum != 0,]) +
#  geom_point(aes(x = log(w_low), y = log(w_high), color = label))
