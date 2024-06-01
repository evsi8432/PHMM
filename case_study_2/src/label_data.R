# add echolocation ethogram
dates0 <- data.frame(Animal = c("I145","D26","I107"), 
                     date0 = c("2020-08-30","2020-08-31","2020-08-25"))

ac_ethogram <- read.csv("../../dat/ethograms/echolocation_ethogram.csv")

ac_ethogram <- left_join(x = ac_ethogram, y = dates0)
ac_ethogram$stime <- as.POSIXct(paste(ac_ethogram$date0,ac_ethogram$time),
                                format="%Y-%m-%d %H:%M:%OS",
                                tz = "UTC")
ac_ethogram$etime <- ac_ethogram$stime + ac_ethogram$duration.sec

# go through acoustic ethogram
for (i in seq_along(ac_ethogram$stime)){
  
  # do clicks
  if (ac_ethogram$broad.type[i] %in% c('clicks','click train')){
    inds <- which((Data_fine$stime <= ac_ethogram$etime[i]) & 
                    (Data_fine$etime >= ac_ethogram$stime[i]))
    Data_fine$echo.steady[inds] <- TRUE
  }
  
  # do buzzes
  if (ac_ethogram$broad.type[i] %in% c('buzz')){
    inds <- which((Data_fine$stime <= ac_ethogram$etime[i]) & 
                    (Data_fine$etime >= ac_ethogram$stime[i]))
    Data_fine$echo.rapid[inds] <- TRUE
  }
  
  # do crunches
  if (ac_ethogram$broad.type[i] %in% c('crunch')){
    inds <- which((Data_fine$stime <= ac_ethogram$etime[i]) & 
                    (Data_fine$etime >= ac_ethogram$stime[i]))
    Data_fine$crunch[inds] <- TRUE
  }
}

# add labels
Data_fine$foraging_signs <- NA
#Data_fine$foraging_signs[Data_fine$echo.steady %in% T] <- "echo.steady"
#Data_fine$foraging_signs[Data_fine$echo.rapid %in% T] <- "echo.rapid"
#Data_fine$foraging_signs[Data_fine$chase %in% T] <- "chase"
Data_fine$foraging_signs[Data_fine$scales %in% T] <- "scales"
Data_fine$foraging_signs[Data_fine$crunch %in% T] <- "crunch"
Data_fine$foraging_signs[Data_fine$foraging %in% T] <- "foraging"
Data_fine$foraging_signs[is.na(Data_fine$foraging_signs) & Data_fine$camon == 1] <- "none"

# label positive dives and crunches based on crunches at the end
dives <- Data$divenum[Data$maxDepth > 30 & !is.na(Data$divenum)]
pos_dives <- c()
neg_dives <- c()
na_dives <- c()

Data_fine$knownState <- NA

for(dive in dives){
  
  # get dive df
  dive_df <- Data_fine[Data_fine$divenum %in% dive,]
  first_crunch_time <- max(dive_df$etime[dive_df$ad > 0.7*max(dive_df$ad)]) - 30
  last_crunch_time <- max(dive_df$etime[dive_df$ad > 0.7*max(dive_df$ad)])
  
  # label first observation as descent
  Data_fine$knownState[(Data_fine$divenum %in% dive) & (Data_fine$stime == dive_df$stime[1])] <- 1
  
  # label dive
  pos_dive <- any((dive_df$crunch | dive_df$foraging) & dive_df$etime > first_crunch_time, na.rm=T) & !(dive %in% c(3161,3957))
  crunch_dive <- pos_dive & any(dive_df$crunch & dive_df$etime < last_crunch_time, na.rm=T)
  neg_dive <- all(dive_df$foraging_signs %in% "none")
  
  if (crunch_dive){
    pos_dives <- c(pos_dives,dive)
    crunch_ind <- (Data_fine$divenum %in% dive) & (Data_fine$crunch %in% T) & (Data_fine$etime > first_crunch_time)
    crunch_ind <- which(crunch_ind)[1]
    
    Data_fine$knownState[crunch_ind] <- 4
    Data_fine$knownState[(Data_fine$divenum %in% dive) & (Data_fine$etime == tail(dive_df$etime,1))] <- 6
  } else if (pos_dive){
    pos_dives <- c(pos_dives,dive)
    Data_fine$knownState[(Data_fine$divenum %in% dive) & (Data_fine$etime == tail(dive_df$etime,1))] <- 6
  } else if (neg_dive){
    neg_dives <- c(neg_dives,dive)
    Data_fine$knownState[(Data_fine$divenum %in% dive) & (Data_fine$etime == tail(dive_df$etime,1))] <- 5
  } else {
    na_dives <- c(na_dives,dive)
  }
}

for(dive in c(3957,pos_dives)){
  
  dive_df <- Data_fine[Data_fine$divenum %in% dive,]
  if(T){
    print(dive)
    dive_df$elev <- -dive_df$ad
    plot0 <- ggplot(dive_df, aes(x = stime, y = -ad)) +
      geom_line() + 
      labs(title = paste0("Whale: ",dive_df$ID[1],
                          "    Date: ",substr(dive_df$stime[1],1,10),
                          "    TDR dive number: ",dive_df$TDR.dive.no[1])) + 
      geom_hline(yintercept = 0) +
      geom_vline(aes(xintercept = stime),
                 data = dive_df %>% dplyr::filter(crunch %in% T),
                 color = "red") + 
      geom_vline(aes(xintercept = stime),
                 data = dive_df %>% dplyr::filter(echo.rapid %in% T),
                 color = "black")
    
    print(plot0)
    
    ggsave(paste0(dive_df$ID[1],"_",
                  dive_df$TDR.dive.no[1],
                  ".png"), 
           plot0, 
           width = 4, height = 4)
  }
}
