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
library(party)
library(numDeriv)
library(data.tree)
library(scales)
library(gridExtra)

# define parameters
sex <- "Male"
date <- "2023-04-03"
holdout_whale <- "None"
nseeds <- 1
N_coarse <- 3
N_fine <- 4
N <- N_coarse * N_fine

### Load Data ###
setwd("~/Documents/Research/PHMM/src")

# load data
Data <- data.frame(fread('../../dat/Final_Data_Beth.csv'))
ethogram <- data.frame(fread('../../dat/Final_ethogram_Beth.csv'))

# get post dive int
#Data[is.infinite(Data$postDiveInt),"postDiveInt"] <- NA

# turn data positive for gamma dists
Data$postDiveInt <- exp(Data$postDiveInt) + runif(nrow(Data))/2
Data$maxDepth <- exp(Data$maxDepth)
Data$diveDuration <- exp(Data$diveDuration) + runif(nrow(Data))/2

# change post-diveinterval to a categorical
Data$postDiveIntCat <- 0
Data$postDiveIntCat <- Data$postDiveIntCat + (Data$postDiveInt < 10)
Data$postDiveIntCat <- Data$postDiveIntCat + (Data$postDiveInt < 5)
Data$postDiveIntCat <- Data$postDiveIntCat + (Data$postDiveInt < 2)
Data$postDiveIntCat <- factor(Data$postDiveIntCat)

# make wiggliness a log-scale
Data$avg_w_low <- log(Data$avg_w_low)
Data$avg_w_high <- log(Data$avg_w_high)

# turn dive Type into a factor
Data$diveType <- factor(Data$diveType,
                        levels = c("Resting", "Travelling", "Foraging", "Logging"))
Data$diveType <- addNA(Data$diveType)
Data[!is.na(Data$postDiveInt) & 
       Data$postDiveInt >= 10,"diveType"] <- "Logging"

### Label Deep vs Shallow Dives ###
shallow_thresh <- c(0,7.5)
medium_thresh <- c(10,30)
deep_thresh <- c(50,Inf)

Data$broadDiveType <- NA
Data$broadDiveType[(shallow_thresh[1] < Data$maxDepth) & (Data$maxDepth < shallow_thresh[2])] <- 1
Data$broadDiveType[(medium_thresh[1] < Data$maxDepth) & (Data$maxDepth < medium_thresh[2])] <- 2
Data$broadDiveType[(deep_thresh[1] < Data$maxDepth) & (Data$maxDepth < deep_thresh[2])] <- 3
Data$broadDiveType[Data$diveType %in% "Logging"] <- 4
Data$broadDiveType[is.na(Data$broadDiveType)] <- 5

### add label to next deep dive for foraging dives

# deep dives mean foraging
Data[Data$broadDiveType %in% 3,"diveType"] <- "Foraging"

# deep-ish dives after labels mean foraging too
for(stime in Data$stime[Data$diveType %in% "Foraging"]){
  future_dives <- Data[(Data$stime >= stime) & (Data$stime < stime + 120),]
  future_deep_dives <- future_dives[future_dives$maxDepth > 30,]
  if(nrow(future_deep_dives) >= 1){
    dive_ind <- rownames(future_deep_dives)[1]
    Data[dive_ind,"diveType"] <- "Foraging"
  }
}

### Label knownStates ###
Data$knownState <- NA
Data$knownState[Data$diveType %in% "Resting"] <- 1
Data$knownState[Data$diveType %in% "Travelling"] <- 2
Data$knownState[Data$diveType %in% "Foraging"] <- 3
Data$knownState[is.na(Data$knownState)] <- 4

span <- 10 # minutes

### change to hierarchical ###
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
    
    # check for repeated knownState
    if(forg){
      tmp2$knownState <- 3
    }
    else if (rest + trav > 1){
      tmp2$knownState <- 4
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
    time <- tail(tmp2,n=1)[1,"etime"]
  }
}
Data <- Data0

rownames(Data) <- 1:nrow(Data)

### add sex to Data ###
Data$Sex <- "Female"
Data[Data0$ID %in% c("I107","D21","L87","L88"),"Sex"]  <- "Male"
Data$Sex <- factor(Data$Sex,levels=c("Male","Female"))
Data <- prepData(Data,
                 coordNames=NULL,
                 hierLevels = c("1", "2i", "2"))

group.colors.coarse <- c("1" = "#F8766D",
                         "2" = "#7CAE00",
                         "3" = "#00BFC4")

group.colors.fine <- c("1" = "#007094",
                       "2" = "#00BE7D",
                       "3" = "#FDE333",
                       "4" = "#C77CFF")

### plot histogram of depths by sex ###
Data$maxDepth[Data$ID == "D21"] <- Data$maxDepth[Data$ID == "D21"] + 0.5

plot0 <- ggplot(Data,
                aes_string(x="maxDepth")) +
  annotate(geom = "rect", 
           xmin = 0, 
           xmax = 7.5, 
           ymin = -Inf, 
           ymax = Inf,
           fill = group.colors.fine[1], 
           alpha = 0.5) +
  annotate(geom = "rect", 
           xmin = 10, 
           xmax = 30, 
           ymin = -Inf, 
           ymax = Inf,
           fill = group.colors.fine[2], 
           alpha = 0.5) +
  annotate(geom = "rect", 
           xmin = 50, 
           xmax = Inf, 
           ymin = -Inf, 
           ymax = Inf,
           fill = group.colors.fine[3], 
           alpha = 0.5) +
  geom_histogram(position="identity",binwidth = 0.05) + 
  facet_wrap(~Sex,
             ncol=1,
             strip.position = "top") +
  labs(fill = "Dive Type", x = "Maximum dive depth (m)", y = "Number of dives") + 
  scale_x_log10(breaks = c(seq(0.5,0.9,0.1),seq(1,9,1),seq(10,90,10),seq(100,400,100)),
                labels = c(rep("",5),1,rep("",8),10,rep("",8),100,rep("",3))) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        text = element_text(size=12, family="Arial")) +
  geom_text(data = data.frame(Sex = c("Female","Male"),
                              x = c(0.5,0.5),
                              y = c(425,425),
                              label = c("Females","Males")),
            aes(x=x,y=y,label=label),
            size = 4, hjust=0, vjust=1, fontface = "bold")

print(plot0)

ggsave(paste0("../plots/",date,"/hist_maxDepth.tiff"),
       plot = plot0,
       width = 5,
       height = 3,
       units = 'in',
       device = 'tiff', 
       dpi=600)

# Do some small EDA
if (sex == "Male"){
  Data <- Data[Data$Sex %in% "Male",]
} else {
  Data <- Data[Data$Sex %in% "Female",]
}

# select the best model
best_hhmm <- NULL
max_like <- -Inf

for(rand_seed in 1:nseeds){
  hhmm <- readRDS(paste0("../params/",date,"/hierHmm_", 
                         holdout_whale,"_",
                         sex,"_",
                         rand_seed,".rds"))
  print(rand_seed)
  print(hhmm)
  if (-hhmm$mod$minimum > max_like){
    best_hhmm <- hhmm
    max_like <- -hhmm$mod$minimum
  }
}

hhmm <- best_hhmm

print(hhmm)

# define dive types
group.names <- c("1" = "Resting (Shallow)",
                 "2" = "Resting (Medium)",
                 "3" = "Resting (Deep)",
                 "4" = "Resting (Logging)",
                 "5" = "Travelling (Shallow)",
                 "6" = "Travelling (Medium)",
                 "7" = "Travelling (Deep)",
                 "8" = "Travelling (Logging)",
                 "9" = "Foraging (Shallow)",
                 "10" = "Foraging (Medium)",
                 "11" = "Foraging (Deep)",
                 "12" = "Foraging (Logging)")

group.names.coarse <- c("1" = "Resting",
                        "2" = "Travelling",
                        "3" = "Foraging")

group.names.fine <- c("1" = "Shallow",
                      "2" = "Medium",
                      "3" = "Deep",
                      "4" = "Logging")

# decode vitirbi 
vstates <- viterbi(hhmm)

vstates_fine <- rep(NA,length(vstates))
vstates_fine[vstates %in% c(1,5,9)] <- 1
vstates_fine[vstates %in% c(2,6,10)] <- 2
vstates_fine[vstates %in% c(3,7,11)] <- 3
vstates_fine[vstates %in% c(4,8,12)] <- 4

vstates_coarse <- rep(NA,length(vstates))
vstates_coarse[vstates %in% c(1,2,3,4)] <- 1
vstates_coarse[vstates %in% c(5,6,7,8)] <- 2
vstates_coarse[vstates %in% c(9,10,11,12)] <- 3

fb <- stateProbs(hhmm)
colnames(fb) <- group.names

# add dive types to unlabelled and labelled Data
Data$vstate1 <- vstates
Data$vstate1_coarse <- vstates_coarse
Data$vstate1_fine <- vstates_fine

Data <- left_join(Data,
                  data.frame(vstate1=1:N,vdivetype=group.names),
                  by="vstate1")
Data <- cbind(Data,fb)

# add the dive types to the RawData
dive_types <- data.frame(vstate1 = vstates,
                         vstate1_coarse = vstates_coarse,
                         vstate1_fine = vstates_fine,
                         divenum = Data$divenum,
                         knownState = Data$knownState)

# load in the rawData
rawData <- data.frame(fread('../../dat/Final_rawData_Beth.csv'))
rawData <- rawData[,!(names(rawData) %in% c("vstate1","label1"))]
rawData <- left_join(rawData,dive_types,by="divenum")

rawData <- rawData[,!(names(rawData) %in% c("diveType"))]
rawData <- left_join(rawData,Data[,c("divenum","diveType")],by="divenum")

rawData$Elevation <- -rawData$p

cols_to_plot <- c("Elevation")#,"w_low","w_high")

labs <- c(Elevation = "Depth (meters)", 
          head = "Heading (degrees)",
          roll = "Roll (degrees)",
          pitch = "Pitch (degrees)",
          jp = "Jerk Peak (m/s^3)",
          htv = "Heading Total Variation",
          w_low = "wiggliness (< 5 Hz)",
          w_high = "wiggliness (>= 5 Hz)",
          avg_w_low = "Average Wiggliness (< 5 Hz)",
          avg_w_high = "Average Wiggliness (>= 5 Hz)",
          Aw_1 = "x-acc (m/s^2)",
          Aw_2 = "y-acc (m/s^2)",
          Aw_3 = "z-acc (m/s^2)",
          maxDepth = "Maximum Depth (m)",
          diveDuration = "Dive Duration (s)",
          postDiveInt = "Post Dive Interval (s)",
          max_bot_jp = "log(Bottom Peak Jerk (m/s^3))",
          avg_bot_htv = "log(Average Heading Variation Rate at Bottom (rad/s))",
          avg_bot_abs_roll = "log(Average Absolute Roll at Bottom (rad))")

rawDataDownLong <- rawData[seq(1,nrow(rawData),10),] %>% 
  pivot_longer(cols = cols_to_plot, 
               names_to = "feature")

# plot the results
group.colors.coarse <- c("1" = "#F8766D",
                         "2" = "#7CAE00",
                         "3" = "#00BFC4")

group.colors.fine <- c("1" = "#007094",
                       "2" = "#00BE7D",
                       "3" = "#FDE333",
                       "4" = "#C77CFF")

### plot data
for (whale in c("D21")){#unique(Data$ID)){
  
  dives = Data[(Data$ID %in% whale),]$divenum
  df <- rawDataDownLong[(rawDataDownLong$divenum %in% dives) &
                          (rawDataDownLong$segnum != 0),]
  df$Time <- (df$Time - min(df$Time))/3600
  
  plot1 <- ggplot(df,aes(x=Time, y=value)) +
    geom_line(aes(color=factor(vstate1_fine),
                  group=divenum),
              linewidth = 0.25) + 
    geom_hline(yintercept = 0) + 
    scale_color_manual(labels=rep(NULL,4),
                       values=group.colors.fine) +
    labs(color=NULL,
         y="Dive depth (m)",
         x=NULL) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme_classic() + 
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          text = element_text(size=12,family = "Arial"),
          legend.key.height = unit(4,"mm"),
          legend.key.width = unit(4,"mm"),
          legend.position=c(.2,.45),
          legend.justification = c(1,1),
          axis.title.y = element_text(vjust = 2.5),
          axis.text.x=element_blank(),
          plot.margin=unit(c(5.5,5.5,2.25,5.5), "points")) + 
    guides(colour = guide_legend(override.aes = list(linewidth=2)))
  
  plot2 <- ggplot(df,aes(x=Time, y=value)) +
    geom_line(aes(color=factor(vstate1_coarse),
                  group=divenum),
              linewidth = 0.25) + 
    geom_hline(yintercept = 0) + 
    scale_color_manual(labels=rep(NULL,3),
                       values=group.colors.coarse) +
    labs(color=NULL,
         y="Dive depth (m)",
         x="Elapsed time (hours)") + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme_classic() + 
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          text = element_text(size=12,family = "Arial"),
          legend.key.height = unit(4,"mm"),
          legend.key.width = unit(4,"mm"),
          legend.position=c(.2,.45),
          legend.justification = c(1,1),
          axis.title.y = element_text(vjust = 2.5),
          axis.title.x = element_text(vjust = -0.4),
          plot.margin = unit(c(2.25,5.5,5.5,5.5), "points")) + 
    guides(colour = guide_legend(override.aes = list(linewidth=2)))
  
  print(plot1)
  #print(plot2)
  
  g <- arrangeGrob(plot1, plot2)
  ggsave(paste0("../plots/",date,"/profile_",whale,"_all.tiff"),
         plot = g,
         width = 5,
         height = 5,
         device='tiff', 
         dpi=600)
}

# plots of features

for(feature in c("maxDepth","diveDuration","postDiveInt")){
  
  shape1 = hhmm$mle[[feature]]['mean',1]^2 / hhmm$mle[[feature]]['sd',1]^2
  scale1 = hhmm$mle[[feature]]['sd',1]^2 / hhmm$mle[[feature]]['mean',1]
  shape2 = hhmm$mle[[feature]]['mean',2]^2 / hhmm$mle[[feature]]['sd',2]^2
  scale2 = hhmm$mle[[feature]]['sd',2]^2 / hhmm$mle[[feature]]['mean',2]
  shape3 = hhmm$mle[[feature]]['mean',3]^2 / hhmm$mle[[feature]]['sd',3]^2
  scale3 = hhmm$mle[[feature]]['sd',3]^2 / hhmm$mle[[feature]]['mean',3]
  shape4 = hhmm$mle[[feature]]['mean',4]^2 / hhmm$mle[[feature]]['sd',4]^2
  scale4 = hhmm$mle[[feature]]['sd',4]^2 / hhmm$mle[[feature]]['mean',4]
  
  plot0 <- ggplot(Data,aes_string(x=feature)) +
    scale_fill_manual(labels=group.names.fine,values=group.colors.fine) +
    geom_histogram(aes(fill=factor(vstate1_fine),y=stat(density)),alpha=0.5,position = "identity") + 
    stat_function(fun = function(x) {dgamma(x,shape=shape1,scale=scale1)}, n = 1000,
                  color = group.colors.fine[1]) +
    stat_function(fun = function(x) {dgamma(x,shape=shape2,scale=scale2)}, n = 1000,
                  color = group.colors.fine[2]) +
    stat_function(fun = function(x) {dgamma(x,shape=shape3,scale=scale3)}, n = 1000,
                  color = group.colors.fine[3]) +
    stat_function(fun = function(x) {dgamma(x,shape=shape4,scale=scale4)}, n = 1000,
                  color = group.colors.fine[4]) +
    labs(fill = "Dive Type", x = labs[feature]) + 
    ggtitle(paste("Histogram for",sex,"Killer Whale Dives from 2020 Field Season"))
  
  #print(plot0)
  
  ggsave(paste0("../plots/",date,"/hist_hier_",
                feature,"_",sex,".tiff"),
         plot = plot0,
         width = 5,
         height = 1.25,
         device='tiff', 
         dpi=500)
}

# only keep level 2 of Data
Data <- Data[Data$level %in% "2",]

### write Dive Data to csv ###
if (sex == "Male"){
  write.csv(Data[,c("ID","boutnum","TDR.dive.no","divenum","stime","etime",
                    "diveDuration","maxDepth","postDiveInt",
                    "vdivetype",group.names)],
            paste0("../params/",date,"/male_hier_dive_types.csv"))
} else {
  write.csv(Data[,c("ID","boutnum","TDR.dive.no","divenum","stime","etime",
                    "diveDuration","maxDepth","postDiveInt",
                    "vdivetype",group.names)],
            paste0("../params/",date,"/female_hier_dive_types.csv"))
}

### print out stationary distribution ###
print(solve(t(diag(N)-hhmm$CIreal$gamma$est+1),rep(1,N)))

### write summary Data to csv ###
sum <- Data %>% 
  group_by(vstate1) %>% 
  summarise(min_Depth_m = min(maxDepth,na.rm=T), max_Depth_m = max(maxDepth,na.rm=T),
            min_Duration_s = min(diveDuration,na.rm=T), max_Duration_s = max(diveDuration,na.rm=T),
            min_postDiveInt_s = min(postDiveInt,na.rm=T), max_postDiveInt_s = max(postDiveInt,na.rm=T))

sum <- left_join(data.frame(vstate1=1:N,diveType=group.names),sum)

sum[,c("mean_Depth_m","sd_Depth_m",
       "se_mean_Depth_m","se_sd_Depth_m")] <- t(rbind(hhmm$mle$maxDepth,hhmm$CIreal$maxDepth$se))

sum[,c("mean_Duration_s","sd_Duration_s",
       "se_mean_Duration_s","se_sd_Duration_s")] <- t(rbind(hhmm$mle$diveDuration,hhmm$CIreal$diveDuration$se))

sum[,c("mean_postDiveInt_s","sd_postDiveInt_s",
       "se_mean_postDiveInt_s","se_sd_postDiveInt_s")] <- t(rbind(hhmm$mle$postDiveInt,hhmm$CIreal$postDiveInt$se))

# get coarse-scale beta values
beta_vals <- hhmm$mle$beta[1,hhmm$mle$beta[1,] != -1e10]
coarse_inds <- rep(NA,N_coarse*(N_coarse-1))
for(i in 1:(N_coarse*(N_coarse-1))){
  coarse_inds[i] <- which(hhmm$mod$wpar == beta_vals[i])[1]
}

# get fine-scale beta values
beta_vals <- hhmm$mle$beta[3,hhmm$mle$beta[3,] != -1e10]
fine_inds <- matrix(rep(NA,N*(N_fine-1)),nrow=N_coarse)
k <- 1
l <- 1
for(i in 1:N_coarse){
  for(j in 1:(N_fine*(N_fine-1))){
    print(which(hhmm$mod$wpar == beta_vals[k]))
    print(l)
    if(length(which(hhmm$mod$wpar == beta_vals[k]))>1){
      fine_inds[i,j] <- which(hhmm$mod$wpar == beta_vals[k])[l]
      l <- l+1
    } else {
      fine_inds[i,j] <- which(hhmm$mod$wpar == beta_vals[k])[1]
    }
    k <- k+1
  }
}

# get diveDuration indices
dd_means <- unique(hhmm$mle$diveDuration[1,])
dd_inds <- rep(NA,N_fine)
for(i in 1:N_fine){
  dd_inds[i] <- which(hhmm$mod$wpar == log(dd_means[i]))[1]
}

# get postDiveInterval indices
pdi_means <- unique(hhmm$mle$postDiveInt[1,])
pdi_inds <- rep(NA,N_fine)
for(i in 1:N_fine){
  pdi_inds[i] <- which(exp(hhmm$mod$wpar) == pdi_means[i])[1]
}

estimate <- hhmm$mod$wpar

est_2_time <- function(estimate){
  
  # get parameters from hmm
  mean_dd <- exp(estimate[dd_inds])
  mean_pdi <- exp(estimate[pdi_inds])
  
  # get coarse-scale transition prob matrix
  Gamma_coarse <- matrix(rep(0.0,N_coarse^2),nrow = N_coarse,ncol = N_coarse)
  for (rownum in 1:N_coarse){
    row <- estimate[coarse_inds[((N_coarse-1)*(rownum-1)+1):((N_coarse-1)*(rownum))]]
    row <- append(row, 0.0, after=(rownum-1))
    Gamma_coarse[rownum,] = exp(row) / sum(exp(row))
  }
  
  # get stationary distribution
  coarse_delta <- solve(t(diag(N_coarse)-Gamma_coarse+1),rep(1,N_coarse))
  
  time <- c(coarse_delta,0)
  
  # subtract out logging from each distribution
  for(i in 1:N_coarse){
    Gamma_fine <- matrix(rep(0.0,N_fine^2), 
                         nrow = N_fine, 
                         ncol = N_fine)
    for (rownum in 1:N_fine){
      row <- estimate[fine_inds[i,((N_fine-1)*(rownum-1)+1):((N_fine-1)*(rownum))]]
      row <- c(0,row)
      Gamma_fine[rownum,] = exp(row) / sum(exp(row))
    }
    
    # get stationary distribution
    fine_delta <- solve(t(diag(N_fine)-Gamma_fine+1),rep(1,N_fine))
    
    # get time
    fine_time <- fine_delta * (mean_dd + mean_pdi) 
    fine_time <- fine_time / sum(fine_time)
    
    # get time in logging
    time[i] <- time[i] - fine_time[N_fine]
    time[N_coarse+1] <- time[N_coarse+1] + fine_time[N_fine]
  }
  return(time)
}

get_se <- function(estimate,hess,func) {
  
  # get the jacobian
  M <- jacobian(func,estimate)
  
  # get rid of indices with no curvature
  inds <- (diag(hess) != 0.0) & !(is.na(diag(hess)))
  hess <- hess[inds,inds]
  M <- M[,inds]
  
  return(sqrt(diag(M %*% solve(hess) %*% t(M))))
}

estimate <- hhmm$mod$wpar
hess <- hhmm$mod$hessian

# get percent of time
#sum$perc_time <- est_2_time(estimate)
#sum$se_perc_time <- get_se(estimate,hess,est_2_time)

# get percent of time combined
temp <- est_2_time(estimate)
temp_se <- get_se(estimate,hess,est_2_time)
sum$perc_time_combined <- rep(c(temp[1],temp[2],temp[3],temp[4]),each=3)
sum$se_perc_time_combined <- rep(c(temp_se[1],temp_se[2],temp_se[3],temp_se[4]),
                                 each=3)

# get total dives
ndives <- data.frame(table(Data$vstate1))
colnames(ndives) <- c("vstate1","num_dives")
ndives$vstate1 <- as.integer(as.character(ndives$vstate1))
sum <- left_join(sum,ndives)

if (sex == "Male"){
  write.csv(sum,paste0("../params/",date,"/male_summary_hier.csv"))
} else {
  write.csv(sum,paste0("../params/",date,"/female_summary_hier.csv"))
}

print("DONE")