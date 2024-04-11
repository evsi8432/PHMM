#.rs.restartR()
library(Rcpp)
library(tools)
#compileAttributes("/Users/evsi8432/Documents/Research/momentuHMM")
#package_native_routine_registration_skeleton("/Users/evsi8432/Documents/Research/momentuHMM")
#install.packages("/Users/evsi8432/Documents/Research/momentuHMM",
#                 repos=NULL,
#                 type="source")
#library(momentuHMM)
#detach("package:momentuHMM", unload=TRUE)
library(momentuHMM)

library(dplyr)
library(tidyr)
#library(mclust)
library(data.table)
library(ggplot2)
library(GGally)
library(signal)
library(oce)

directory <- "/Users/evsi8432/Documents/Research/PHMM/case_study_2"
setwd(directory)

# set options
plot <- T

# create directories
dir.create(directory, showWarnings = FALSE)
dir.create(paste0(directory,"/params"), showWarnings = FALSE)
dir.create(paste0(directory,"/plt"), showWarnings = FALSE)

# set seed
set.seed(1)

# load in data
source("src/load_data_fine.R")

# load in best hmm
directory_old <- "../exp/logMDDD_1-1-1_dd_30_2023-10-23/params/"
files <- Sys.glob(paste0(directory_old,"1-1-1-logMDDD_all_fixed-0.049-none-*-hmm.rds"))

best_hmm <- NULL
best_nll <- Inf
for(file in files){
  hmm <- readRDS(file)
  if(hmm$mod$minimum < best_nll){
    best_hmm <- hmm
    best_nll <- hmm$mod$minimum
  }
}
hmm <- best_hmm

# find the vitirbi dive types from the previous case study
sex <- c("Male","Female")
dd_thresh <- 30
md_thresh <- 0.5
md_threshs <- c(5,10,30,50)
hier = F
source("src/load_data_coarse.R")
Data$vstate <- viterbi(hmm)
dive_types <- Data[,c("divenum","vstate")]

# attach dive types to current dive data
sex <- c("Male","Female")
dd_thresh <- 2
md_thresh <- 0.5
md_threshs <- c(5,10,30,50)
hier = F
source("src/load_data_coarse.R")
Data <- left_join(x = Data, y = dive_types)
Data$vstate <- nafill(Data$vstate, "locf")

# add vstate to the fine Data
Data_fine <- left_join(x = Data_fine, y = Data[,c("divenum","vstate")])

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

# add features
Data_fine$delt_d <- c(NA,diff(Data_fine$ad))
Data_fine$logHtv <- log(Data_fine$htv)
Data_fine$logJpNorm <- log(Data_fine$jp_normed)

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
  
  # label dive
  pos_dive <- any(dive_df$crunch & dive_df$etime > first_crunch_time, na.rm=T) & (dive != 3957)
  crunch_dive <- pos_dive & any(dive_df$crunch & dive_df$etime < last_crunch_time, na.rm=T)
  neg_dive <- all(dive_df$foraging_signs %in% "none")
  
  if (crunch_dive){
    pos_dives <- c(pos_dives,dive)
    crunch_ind <- (Data_fine$divenum %in% dive) & (Data_fine$crunch %in% T) & (Data_fine$etime > first_crunch_time)
    crunch_ind <- which(crunch_ind)[1]
    Data_fine$knownState[crunch_ind] <- 4
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

if(plot){
  
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
  for(dive in 403){#pos_dives){
    print(dive)
    dive_df <- rawData[rawData$divenum %in% dive,]
    start_ascent <- max(dive_df$Time[dive_df$elev < 0.7*min(dive_df$elev)])
    dive_df <- dive_df[(dive_df$Time > start_ascent - 60) & (dive_df$Time < start_ascent+60),]
    dive_df$label[is.na(dive_df$label)] <- "NA"
    dive_df$abs_rajp <- abs(dive_df$rajp)
    
    dive_df0 <- dive_df %>%
      pivot_longer(cols = c("elev","jp_normed","abs_rajp","htv"),
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
}

# define states
bhavs <- c("descent","bottom",
           "chase","crunch",
           "ascent1","ascent2")
N <- length(bhavs)

# set distributions
dist <- list()
dist[["knownState"]] <- paste0("cat",N+1)
dist[["delt_d"]] <- "norm"
#dist[["logHtv"]] <- "norm"
#dist[["logJpNorm"]] <- "norm"
#dist[["rajp"]] <- "vm"
dist[["logWLow"]] <- "norm"
#dist[["eating"]] <- "cat3" # crunch, scales, nothing
#dist[["echo"]] <- "cat3" # steady, rapid, nothing
#dist[["forage"]] <- "cat3" # chase, capture, nothing

# make the DM matrix (share features for ascending with and without a fish)
DM <- list()
for(feature in names(dist)){
  if(dist[[feature]] == "norm"){
    DM0 <- kronecker(diag(2),
                     rbind(diag((N-1)),c(rep(0,(N-2)),1)))
  } else if (dist[[feature]] == "vm"){
    DM0 <- rbind(diag((N-1)),c(rep(0,(N-2)),1))
  } else if (feature == "knownState") {
    DM0 <- diag(N^2)
  }
  DM[[feature]] <- DM0
}

# Set fixed Parameters
fixPar <- list()

eps <- 1e-50

fixPar$delt_d <- c(NA, 0, 0,NA,NA, # means each state
                   NA,NA,NA,NA,NA) # sds each state

fixPar$knownState <- c(  0,-99,-99,-99,-99,-99, # prob desc label, each state
                       -99,  0,-99,-99,-99,-99, # prob bot label, each state
                       -99,-99,  0,-99,-99,-99, # prob chase label, each state
                       -99,-99,-99,  0,-99,-99, # prob capture label, each state
                       -99,-99,-99,-99,  0,-99, # prob asc 1 label, each state
                       -99,-99,-99,-99,-99,  0) # prob asc 2 label, each state

# fix beta
fixPar$beta <- c(       NA,  NA,  NA,  NA,-1e2,  # descent
                 -1e2,       NA,  NA,  NA,-1e2,  # bottom
                 -1e2,  NA,       NA,  NA,-1e2,  # chase
                 -1e2,-1e2,-1e2,     -1e2, 1e2,  # crunch
                 -1e2,-1e2,-1e2,-1e2,     -1e2,  # ascent 1
                 -1e2,-1e2,-1e2,-1e2,-1e2      ) # ascent 2

fixPar$delta <- c(1.0-5*eps,eps,eps,eps,eps,eps)

# set initial parameters
Par0 <- list()
Par0[["delt_d"]] <- c(c(2, 0, 0, 0,-2), # means
                  log(c(2, 0.25, 0.5, 2, 2))) # sds

Par0[["logWLow"]] <- c(c(2, 3, 5, 7, 2), # means
                   log(c(2, 2, 3, 1, 2))) # sds

Par0[["logHtv"]] <- c(c(-2,-2,-1, 0  ,-2), # means
                  log(c( 1, 1, 2, 0.5, 1))) # sds

Par0[["rajp"]] <- log(c(2, 2, 1, 0.1, 2)) # concentrations

Par0[["logJpNorm"]] <- c(c(1, 1, 2, 4, 1), # means
                     log(c(1, 1, 1, 1, 1))) # sds

Par0[["knownState"]] <- c(  0,-99,-99,-99,-99,-99, # prob desc label, each state
                          -99,  0,-99,-99,-99,-99, # prob bot label, each state
                          -99,-99,  0,-99,-99,-99, # prob chase label, each state
                          -99,-99,-99,  0,-99,-99, # prob capture label, each state
                          -99,-99,-99,-99,  0,-99, # prob asc 1 label, each state
                          -99,-99,-99,-99,-99,  0) # prob asc 2 label, each state

# pick initial beta
beta0  <- c(       -5,-1e2,-1e2,  -1,-1e2,  # descent
            -1e2,        0,-1e2,  -3,-1e2,  # bottom
            -1e2,   0,       -3,  -3,-1e2,  # chase
            -1e2,-1e2,-1e2,     -1e2, 1e2,  # crunch
            -1e2,-1e2,-1e2,-1e2,      -1e2, # ascent 1
            -1e2,-1e2,-1e2,-1e2,-1e2      ) # ascent 2

beta0 <- matrix(beta0,nrow=1)

# pick initial delta
delta0 <- matrix(c(1.0-5*eps,eps,eps,eps,eps,eps),nrow=1)

# prep data
Data_fine$whale <- Data_fine$ID
Data_fine$ID <- Data_fine$divenum
Data_fine_final <- prepData(Data_fine[Data_fine$divenum %in% dives,
                                      c("ID","divenum","stime","ad","knownState",
                                        "delt_d","logHtv","rajp","logJpNorm","logWLow")],
                            coordNames=NULL)

checkPar0(data=Data_fine_final,
          nbStates=N,
          dist=dist,
          DM=DM,
          beta0=beta0,
          delta0=delta0,
          Par0=Par0,
          fixPar=fixPar,
          #userBounds=userBounds,
          #workBounds=workBounds
          )

# 5-fold cross validation (T = foraging)
test_dives <- c()
#test_dives <- c(4585, 5975, 3059, 3177, 3169)
#test_labs <- c(T,T,F,F,F)

#test_dives <- c(3818, 5503, 3231, 5566)
#test_labs <- c(T,T,F,F)

#test_dives <- c(5553, 3161, 5579, 5489)
#test_labs <- c(T,T,F,F)

#test_dives <- c(4264, 5541, 3157, 5585)
#test_labs <- c(T,F,F,F)

#test_dives <- c(3932, 3260, 5572, 3194)
#test_labs <- c(T,F,F,F)

#train_dives <- setdiff(c(pos_dives,neg_dives),test_dives)
train_dives <- setdiff(dives,test_dives)
#lambda <- 1
lambda <- 1e-4

# get knownStates for lambda purposes
knownStates <- Data_fine_final$knownState[Data_fine_final$ID %in% train_dives]
knownStates[!(knownStates %in% 4)] <- NA

# add category for knownStates Data_fine_final
Data_fine_final$knownState[is.na(Data_fine_final$knownState)] <- N+1

# fit the HMM
hmm <- fitHMM(data=Data_fine_final[Data_fine_final$ID %in% train_dives,],
              nbStates=N,
              dist=dist,
              DM=DM,
              beta0=beta0,
              delta0=delta0,
              Par0=Par0,
              fixPar=fixPar,
              knownStates=knownStates,
              lambda=lambda,
              stateNames = bhavs,
              #userBounds=userBounds,
              #workBounds=workBounds,
              nlmPar = list('iterlim'=1000,
                            'print.level'=2))

plot(hmm)

model <- paste0("EDA_hmm_",
                str(lambda),"_",
                str(test_dives[1]),".rds")
saveRDS(hmm,paste0(directory,"/params/",model))

### plot results ###

Data_fine0 <- Data_fine_final[Data_fine_final$ID %in% train_dives,]
Data_fine0$viterbi <- viterbi(hmm)

labs <- c(Elevation = "Depth (meters)",
          maxDepth = "Maximum Depth (m)",
          diveDuration = "Dive Duration (s)",
          w_low = "Wiggliness (Low Frequency)",
          w_high = "Wiggliness (High Frequency)",
          logWLow = "Wiggliness (Low Frequency) (log10)",
          logWHigh = "Wiggliness (High Frequency) (log10)",
          logWTotal = "Wiggliness (Total) (log10)",
          postDiveInt = "Post Dive Interval (s)",
          htv = "Heading Total Variation (m)",
          logHtv = "Heading Total Variation (m) (log10)",
          logJpNorm = "Jerk Peak (normalized, log10)",
          delt_d = "Change in Depth (m)",
          roll = "Roll (rad)",
          rajp = "Roll @ Peak Jerk (rad)",
          VeDBA = "VeDBA")

plot_dives <- function(dives,title0){
  for(dive in dives){
    df <- Data_fine0[Data_fine0$ID == dive,]
    df$Elevation <- -df$ad
    stime <- min(df$stime,na.rm=T)
    df$stime <- (df$stime - min(df$stime,na.rm=T))/60
    cols_to_plot <- c("Elevation","delt_d","logWLow")#,"logJpNorm","rajp","logHtv")
    df_long <- df %>%
      pivot_longer(cols = cols_to_plot,
                   names_to = "feature")
    
    
    colors <- hcl(h = seq(15, 375, length = 6 + 1),
                  l = 65, c = 100)[1:6]
    colors <- c("1" = colors[1],
                "2" = colors[2],
                "3" = colors[3],
                "4" = colors[4],
                "5" = colors[5],
                "6" = colors[6])
    
    plot0 <- ggplot(df_long,aes(x=stime, y=value)) +
      geom_line() +
      geom_point(aes(color=factor(viterbi,levels=1:6),
                     group=divenum)) +
      geom_hline(yintercept = 0) +
      labs(title = paste(df$ID[1],
                         title0),
           color="", y="",
           x="Elapsed time (minutes)") +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
      scale_color_manual(labels = c("1" = "descent",
                                    "2" = "bottom",
                                    "3" = "chase",
                                    "4" = "capture",
                                    "5" = "ascent w/o fish",
                                    "6" = "ascent w/ fish"),
                           values = colors) +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            text = element_text(size=16)) +
      facet_wrap(~feature, ncol = 1, labeller = as_labeller(labs), scales = "free_y")
    
    ggsave(paste0(directory,"/plt/",dive,".png"), 
           plot0, 
           width = 8, height = 6)
    print(plot0)
    print(dive)
  }
}

#plot_dives(setdiff(dives,c(pos_dives,neg_dives)),"NA")
plot_dives(neg_dives,"neg")
plot_dives(pos_dives,"pos")
plot_dives(na_dives,"NA")

plot(hmm)
