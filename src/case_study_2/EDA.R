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
library(mclust)
library(data.table)
library(ggplot2)
library(GGally)
library(signal)
library(oce)

setwd("/Users/evsi8432/Documents/Research/PHMM/src/case_study_2")

# set seed
set.seed(1)

# load in data
source("load_data_fine.R")
Data_fine <- Data

# load in best hmm
directory <- "../../exp/logMDDD_1-1-1_dd_30_2023-10-23/params/"
files <- Sys.glob(paste0(directory,"1-1-1-logMDDD_all_fixed-0.049-none-*-hmm.rds"))

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

# find the right dive types
sex <- c("Male","Female")
dd_thresh <- 30
md_thresh <- 0.5
md_threshs <- c(5,10,30,50)
hier = F

source("../preprocessing/load_data.R")
Data$vstate <- viterbi(hmm)

# add vstate to the fine Data
Data_fine <- left_join(x = Data_fine, y = Data[,c("divenum","vstate")])

# only keep foraging dives 
#Data_fine <- Data_fine[Data_fine$vstate %in% 3,]
#Data_fine <- Data_fine[Data_fine$bottom %in% T,]

# add labels
Data_fine$foraging_signs <- NA
#Data_fine$foraging_signs[Data_fine$echo.steady %in% T] <- "echo.steady"
Data_fine$foraging_signs[Data_fine$echo.rapid %in% T] <- "echo.rapid"
Data_fine$foraging_signs[Data_fine$chase %in% T] <- "chase"
Data_fine$foraging_signs[Data_fine$scales %in% T] <- "scales"
Data_fine$foraging_signs[Data_fine$crunch %in% T] <- "crunch"
Data_fine$foraging_signs[Data_fine$foraging %in% T] <- "forgaing"

# add features
Data_fine$delt_d <- c(NA,diff(Data_fine$ad))
Data_fine$logHtv <- log(Data_fine$htv)

Data_fine$eating <- 3
Data_fine$eating[Data_fine$crunch %in% T] <- 1
Data_fine$eating[Data_fine$scales %in% T] <- 2
Data_fine$eating[Data_fine$camon %in% F] <- NA

Data_fine$echo <- 3
Data_fine$echo[Data_fine$echo.steady %in% T] <- 1
Data_fine$echo[Data_fine$echo.rapid %in% T] <- 2
Data_fine$echo[Data_fine$camon %in% F] <- NA

Data_fine$forage <- 3
Data_fine$forage[Data_fine$chase %in% T] <- 1
Data_fine$forage[Data_fine$foraging %in% T] <- 2
Data_fine$forage[Data_fine$camon %in% F] <- NA

dives <- unique(Data_fine$divenum[!is.na(Data_fine$foraging_signs)])
dives <- unique(Data$divenum[Data$maxDepth > 50 & Data$vstate == 3])

plot <- F

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
          delt_d = "Change in Depth (m)",
          roll = "Roll (radians)",
          VeDBA = "VeDBA")

if(plot){
  
  # plot scatterplots
  plot1 <- ggplot(Data_fine,
                  aes(x=w_low,
                      y=w_high,
                      color=as.factor(foraging_signs))) +
    geom_point() +
    stat_density_2d(color="white") +
    scale_x_log10() + scale_y_log10() +
    labs(color="Labelled Behaviour") +
    facet_wrap(~foraging_signs , ncol = 3)
  
  print(plot1)
  
  plot1 <- ggplot(Data_fine,
                aes(x=htv,
                    y=delt_d,
                    color=as.factor(foraging_signs))) +
  geom_point() +
  stat_density_2d(color="white") +
  scale_x_log10() +
  labs(color="Labelled Behaviour") +
  facet_wrap(~foraging_signs , ncol = 3)
  
  print(plot1)
  
  plot1 <- ggplot(Data_fine,
                  aes(x=abs_roll,
                      y=jp,
                      color=as.factor(foraging_signs))) +
    geom_point() +
    stat_density_2d(color="white") +
    scale_y_log10() +
    labs(color="Labelled Behaviour") + 
    facet_wrap(~foraging_signs , ncol = 3)
  
  print(plot1)
  
  plot1 <- ggplot(Data_fine,
                  aes(x=roll,
                      y=pitch,
                      color=as.factor(foraging_signs))) +
    geom_point() +
    stat_density_2d(color="white") +
    labs(color="Labelled Behaviour") + 
    facet_wrap(~foraging_signs , ncol = 3)
  
  print(plot1)
  
  plot1 <- ggplot(Data_fine,
                  aes(x=rtv,
                      y=ptv,
                      color=as.factor(foraging_signs))) +
    geom_point() +
    stat_density_2d(color="white") +
    scale_x_log10() + scale_y_log10() +
    labs(color="Labelled Behaviour") +
    facet_wrap(~foraging_signs , ncol = 3)
  
  print(plot1)
  
  for(dive in dives){
    df <- Data_fine[Data_fine$divenum == dive,]
    
    df$foraging_signs[is.na(df$foraging_signs)] <- "none"
    df$Elevation <- -df$ad
    stime <- min(df$stime,na.rm=T)
    df$stime <- (df$stime - min(df$stime,na.rm=T))/60
    
    cols_to_plot <- c("Elevation","delt_d","htv","VeDBA")
    df_long <- df %>%
      pivot_longer(cols = cols_to_plot,
                   names_to = "feature")
    
    plot0 <- ggplot(df_long,aes(x=stime, y=value)) +
      geom_line(aes(color=as.factor(foraging_signs),
                    group=divenum)) +
      geom_hline(yintercept = 0) +
      labs(title = paste(df$ID[1],
                         "camon:",df$camon[3],
                         "foraging:",sum(df$forage %in% c(1,2)),
                         "eating:",sum(df$eating %in% c(1,2))),
           color="", y="",
           x="Elapsed time (minutes)") +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            text = element_text(size=16)) +
      facet_wrap(~feature, ncol = 1, labeller = as_labeller(labs), scales = "free_y")
    
    print(plot0)
    print(dive)
  }
}

# make an HMM (treat dives as independent)
Data_fine <- Data_fine[Data_fine$divenum %in% dives,]
Data_fine <- Data_fine[Data_fine$level == 2,]

# define states
bhavs <- c("descent","bottom",
           "chase","capture",
           "ascent1","ascent2")
N <- length(bhavs)

# add labels
Data_fine$knownState <- NA
Data_fine$knownState[Data_fine$chase %in% T] <- 3
Data_fine$knownState[Data_fine$foraging %in% T] <- 4

# set distributions
dist <- list()
dist[["delt_d"]] <- "norm"
#dist[["logHtv"]] <- "norm"
#dist[["roll"]] <- "vm"
dist[["logWTotal"]] <- "norm"
dist[["eating"]] <- "cat3" # crunch, scales, nothing
dist[["echo"]] <- "cat3" # steady, rapid, nothing
dist[["forage"]] <- "cat3" # chase, capture, nothing

# make the DM matrix (share features for ascending with and without a fish)
DM <- list()
for(feature in names(dist)){
  if(dist[[feature]] == "norm"){
    DM0 <- kronecker(diag(2),
                     rbind(diag((N-1)),c(rep(0,(N-2)),1)))
    DM[[feature]] <- DM0
  } else if (dist[[feature]] == "vm"){
    DM0 <- rbind(diag((N-1)),c(rep(0,(N-2)),1))
    DM[[feature]] <- DM0
  }
}

#else if (dist[[feature]] == "norm"){
#  DM[[feature]] <- list(mean = ~1, sd = ~1)
#} else if (substring(dist[[feature]], 1,3) == "cat"){
#  ncats <- as.integer(substring(dist[[feature]], 4))
#  DM[[feature]] <- list()
#  for(i in 1:(ncats-1)){
#    DM[[feature]][[paste0("prob",i)]] = ~1
#  }
#}

# Set fixed Parameters
fixPar <- list()
eps <- 1e-50
fixPar$eating <- c(eps,eps,eps, NA, eps, NA, # prob crunch each state
                   eps,eps,eps, NA, eps, NA) # prob scales each state

fixPar$echo <- c(  NA,  NA,  NA,  NA,  NA, NA,  # prob steady each state
                  eps, eps,  NA,  NA, eps, eps) # prob rapid each state

fixPar$forage <- c(eps,eps,  NA,eps,eps,eps,  # prob chase each state
                   eps,eps,eps,  NA,eps,eps) # prob capture each state

fixPar$delt_d <- c(NA, 0, 0,NA,NA, # means each state
                   NA,NA,NA,NA,NA) # sds each state

# fix beta
fixPar$beta <- c(       NA,-1e2,-1e2,  NA,-1e2,  # descent
                 -1e2,       NA,-1e2,  NA,-1e2,  # bottom
                 -1e2,  NA,       NA,  NA,-1e2,  # chase
                 -1e2,-1e2,-1e2,     -1e2, 1e2,  # catch
                 -1e2,-1e2,-1e2,-1e2,     -1e2,  # ascent 1
                 -1e2,-1e2,-1e2,-1e2,-1e2      ) # ascent 2

fixPar$delta <- c(1.0-5*eps,eps,eps,eps,eps,eps)

# set initial parameters
Par0 <- list()
Par0[["delt_d"]] <- c(c(2, 0, 0, 0,-2), # means
                  log(c(2, 1, 2, 2, 2))) # sds

Par0[["logHtv"]] <- c(c(-2,-2,-1, 0,-2), # means
                  log(c(1, 1, 2, 1, 1))) # sds

Par0[["roll"]] <- log(c(1, 1, 1, 1, 1)) # concentrations

Par0[["logWTotal"]] <- c(c(1, 1, 1, 7, 1), # means
                     log(c(2, 2, 2, 1, 2))) # sds

Par0[["eating"]] <- c(eps,eps,eps,0.1,eps,0.1, # prob crunch each state
                      eps,eps,eps,0.1,eps,0.1) # prob scales each state

Par0[["echo"]] <- c(0.1,0.1,0.1,0.1,0.1,0.1, # prob steady each state
                    eps,eps,0.1,0.1,eps,eps) # prob rapid each state

Par0[["forage"]] <- c(eps,eps,0.1,eps,eps,eps,  # prob chase each state
                      eps,eps,eps,0.1,eps,eps)  # prob capture each state

# pick initial beta
beta0  <- c(        -1,-1e2,-1e2,  -1,-1e2,  # descent
            -1e2,        -1,-1e2,  -1,-1e2,  # bottom
            -1e2,  -1,        -1,  -1,-1e2,  # chase
            -1e2,-1e2,-1e2,      -1e2, 1e2,  # catch
            -1e2,-1e2,-1e2,-1e2,      -1e2,  # ascent 1
            -1e2,-1e2,-1e2,-1e2,-1e2      ) # ascent 2

beta0 <- matrix(beta0,nrow=1)

# pick initial delta
delta0 <- matrix(c(1.0-5*eps,eps,eps,eps,eps,eps),nrow=1)

# prep data
Data_fine$whale <- Data_fine$ID
Data_fine$ID <- Data_fine$divenum
Data_fine <- prepData(Data_fine,coordNames=NULL)

# add labels using new categorical dists
source("add_labels_fine.R")

plot_raw = T

if(plot_raw){
  
  # import raw Data
  rawData <- data.frame(fread('../../../dat/Final_rawData_Beth.csv'))
  rawData <- rawData[rawData$divenum %in% dives,]
  rawData$elev <- -rawData$p
  
  # make seperate dataframe for ggplot
  
  #cols_to_plot <- c("head","pitch","roll","elev")
  #cols_to_plot <- c("Aw_1","Aw_2","Aw_3","elev")
  #cols_to_plot <- c("dyn.aw1","dyn.aw2","dyn.aw3","elev")
  cols_to_plot <- c("VeDBA","elev","d_elev")
  
  for(i in c(1,2,3)){
    
    if(i == 1){
      dives0 = neg_dives
      title0 = "neg dive"
    } else if (i == 2) {
      dives0 = pos_dives
      title0 = "pos dive"
    } else {
      dives0 = setdiff(dives, c(pos_dives,neg_dives))
      title0 = "NA dive"
    }
    
    for(dive in dives0){
      
      df <- rawData[rawData$divenum == dive,]
      df$d_elev <- c(NA,diff(stats::filter(df$elev, rep(1 / 50, 50), sides = 2)))
      
      df <- df %>% pivot_longer(cols = cols_to_plot,
                                names_to = "feature")
      
      df$Time <- df$Time - min(df$Time)
      
      plot1 <- ggplot(df,aes(x=Time, y=value)) +
        geom_line() +
      geom_hline(yintercept = 0) +
      labs(color="", y="",
           x="Elapsed time (seconds)",
           title=paste(title0,dive)) +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            text = element_text(size=16)) +
      facet_wrap(~feature, ncol = 1, scales = "free_y")
      
      print(plot1)
        
      for(axs in c()){#cols_to_plot){
        
        acc <- rawData[rawData$divenum == dive,axs]
        
        # create spectrogram
        spec = specgram(x = acc,
                        n = 256,
                        Fs = 50
        )
        
        # discard phase information
        P = abs(spec$S)
        
        # normalize
        max_P <- max(P)
        #P = P/max_P
        
        # convert to dB
        P = 10*log10(P)
        
        # config time axis
        t = spec$t
        
        max_freq = 25
        eind <- ceiling(nrow(P) * max_freq / 25)
        
        # plot spectrogram
        imagep(x = t,
               y = spec$f[1:eind],
               z = t(P[1:eind,]),
               col = oce.colorsViridis,
               ylab = 'Frequency [Hz]',
               xlab = 'Time [s]',
               zlab = paste(title0,axs,dive,"max intensity",max_P),
               zlim = c(-12,12),
               drawPalette = T,
               decimate = F
        )
      }
    }
  }
}

checkPar0(data=Data_fine,
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

hmm <- fitHMM(data=Data_fine[Data_fine$ID %in% c(pos_dives,
                                                 neg_dives),],
              nbStates=N,
              dist=dist,
              DM=DM,
              beta0=beta0,
              delta0=delta0,
              Par0=Par0,
              fixPar=fixPar,
              knownStates=knownStates[Data_fine$ID %in% c(pos_dives,
                                                          neg_dives)],
              #lambda=1.0,
              #userBounds=userBounds,
              #workBounds=workBounds,
              nlmPar = list('iterlim'=1000,
                            'print.level'=2))

directory <- "../../exp/case_study_2/"
model <- "EDA_hmm.rds"

dir.create(paste0(directory), showWarnings = FALSE)
dir.create(paste0(directory,"params"), showWarnings = FALSE)
dir.create(paste0(directory,"plt"), showWarnings = FALSE)

saveRDS(hmm,paste0(directory,"params/",model))


### plot results

Data_fine0 <- Data_fine[Data_fine$ID %in% c(pos_dives,
                                            neg_dives),]

Data_fine0$viterbi <- viterbi(hmm)

plot_dives <- function(dives,title0){
  for(dive in dives){
    df <- Data_fine0[Data_fine0$ID == dive,]
    df$Elevation <- -df$ad
    stime <- min(df$stime,na.rm=T)
    df$stime <- (df$stime - min(df$stime,na.rm=T))/60
    cols_to_plot <- c("Elevation","delt_d","logWTotal")
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
                         "camon:",df$camon[3],
                         "foraging:",sum(df$forage %in% c(1,2)),
                         "eating:",sum(df$eating %in% c(1,2)),
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
    
    print(plot0)
    print(dive)
  }
}

#plot_dives(setdiff(dives,c(pos_dives,neg_dives)),"NA")
#plot_dives(neg_dives,"neg")
plot_dives(pos_dives,"pos")
