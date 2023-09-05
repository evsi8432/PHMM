library(momentuHMM)
library(dplyr)
library(tidyr)
library(mclust)
library(data.table)
library(ggplot2)
library(GGally)

#setwd("/Users/evsi8432/Documents/Research/PHMM/src/bash")

# set seed
set.seed(1)

# get command-line arguments
args <- commandArgs(trailingOnly=TRUE)

opt_file <- args[1]
args <- as.integer(args[2])

# get options
source(paste0('../opt/',opt_file))

# do the hierarchical thing
if(hier){
  source("../HHMM/plot_PHHMM.R")
  quit()
}

# define way to make titles
make_title <- function(start,end){
  title <- paste0(start,statesPerBehaviour[1])
  for(nstates in statesPerBehaviour[2:3]){
    title <- paste0(title,"-",nstates)
  }
  for(feature in features1){
    title <- paste0(title,"-",feature)
  }
  if(length(sex) > 1){
    title <- paste0(title,"_all")
  } else {
    title <- paste0(title,"_",sex)
  }
  title <- paste0(title,"_",end)
  return(title)
}

sind <- 0
if(is.na(args)){
  args_list <- sind:59
} else {
  args_list <- c(args)
}

# only load in data once
source("../HMM/load_data.R")
DataBackup <- Data

# load in the rawData
rawData <- data.frame(fread('../../../dat/Final_rawDataDown_Beth.csv'))

# add columns
rawData$Elevation <- -rawData$p
rawData$log_w_low <- log10(rawData$w_low)
rawData$log_w_high <- log10(rawData$w_high)
rawData$log_w_total <- log10(rawData$w_low + rawData$w_high)

rawDataBackup <- rawData

for(args in args_list){

# Select Model
model_ind <- (args[1] %% 5) + 1
models <- c("no","fixed_1","fixed_2","half_random","random")
model <- models[model_ind]

# select holdout whale
whale_ind <- floor(args[1] / 5) + 1
whales <- c("none","A100","A113","D21","D26","I107","I129","I145","L87","L88","R48","R58")
holdout_whale <- whales[whale_ind]

print(model)
print(holdout_whale)

### BEGIN COMPUTATION ###

# load in data
Data <- DataBackup

if(holdout_whale == "none"){
  whales = unique(Data$ID)
} else {
  whales <- holdout_whale
}
Data <- Data[Data$ID %in% whales,]
if(holdout_whale != "none"){
  Data$label <- 7
}
Data <- prepData(Data,coordNames=NULL)

# load in best hmm
files <- Sys.glob(make_title(paste0(directory,"/params/"),
                             paste0(model,"-",holdout_whale,"-*-hmm.rds")))

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

# make hmm only on heldout whale
Par0 <- getPar0(hmm)
eps <- 1e-8
if(model == "no"){
  hmm0 <- fitHMM(data=Data,
                 nbStates=N,
                 dist=hmm$conditions$dist,
                 DM=hmm$conditions$DM,
                 beta0=Par0$beta,
                 delta0=(1-eps)*(Par0$delta)+eps*rep(1/N,N),
                 Par0=Par0$Par,
                 userBounds=userBounds,
                 workBounds=workBounds,
                 nlmPar = list('print.level'=2,
                               'stepmax'=1e-100,
                               'iterlim'=1))
} else {
  hmm0 <- fitHMM(data=Data,
                 nbStates=N,
                 dist=hmm$conditions$dist,
                 fixPar=list(label = hmm$conditions$fixPar$label),
                 DM=hmm$conditions$DM,
                 beta0=Par0$beta,
                 delta0=(1-eps)*(Par0$delta)+eps*rep(1/N,N),
                 Par0=Par0$Par,
                 userBounds=userBounds,
                 workBounds=workBounds,
                 nlmPar = list('print.level'=2,
                               'stepmax'=1e-100,
                               'iterlim'=1))
}
# add pairs to hmm0
hmm0$pairs <- hmm$pairs

# find states
inverse_pairs <- order(hmm0$pairs[,2])
Data$vstates <- inverse_pairs[viterbi(hmm0)]
Data$vstates <- as.factor(Data$vstates)

# prepare colors and labels
behaviours <- c("Resting","Travelling","Foraging")
colors <- hcl(h = seq(15, 375, length = N + 1),
              l = 65, c = 100)[1:N]

group.names1 <- c()
group.names2 <- c()
group.colors1 <- c()
group.colors2 <- c()
statenum <- 1

for(i in 1:3){
  group.names1 <- c(group.names1, paste(behaviours[i],1))
  group.colors1 <- c(group.colors1, colors[statenum])
  for(j in 1:statesPerBehaviour[i]){
    group.names2 <- c(group.names2, paste(behaviours[i],j))
    group.colors2 <- c(group.colors2, colors[statenum])
    statenum <- statenum + 1
  }
}

names(group.names1) <- as.character(1:3)
names(group.colors1) <- as.character(1:3)
names(group.names2) <- as.character(1:N)
names(group.colors2) <- as.character(1:N)

labs <- c(Elevation = "Depth (meters)",
          maxDepth = "Maximum Depth (m)",
          diveDuration = "Dive Duration (s)",
          w_low = "Wiggliness (Low Frequency)",
          w_high = "Wiggliness (High Frequency)",
          log_w_low = "Wiggliness (Low Frequency) (log10)",
          log_w_high = "Wiggliness (High Frequency) (log10)",
          log_w_total = "Wiggliness (Total) (log10)",
          postDiveInt = "Post Dive Interval (s)")

if(holdout_whale == "none"){
  # plot scatterplots
  plot0 <- ggplot(Data,
                  aes(x=diveDuration,
                      y=maxDepth,
                      color=as.factor(knownState))) +
    geom_point() +
    scale_color_manual(labels = group.names1,
                       values = group.colors1) +
    stat_density_2d(color="white") +
    scale_x_log10() + scale_y_log10() +
    labs(color="Labelled Behaviour",
         y="Dive Depth (m)",
         x="Dive Duration (s)") +
    facet_wrap(~vstates,
               ncol = 2,
               labeller = as_labeller(group.names2))
  
  ggsave(make_title(paste0(directory,"/plt/"),
                    paste0(model,"_scatterplot-mddd.png")),
         plot = plot0,
         width = 8,
         height = 8,
         device='png',
         dpi=500)
  
  plot1 <- ggplot(Data,
                  aes(x=avg_w_low,
                      y=avg_w_high,
                      color=as.factor(knownState))) +
    geom_point() +
    scale_color_manual(labels = group.names1,
                       values = group.colors1) +
    stat_density_2d(color="white") +
    scale_x_log10() + scale_y_log10() +
    labs(color="Labelled Behaviour",
         x="W low",
         y="W high") +
    facet_wrap(~vstates,ncol = 2, labeller = as_labeller(group.names2))
  
  ggsave(make_title(paste0(directory,"/plt/"),
                    paste0(model,"_scatterplot-w.png")),
         plot = plot1,
         width = 8,
         height = 8,
         device='png',
         dpi=500)
  
  plot2 <- ggpairs(Data[,c(features2,"vstates")],
                   aes(colour = vstates, alpha = 0.4)) +
    scale_color_manual(labels = group.names2,
                       values = group.colors2)
  
  ggsave(make_title(paste0(directory,"/plt/"),
                    paste0(model,"_ggpairs.png")),
         plot = plot2,
         width = 8,
         height = 8,
         device='png',
         dpi=500)
}

rawData <- rawDataBackup

# add the dive types to the RawData
dive_types <- data.frame(vstate1    = Data$vstates,
                         knownState = Data$knownState,
                         divenum    = Data$divenum,
                         maxDepthCat= Data$maxDepthCat)

rawData <- left_join(rawData,dive_types,by="divenum")

cols_to_plot <- c("Elevation","log_w_total")

rawDataDownLong <- rawData %>%
  pivot_longer(cols = cols_to_plot,
               names_to = "feature")

# plot data
for (whale in whales){

  dives = Data[(Data$ID %in% whale),]$divenum
  df <- rawDataDownLong[(rawDataDownLong$divenum %in% dives) &
                          (rawDataDownLong$segnum != 0),]
  df$Time <- (df$Time - min(df$Time))/3600

  plot1 <- ggplot(df,aes(x=Time, y=value)) +
    geom_line(aes(color=factor(vstate1),
                  group=divenum)) +
    geom_hline(yintercept = 0) +
    labs(color="", y="",
         x="Elapsed time (hours)") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_color_manual(labels=group.names2,
                       values=group.colors2) +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          text = element_text(size=16)) +
    guides(colour = guide_legend(override.aes = list(linewidth=3))) +
    facet_wrap(~feature, ncol = 1, labeller = as_labeller(labs),scales = "free_y")

  ggsave(make_title(paste0(directory,"/plt/"),
                    paste0(model,"-",
                           holdout_whale,"-",
                           "profile-",whale,".png")),
         plot = plot1,
         width = 8,
         height = 8,
         device='png',
         dpi=500)

  if(model == "no"){
    depth.names <- c("shallow","med-shallow","medium","med-deep","deep")
    depth.colors <- hcl(h = seq(15, 375, length = 6),
                        l = 65, c = 100)[1:5]

    plot2 <- ggplot(df,aes(x=Time, y=value)) +
      geom_line(aes(color=factor(maxDepthCat),
                    group=divenum)) +
      geom_hline(yintercept = 0) +
      labs(color="", y="",
           x="Elapsed time (hours)") +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
      scale_color_manual(labels=depth.names,
                         values=depth.colors) +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            text = element_text(size=16)) +
      guides(colour = guide_legend(override.aes = list(linewidth=3))) +
      facet_wrap(~feature, ncol = 1, labeller = as_labeller(labs), scales = "free_y")

    ggsave(make_title(paste0(directory,"/plt/"),
                      paste0(model,"-",
                             holdout_whale,"-",
                             "profile-",whale,"-dive_types.png")),
           plot = plot2,
           width = 8,
           height = 8,
           device='png',
           dpi=500)
  }
}
}