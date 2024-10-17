library(momentuHMM)
library(dplyr)
library(tidyr)
library(mclust)
library(data.table)
library(ggplot2)
library(GGally)

# set seed
rand_seed <- 1
set.seed(rand_seed)

# get command-line arguments
args <- commandArgs(trailingOnly=TRUE)
args <- c("logMDDD_1-1-1_dd-30_2023-10-23.R",NA)

opt_file <- args[1]
args <- as.integer(args[2])

# get options
source(paste0('../opt/',opt_file))

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

# only load in data once
Data <- data.frame(fread("../../dat/case_study_1_data.csv"))
DataBackup <- Data

# load in the rawData
rawData <- data.frame(fread("../../dat/case_study_1_rawdata.csv"))

# add columns
rawData$Elevation <- -rawData$p
rawDataBackup <- rawData

# define whales
whales <- c("none","A100a","A100b","A113a","A113b",
            "D21a","D21b","D26a","D26b",
            "I107a","I107b","I129","I145a","I145b",
            "L87","L88","R48a","R48b","R58a","R58b")
n_whales <- length(whales)

# define models
Data <- data.frame(fread("../../dat/case_study_1_data.csv"))
Data <- Data[,c("ID","divenum","stime","logMDDD.x","logMDDD.y","knownState")]

ratio <- sum(!(Data$knownState %in% 4)) / nrow(Data)
ratio <- round(ratio,3)
models <- list()
models[[1]] <- c("fixed",  0.0)       # no weight
models[[2]] <- c("fixed",  0.5*ratio) 
models[[3]] <- c("fixed",  ratio)     # equal weight
models[[4]] <- c("fixed",  0.5 + 0.5*ratio)
models[[5]] <- c("fixed",  1.0)       # natural weight
n_models <- length(models) 

sind <- 0
args_list <- sind:(n_whales*n_models-1)

for(args in args_list){

# Set Model
model_ind <- (floor(args) %% n_models) + 1
model <- models[[model_ind]][1]
lamb  <- as.numeric(models[[model_ind]][2])

# Select Holdout Whale
whale_ind <- (floor(args / n_models) %% n_whales) + 1
holdout_whale <- whales[whale_ind]

print(model)
print(lamb)
print(holdout_whale)

### BEGIN COMPUTATION ###

# load in data
Data <- DataBackup

if(holdout_whale == "none"){
  whales_to_keep = unique(Data$ID)
} else {
  whales_to_keep <- holdout_whale
}
Data <- Data[Data$ID %in% whales_to_keep,]
if(holdout_whale != "none"){
  Data$label <- 4
}
Data <- prepData(Data,coordNames=NULL)

# load in best hmm
files <- Sys.glob(make_title(paste0(directory,"/params/"),
                             paste0(model,"-",
                                    lamb,"-",
                                    holdout_whale,"-",
                                    "*-hmm.rds")))

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
               nlmPar = list('stepmax'=1e-100,
                             'iterlim'=1))

# find states
Data$vstates <- viterbi(hmm0)
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
  group.names1 <- c(group.names1, behaviours[i])#paste(behaviours[i],1))
  group.colors1 <- c(group.colors1, colors[statenum])
  for(j in 1:statesPerBehaviour[i]){
    group.names2 <- c(group.names2, behaviours[i])#paste(behaviours[i],j))
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

rawData <- rawDataBackup

# add the dive types to the RawData
dive_types <- data.frame(vstate1    = Data$vstates,
                         knownState = Data$knownState,
                         divenum    = Data$divenum)

rawData <- left_join(rawData,dive_types,by="divenum")

cols_to_plot <- c("Elevation")

rawDataDownLong <- rawData %>%
  pivot_longer(cols = cols_to_plot,
               names_to = "feature")

# plot data
for (whale in whales_to_keep){

  dives = Data[(Data$ID %in% whale),]$divenum
  df <- rawDataDownLong[(rawDataDownLong$divenum %in% dives) &
                        (rawDataDownLong$segnum != 0),]
  stime <- as.character(min(df$Time))
  df$Time <- (df$Time - min(df$Time))/3600
  
  plot1 <- ggplot(df,aes(x=Time, y=value)) +
    geom_line(aes(color=factor(vstate1),
                  group=divenum)) +
    geom_hline(yintercept = 0) +
    labs(color="", y="",
         x=paste("Elapsed time (hours after",stime,")")) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_color_manual(labels=group.names2,
                       values=group.colors2) +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          text = element_text(size=16)) +
    guides(colour = guide_legend(override.aes = list(linewidth=3))) +
    facet_wrap(~feature, ncol = 1, 
               labeller = as_labeller(labs), scales = "free_y",
               strip.position = "left")

  plot1
  
  ggsave(paste0(directory,"/plt/",
                holdout_whale,"-",
                "profile-",whale,"-",
                 model,"-",
                 lamb,"-",
                 rand_seed,".png"),
         plot = plot1,
         width = 8,
         height = 4,
         units = "in",
         device='png',
         dpi=500)
  
  plot_known_states = TRUE
  if(plot_known_states & lamb == 0.0){
      plot2 <- ggplot(df,aes(x=Time, y=value)) +
        geom_line(aes(#color=factor(knownState),
                      group=divenum)) +
        geom_hline(yintercept = 0) +
        labs(color="", y="",
             x=paste("Elapsed time (hours after",stime,")")) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
        #scale_color_manual(labels=group.names1,
        #                   values=group.colors1) +
        theme_classic() +
        theme(strip.background = element_blank(),
              strip.placement = "outside",
              text = element_text(size=16)) +
        #guides(colour = guide_legend(override.aes = list(linewidth=3))) +
        facet_wrap(~feature, ncol = 1, 
                   labeller = as_labeller(labs), scales = "free_y",
                   strip.position = "left")
    
    ggsave(paste0(directory,"/plt/",
                  holdout_whale,"-",
                  "profile-",whale,
                  "-known_states.png"),
           plot = plot2,
           width = 8,
           height = 4,
           units = "in",
           device='png',
           dpi=500)
  }
}
}
