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
library(data.table)
library(ggplot2)
library(GGally)
library(signal)
library(oce)
library(pROC)
library(latex2exp)
library(zoo)

directory <- "/Users/evsi8432/Documents/Research/PHMM/case_study_2"
setwd(directory)

# set options

plot <- T # whether to plot results
load_raw <- F # whether to load raw Data from scratch

args = commandArgs(trailingOnly=TRUE)

K <- 1 # number of cross-validations (one means just do all the data)
lambda <- 0.01 # lambda for paper
num_seeds <- 10 # number of random seeds

# set distributions
dist <- list()
dist[["delt_d"]] <- "norm"
#dist[["rajp"]] <- "vm"
dist[["htv"]] <- "gamma"
dist[["jp_normed"]] <- "gamma"

print("fitting PHMM...")
best_hmm <- NULL
max_ll <- -Inf

for(rand_seed in 1:num_seeds){
  
  directory <- "/Users/evsi8432/Documents/Research/PHMM/case_study_2"
  
  model_name <- paste0("hmm_",
                       paste0(setdiff(names(dist),"knownState"),collapse = "_"),"_",
                       log10(lambda),"_",
                       1,"_",1,"_", # fold 1 of 1
                       rand_seed,".rds")
  
  hmm <- readRDS(paste0(directory,"/params/",model_name))
  
  print(-hmm$mod$minimum)
  
  if(-hmm$mod$minimum > max_ll){
    best_hmm <- hmm
    max_ll <- -hmm$mod$minimum
    print("new best hmm")
  }
}
hmm <- best_hmm

# get threshold
Par0 <- getPar0(hmm)
eps <- 1e-10
N <- 6

hmm0 <- fitHMM(data=hmm$data,
               nbStates=N,
               dist=hmm$conditions$dist,
               DM=hmm$conditions$DM,
               beta0=Par0$beta,
               delta0=(1-eps)*(Par0$delta)+eps*rep(1/N,N),
               Par0=Par0$Par,
               stateNames = hmm$stateNames,
               nlmPar = list('stepmax'=1e-100,
                             'iterlim'=1))

hmm0$data$viterbi <- viterbi(hmm0)
hmm0$data$p_catch <- stateProbs(hmm0)[,4]
hmm0$data$p_no_fish <- stateProbs(hmm0)[,5]
hmm0$data$p_fish <- stateProbs(hmm0)[,6]

probs <- c()
labs <- c()
divenums <- c()

for(divenum in levels(hmm0$data$ID)){
  
  divenums <- c(divenums,divenum)
  
  if(divenum %in% pos_dives){
    labs <- c(labs,TRUE)
  } else if (divenum %in% neg_dives){
    labs <- c(labs,FALSE)
  } else {
    labs <- c(labs,NA)
  }
  
  p_fish <- tail(hmm0$data$p_fish[hmm0$data$ID == divenum],1)
  p_fish <- p_fish + tail(hmm0$data$p_catch[hmm0$data$ID == divenum],1)
  print(p_fish)
  
  probs <- c(probs,p_fish)
}

thresh <- 0.5#min(probs[labs %in% T])
pred_pos_dives <- divenums[probs > thresh]

# get the killer whale ID for each positive dive
pred_pos_ids <- Data$ID[Data$divenum %in% pred_pos_dives]

# plot the predicted positive dives
### plot results ###

labs <- c(Elevation = "Depth (m)",
          maxDepth = "Maximum Depth (m)",
          diveDuration = "Dive Duration (s)",
          w_low = "Wiggliness (Low Frequency)",
          w_high = "Wiggliness (High Frequency)",
          logWLow = "Wiggliness (Low Frequency) (log10)",
          logWHigh = "Wiggliness (High Frequency) (log10)",
          logWTotal = "Wiggliness (Total) (log10)",
          postDiveInt = "Post Dive Interval (s)",
          htv = "Heading Total Variation (rad / s)",
          logHtv = "Heading Total Variation (rad) (log10)",
          logJpNorm = "Jerk Peak (normalized, log10)",
          jp_normed = "Jerk Peak (normalized)",
          delt_d = "Change in Depth (m)",
          roll = "Roll (rad)",
          rajp = "Roll @ Peak Jerk (rad)",
          p_catch = "Probability of Catch",
          VeDBA = "VeDBA")

plot_dives <- function(dives,df){
  
  for(dive in dives){
    
    if(dive %in% pos_dives){
      title0 <- "pos"
    } else if (dive %in% neg_dives){
      title0 <- "neg"
    } else {
      title0 <- "NA"
    }
    
    title1 <- tail(df$p_fish[df$ID == dive],1)
    title1 <- title1 + tail(df$p_catch[df$ID == dive],1)
    title1 <- round(title1, 3)
    
    dive_df <- df[df$ID == dive,]
    dive_df$Elevation <- -dive_df$ad
    stime <- min(dive_df$stime,na.rm=T)
    dive_df$stime <- (dive_df$stime - min(dive_df$stime,na.rm=T))/60
    cols_to_plot <- c("Elevation","p_catch",setdiff(names(dist),c("knownState")))
    dive_df_long <- dive_df %>%
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
    
    if(lambda == 1.0){
      legend.pos <- "right"
      wid <- 8
    } else {
      legend.pos <- "none"
      wid <- 6
    }
    
    plot0 <- ggplot(dive_df_long,aes(x=stime, y=value)) +
      geom_line() +
      geom_point(aes(color=factor(viterbi,levels=1:6),
                     group=divenum,
                     shape=viterbi==4),
                 size = 2) +
      geom_hline(yintercept = 0) +
      labs(title = TeX(paste("$\\alpha =", lambda, "$, ",
                             "$P(X_{s,T_s} \\in \\{4,6\\} \\ | \\ Y_s) =", title1, "$",
                             Data$ID[Data$divenum %in% dive])),
           color="", y="",
           x="Time (min)") +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
      scale_color_manual(limits = c("1","2","3","4","5","6"),
                         labels = c("1" = "descent",
                                    "2" = "bottom",
                                    "3" = "chase",
                                    "4" = "capture",
                                    "5" = "ascent w/o fish",
                                    "6" = "ascent w/ fish"),
                         values = colors) +
      scale_shape_manual(values = c(19,8),
                         guide = "none") +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            text = element_text(size=16),
            legend.position = legend.pos,
            plot.title = element_text(hjust = 0.5)) +
      facet_wrap(~feature, ncol = 1, labeller = as_labeller(labs), scales = "free_y")
    
    ggsave(paste0(directory,"/plt/profile_",
                  paste0(setdiff(names(dist),"knownState"),collapse = "_"),"_",
                  round(log10(lambda),3),"_",
                  K,"_",
                  title0,"_",
                  dive,"_",
                  Data$ID[Data$divenum %in% dive],".png"), 
           plot0, 
           width = wid, height = 8)
    print(plot0)
    print(dive)
  }
}

plot_dives(pred_pos_dives,hmm0$data)

# find number of catches for the southerns
n_SRKW_catches <- sum(pred_pos_ids %in% c("L87","L88"))
n_NRKW_catches <- sum(!(pred_pos_ids %in% c("L87","L88")))

# load in old HMM for number of foraging dives
lambda <- 0.049
num_seeds <- 10

for(rand_seed in 1:num_seeds){
  
  directory <- "/Users/evsi8432/Documents/Research/PHMM/exp/logMDDD_1-1-1_dd_30_2023-10-23"
  
  model_name <- paste0("1-1-1-logMDDD_all_fixed-",
                       lambda,
                       "-none-",
                       rand_seed,
                       "-hmm.rds")
  
  hmm_coarse <- readRDS(paste0(directory,"/params/",model_name))
  
  print(-hmm_coarse$mod$minimum)
  
  if(-hmm_coarse$mod$minimum > max_ll){
    best_hmm_coarse <- hmm_coarse
    max_ll <- -hmm_coarse$mod$minimum
    print("new best hmm")
  }
}
hmm_coarse <- best_hmm_coarse
hmm_coarse$data$p_forg <- stateProbs(hmm_coarse)[,3]
Data <- left_join(Data,hmm_coarse$data[,c("divenum","p_forg")])
Data$p_forg <- na.locf(Data$p_forg)

# find total time for each population
whale_time <- function(whale){ 
  stime <- head(Data$stime[Data$ID %in% whale],1)
  etime <- tail(Data$stime[Data$ID %in% whale],1)
  return(etime - stime)
}
SRKW_time <- whale_time("L87") + whale_time("L88")
NRKW_time <- whale_time("A100a") + whale_time("A100b") 
NRKW_time <- NRKW_time + whale_time("A113a") + whale_time("A113b") 
NRKW_time <- NRKW_time + whale_time("D21a") + whale_time("D21b") 
NRKW_time <- NRKW_time + whale_time("D26a") + whale_time("D26b") 
NRKW_time <- NRKW_time + whale_time("I107a") + whale_time("I107b") 
NRKW_time <- NRKW_time + whale_time("I129")
NRKW_time <- NRKW_time + whale_time("I145a") + whale_time("I145b") 
NRKW_time <- NRKW_time + whale_time("R48a") + whale_time("R48b") 
NRKW_time <- NRKW_time + whale_time("R58a") + whale_time("R58b") 

# find foraging time for each population
whale_foraging_time <- function(whale){ 
  foraging_time <- 0
  whale_data <- Data[Data$ID %in% whale,]
  for(i in 1:(nrow(whale_data)-1)){
    if(whale_data$p_forg[i] > 0.5){
      foraging_time <- foraging_time + (whale_data$stime[i+1] - whale_data$stime[i])
    }
  }
  if(whale_data$p_forg[nrow(whale_data)] > 0.5){
    foraging_time <- foraging_time + (whale_data$etime[nrow(whale_data)] - whale_data$stime[nrow(whale_data)])
  }
  return(foraging_time)
}

SRKW_forg_time <- whale_foraging_time("L87") + whale_foraging_time("L88")
NRKW_forg_time <- whale_foraging_time("A100a") + whale_foraging_time("A100b") 
NRKW_forg_time <- NRKW_forg_time + whale_foraging_time("A113a") + whale_foraging_time("A113b") 
NRKW_forg_time <- NRKW_forg_time + whale_foraging_time("D21a") + whale_foraging_time("D21b") 
NRKW_forg_time <- NRKW_forg_time + whale_foraging_time("D26a") + whale_foraging_time("D26b") 
NRKW_forg_time <- NRKW_forg_time + whale_foraging_time("I107a") + whale_foraging_time("I107b") 
NRKW_forg_time <- NRKW_forg_time + whale_foraging_time("I129")
NRKW_forg_time <- NRKW_forg_time + whale_foraging_time("I145a") + whale_foraging_time("I145b") 
NRKW_forg_time <- NRKW_forg_time + whale_foraging_time("R48a") + whale_foraging_time("R48b") 
NRKW_forg_time <- NRKW_forg_time + whale_foraging_time("R58a") + whale_foraging_time("R58b") 


# get percentage of time foraging

perc_SRKW <- as.numeric(SRKW_forg_time, units = "secs") / as.numeric(SRKW_time, units = "secs") 
perc_NRKW <- as.numeric(NRKW_forg_time, units = "secs") / as.numeric(NRKW_time, units = "secs") 
print(perc_SRKW)
print(perc_NRKW)

catch_per_hr_SRKW <- n_SRKW_catches / as.numeric(SRKW_time, units = "hours")
catch_per_hr_NRKW <- n_NRKW_catches / as.numeric(NRKW_time, units = "hours")
print(catch_per_hr_SRKW)
print(catch_per_hr_NRKW)

catch_per_hr_forg_SRKW <- n_SRKW_catches / as.numeric(SRKW_forg_time, units = "hours")
catch_per_hr_forg_NRKW <- n_NRKW_catches / as.numeric(NRKW_forg_time, units = "hours")
print(catch_per_hr_forg_SRKW)
print(catch_per_hr_forg_NRKW)
