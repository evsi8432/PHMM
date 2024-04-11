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
library(mlbench)
library(caret)
library(rattle)

directory <- "/Users/evsi8432/Documents/Research/PHMM/case_study_2"
setwd(directory)

# create directories
dir.create(directory, showWarnings = FALSE)
dir.create(paste0(directory,"/params"), showWarnings = FALSE)
dir.create(paste0(directory,"/plt"), showWarnings = FALSE)

# set seed
set.seed(1)

# load in data
source("src/load_data_fine.R")
Data_fine <- Data

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

# only keep foraging dives 
Data_fine <- Data_fine[Data_fine$vstate %in% 3,]

# add dive duration in Data_fine
Data_fine <- Data_fine[,names(Data_fine) != "diveDuration"]
Data_fine <- left_join(x = Data_fine, y = Data[,c("divenum","diveDuration")])

# add foraging signs
Data_fine$foraging_signs <- NA
Data_fine$foraging_signs[Data_fine$echo.steady %in% T] <- "echo.steady"
Data_fine$foraging_signs[Data_fine$echo.rapid %in% T] <- "echo.rapid"
Data_fine$foraging_signs[Data_fine$chase %in% T] <- "chase"
Data_fine$foraging_signs[Data_fine$scales %in% T] <- "scales"
Data_fine$foraging_signs[Data_fine$crunch %in% T] <- "crunch"
Data_fine$foraging_signs[Data_fine$foraging %in% T] <- "forgaing"

# add features
Data_fine$delt_d <- c(NA,diff(Data_fine$ad))
Data_fine$logHtv <- log(Data_fine$htv)
Data_fine$logRtv <- log(Data_fine$rtv)
Data_fine$logPtv <- log(Data_fine$ptv)

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

# add confirmed biting to Data_fine. Note that F means "unknown"
Data_fine$label <- F

lab_df <- read.csv("dat/captures.csv")

lab_df$time[lab_df$Whale == "A113"] <- paste("2020-08-22",lab_df$time[lab_df$Whale == "A113"],"UTC")
lab_df$time[lab_df$Whale == "D26"] <- paste("2020-08-31",lab_df$time[lab_df$Whale == "D26"],"UTC")
lab_df$time[lab_df$Whale == "I107"] <- paste("2020-08-25",lab_df$time[lab_df$Whale == "I107"],"UTC")
lab_df$time[lab_df$Whale == "I129"] <- paste("2020-08-30",lab_df$time[lab_df$Whale == "I129"],"UTC")
lab_df$time[lab_df$Whale == "I145"] <- paste("2020-08-30",lab_df$time[lab_df$Whale == "I145"],"UTC")
lab_df$time[lab_df$Whale == "L87"] <- paste("2020-09-10",lab_df$time[lab_df$Whale == "L87"],"UTC")

lab_df$time <- as.POSIXct(lab_df$time,tz="UTC")

#lab_df <- lab_df[!grepl("?",lab_df$notes,fixed=T),]

for(rownum in 1:nrow(lab_df)){
  inds <- Data_fine$stime < lab_df$time[rownum] 
  inds <- inds & Data_fine$etime > lab_df$time[rownum]
  inds <- inds & Data_fine$ID == lab_df$Whale[rownum]
  inds[is.na(inds)] <- F
  Data_fine$label[inds] <- T
}

# do inference on the level of segments
Data_fine <- Data_fine[Data_fine$level == 2,]

# do some feature engineering
#features <- c("jp","rajp","ad",
#              "htv","ptv","rtv",
#              "aw1","aw2","aw3",
#              "w_low","w_high",
#              "pitch","VeDBA","abs_roll",
#              "label","camon")

features <- c("jp","rajp",
              "ad","delt_d",
              "logHtv","logPtv","logRtv",
              "pitch","abs_roll",
              "logWLow","logWHigh","logWTotal","VeDBA",
              "label","camon")

Data_fine_complete <- Data_fine[,features]
Data_fine_complete <- Data_fine_complete[complete.cases(Data_fine_complete),]
Data_fine_complete$label <- factor(Data_fine_complete$label)

# subset model
neg_inds <- which(Data_fine_complete$label != T & Data_fine_complete$camon)
pos_inds <- which(Data_fine_complete$label == T)
Data_fine_complete <- Data_fine_complete[c(neg_inds,pos_inds),
                                         names(Data_fine_complete) != "camon"]

# pick features
rfe_results <- rfe(Data_fine_complete[,names(Data_fine_complete) != "label"], 
                factor(Data_fine_complete$label), 
                rfeControl=rfeControl(functions=rfFuncs, method="cv", number=3))
print(rfe_results)
predictors(rfe_results)

# subset data for important features
Data_fine_complete <- Data_fine_complete[,c(predictors(rfe_results),"label")]

# pick weights
n <- sum(Data_fine_complete$label == T)
m <- sum(Data_fine_complete$label == F)
wT <- (n+m)/(2*n)
wF <- (n+m)/(2*m)
weights <- rep(wF,nrow(Data_fine_complete))
weights[Data_fine_complete$label == T] <- wT

# train rpart model
model <- train(Data_fine_complete[,names(Data_fine_complete) != "label"], 
               factor(Data_fine_complete$label), 
               method = "rpart",
               tuneGrid = data.frame(cp = c(0.01,0.1,1.0)),
               weights = weights)
print(varImp(model))
fancyRpartPlot(model$finalModel)

# train rf model
#model <- train(Data_fine_complete[,names(Data_fine_complete) != "label"], 
#               factor(Data_fine_complete$label), 
#               ntree = 100,
#               method = "rf")
#print(varImp(model))

for(feature in predictors(rfe_results)){
  print(ggplot(Data_fine_complete,
         aes_string(x = feature,
                    y = "..density..",
                    fill = "label")) +
    geom_histogram(alpha = 0.5, 
                   position = "identity"))
}

ggplot(Data_fine,
       aes(x = htv, y = w_low, fill = label)) +
  scale_x_continuous(limits = c(0, 2)) +
  scale_y_continuous(limits = c(0, 1500)) +
  scale_alpha_continuous(limits = c(-4, 0)) +
  stat_density_2d(aes(alpha = log10(..level..)), geom = "polygon", 
                  contour = T, bins = 200)

ggplot(Data_fine[Data_fine$label,],
       aes(x = htv, y = w_low)) +
  scale_x_continuous(limits = c(0, 2)) +
  scale_y_continuous(limits = c(0, 1500)) +
  scale_fill_continuous(limits = c(-6, 0)) +
  stat_density_2d(aes(fill = log10(..level..)), geom = "polygon", 
                  color = "white",contour = T)

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

# make camon with NA "none"
Data_fine$foraging_signs[is.na(Data_fine$foraging_signs) & 
                           Data_fine$camon == 1] <- "none"

if(plot_raw){
  
  # import raw Data
  rawData <- data.frame(fread('../../dat/Final_rawData_Beth.csv'))
  rawData <- rawData[rawData$divenum %in% dives,]
  rawData$elev <- -rawData$p
  
  # define color scheme
  colors <- hcl(h = seq(15, 375, length = 6),
                l = 65, c = 100)[1:5]
  colors <- c(colors,"#000000")
  
  colors <- c("echo.rapid" = colors[1],
              "chase" = colors[2],
              "crunch" = colors[3],
              "scales" = colors[4],
              "foraging" = colors[5],
              "none" = colors[6])
  
  cols_to_plot <- c("head","pitch","roll","elev")
  col_title <- "hpr"
  #cols_to_plot <- c("Aw_1","Aw_2","Aw_3","elev")
  #col_title <- "Aw"
  #cols_to_plot <- c("dyn.aw1","dyn.aw2","dyn.aw3","elev")
  #col_title <- "dyn.w"
  #cols_to_plot <- c("VeDBA","elev","d_elev")
  #col_title <- "elev"
  
  for(dive in dives){
    
    df <- rawData[rawData$divenum == dive,]
    
    w <- min(50,floor(nrow(df)/2))
    df$d_elev <- c(NA,diff(stats::filter(df$elev, rep(1 / w, w), sides = 2)))
    
    # add foraging signs to raw Data
    df$foraging_signs <- NA
    Data_fine0 <- Data_fine[Data_fine$divenum == dive & 
                              !is.na(Data_fine$foraging_signs),
                            c("stime","etime","foraging_signs")]
    
    for(seg_ind in 1:nrow(Data_fine0)){
      stime <- Data_fine0$stime[seg_ind]
      etime <- Data_fine0$etime[seg_ind]
      foraging_sign <- Data_fine0$foraging_signs[seg_ind]
      df$foraging_signs[df$Time >= stime & df$Time <= etime] <- foraging_sign
    }
    
    df <- df %>% pivot_longer(cols = cols_to_plot,
                              names_to = "feature")
    
    stime <- min(df$Time)
    #df$Time <- df$Time - stime
    
    plot1 <- ggplot(df,aes(x=Time, y=value, color=foraging_signs)) +
      geom_path(aes(group = 1)) +
      geom_hline(yintercept = 0) +
      scale_color_manual(values = colors) +
      scale_x_datetime(labels = scales::label_date("%H:%M:%S")) +
      labs(color="", y="",
           x="Local Time",
           title=paste(df$ID[1],dive,stime)) +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            text = element_text(size=16)) +
      facet_wrap(~feature, ncol = 1, scales = "free_y")
    
    ggsave(paste0("plt/",
                  df$ID[1],"_",
                  dive,"_",
                  col_title,
                  ".png"),
           plot1,
           width = 8, height = 6)
    
    print(plot1)
    
    for(axs in c()){
      
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
             zlab = paste(axs,dive,"max intensity",max_P),
             zlim = c(-12,12),
             drawPalette = T,
             decimate = F)
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

hmm <- fitHMM(data=Data_fine[Data_fine$ID %in% train_dives,],
              nbStates=N,
              dist=dist,
              DM=DM,
              beta0=beta0,
              delta0=delta0,
              Par0=Par0,
              fixPar=fixPar,
              knownStates=knownStates[Data_fine$ID %in% train_dives],
              lambda=0.5,
              #userBounds=userBounds,
              #workBounds=workBounds,
              nlmPar = list('iterlim'=1000,
                            'print.level'=2))

model <- paste0("EDA_hmm_",
                str(lambda),"_",
                str(test_dives[1]),".rds")
saveRDS(hmm,paste0(directory,"/params/",model))

### plot results

Data_fine0 <- Data_fine[Data_fine$ID %in% train_dives,]

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

# create a confusion matrix for the holdout whales 
