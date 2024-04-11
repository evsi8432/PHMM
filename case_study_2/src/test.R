#.rs.restartR()
library(Rcpp)
library(tools)
compileAttributes("/Users/evsi8432/Documents/Research/momentuHMM")
package_native_routine_registration_skeleton("/Users/evsi8432/Documents/Research/momentuHMM")
install.packages("/Users/evsi8432/Documents/Research/momentuHMM",
                 repos=NULL,
                 type="source")
library(momentuHMM)
#detach("package:momentuHMM", unload=TRUE)
#library(momentuHMM)

if(F){
#library(dplyr)
#library(tidyr)
#library(mclust)
#library(data.table)
#library(ggplot2)
#library(GGally)
#library(signal)
#library(oce)

setwd("/Users/evsi8432/Documents/Research/PHMM/src/case_study_2")

# set seed
set.seed(1)

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
}

T0 <- 3

x <- c(1,2,3,  # x coordinates
       1,2,3) # y coordinates

mu <- rbind(rep(2,T0), # x mean
            rep(2,T0)) # y mean

sig <- rbind(rep(1,T0), #sig xx
             rep(0,T0), #sig xy
             rep(0,T0), #sig yx
             rep(1,T0)) #sig yy

momentuHMM:::dmvnorm_rcpp(x,mu,sig)
