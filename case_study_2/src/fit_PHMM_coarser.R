# define states
bhavs <- c("descent","bottom",
           "chase","crunch",
           "ascent1","ascent2")
N <- length(bhavs)

# set distributions
dist <- list()
if(lambda == -1){
  dist[["knownState"]] <- paste0("cat",N+1)
}
dist[["delt_d"]] <- "norm"
#dist[["logHtv"]] <- "norm"
#dist[["logJpNorm"]] <- "norm"
dist[["rajp"]] <- "vm"
dist[["htv"]] <- "gamma"
dist[["jp_normed"]] <- "gamma"
#dist[["logWLow"]] <- "norm"
#dist[["logWHigh"]] <- "norm"
#dist[["eating"]] <- "cat3" # crunch, scales, nothing
#dist[["echo"]] <- "cat3" # steady, rapid, nothing
#dist[["forage"]] <- "cat3" # chase, capture, nothing

# make the DM matrix (share features for ascending with and without a fish)
DM <- list()
for(feature in names(dist)){
  if(feature %in% c("delt_d")){
    DM0 <- kronecker(diag(2),
                     rbind(diag((N-1)),c(rep(0,(N-2)),1)))
  } else if (feature == "knownState") {
    DM0 <- diag(N^2)
  } else if(dist[[feature]] %in% c("norm","gamma")){
    DM0 <- kronecker(diag(2),
                     rbind(c(1,0,0),diag(3),c(1,0,0),c(1,0,0)))
    #DM0 <- kronecker(diag(2),
    #                 rbind(diag((N-1)),c(rep(0,(N-2)),1)))
  } else if (dist[[feature]] == "vm"){
    DM0 <- rbind(c(1,0,0),diag(3),c(1,0,0),c(1,0,0))
    #DM0 <- rbind(diag((N-1)),c(rep(0,(N-2)),1))
  } 
  DM[[feature]] <- DM0
}

# Set fixed Parameters
fixPar <- list()

eps <- 1e-50

fixPar$delt_d <- c(NA, 0, 0, 0,NA, # means each state
                   NA,NA,NA,NA,NA) # sds each state

if(lambda == -1){
  fixPar$knownState <- c(  -99,-99,-99,-99,-99,-99, # prob desc label, each state
                           -99,-99,-99,-99,-99,-99, # prob bot label, each state
                           -99,-99,-99,-99,-99,-99, # prob chase label, each state
                           -99,-99, NA, NA,-99, NA, # prob capture label, each state
                           -99,-99,-99,-99, NA,-99, # prob asc 1 label, each state
                           -99,-99,-99,-99,-99, NA) # prob asc 2 label, each state
}

# fix beta
fixPar$beta <- c(       NA,-1e2,-1e2,  NA,-1e2,  # descent
                        -1e2,       NA,-1e2,  NA,-1e2,  # bottom
                        -1e2,  NA,       NA,  NA,-1e2,  # chase
                        -1e2,-1e2,-1e2,     -1e2, 1e2,  # crunch
                        -1e2,-1e2,-1e2,-1e2,     -1e2,  # ascent 1
                        -1e2,-1e2,-1e2,-1e2,-1e2      ) # ascent 2

fixPar$delta <- c(1.0-5*eps,eps,eps,eps,eps,eps)

# set initial parameters
Par0 <- list()
Par0[["delt_d"]] <- c(c(20,  0,   0,  0,-20), # means
                  log(c(10, 2.0,2.0,2.0,10))) # sds

Par0[["logWLow"]] <- c(c(2, 5, 6), # means
                   log(c(2, 1, 1))) # sds

Par0[["logWHigh"]] <- c(c(1, 4, 6), # means
                    log(c(2, 1, 1))) # sds

Par0[["htv"]] <- c(log(c(0.2,0.5,2.0)), # means
                   log(c(0.1,0.4,0.5))) # sds

#Par0[["htv"]] <- c(log(c(0.2,0.2,0.5,2.0,0.2)), # means
#                   log(c(0.1,0.1,0.4,0.5,0.1))) # sds

Par0[["rajp"]] <- c(log(c(2, 1, 0.1)))  # concentrations

Par0[["logJpNorm"]] <- c(c(1, 1, 2,   4, 1), # means
                     log(c(1, 1, 1, 0.5, 1))) # sds

Par0[["jp_normed"]] <- c(log(c(3,5,25)), # means
                         log(c(2,15,5))) # sds
#Par0[["jp_normed"]] <- c(log(c(3,3,5,25,3)), # means
#                         log(c(2,2,15,5,2))) # sds

if(lambda == -1){
  Par0[["knownState"]] <- c(-99,-99,-99,-99,-99,-99, # prob desc label, each state
                            -99,-99,-99,-99,-99,-99, # prob bot label, each state
                            -99,-99,-99,-99,-99,-99, # prob chase label, each state
                            -99,-99,  0,  0,-99,  0, # prob capture label, each state
                            -99,-99,-99,-99,  0,-99, # prob asc 1 label, each state
                            -99,-99,-99,-99,-99,  0) # prob asc 2 label, each state
}


# pick initial beta
beta0  <- c(       -3,-1e2,-1e2,  -3,-1e2,  # descent
                   -1e2,       -3,-1e2,  -3,-1e2,  # bottom
                   -1e2,  -3,       -3,  -3,-1e2,  # chase
                   -1e2,-1e2,-1e2,     -1e2, 1e2,  # crunch
                   -1e2,-1e2,-1e2,-1e2,     -1e2,  # ascent 1
                   -1e2,-1e2,-1e2,-1e2,-1e2      ) # ascent 2

beta0 <- matrix(beta0,nrow=1)

# pick initial delta
delta0 <- matrix(c(1.0-5*eps,eps,eps,eps,eps,eps),nrow=1)

# prep data
Data_fine$whale <- Data_fine$ID
Data_fine$ID <- Data_fine$divenum
Data_fine_final <- prepData(Data_fine[Data_fine$divenum %in% c(train_dives,test_dives),],
                            coordNames=NULL)

# turn Data_fine into 10 second windows

window_size <- 10 # seconds- must be multiple of 2
n <- window_size / 2

knownStates <- c()
jps_normed <- c()
htvs <- c()
rajps <- c()
delt_ds <- c()
divenums <- c()
stimes <- c()
ads <- c()

logWLows <- c()
logWHighs <- c()

for(divenum in unique(Data_fine_final$divenum)){
  print(divenum)
  dive_df <- Data_fine_final[Data_fine_final$divenum %in% divenum,]
  
  if(nrow(dive_df) < 2*n){
    next
  }
  
  inds <- seq(1, nrow(dive_df)-n, by = n)
  for(rownum in inds){
    
    if(rownum == tail(inds,1)){
      seg_df <- dive_df[rownum:nrow(dive_df),]
    } else {
      seg_df <- dive_df[rownum:(rownum+n-1),]
    }
    
    if(4 %in% seg_df$knownState){
      knownStates <- c(knownStates,4)
    } else if (5 %in% seg_df$knownState){
      knownStates <- c(knownStates,5)
    } else if (6 %in% seg_df$knownState){
      knownStates <- c(knownStates,6)
    } else {
      knownStates <- c(knownStates,NA)
    }
    
    max_jp_normed <- max(seg_df$jp_normed,na.rm=T)
    jps_normed <- c(jps_normed,max_jp_normed)
    jp_ind <- which(seg_df$jp_normed == max_jp_normed)[1]
    rajps <- c(rajps,seg_df$rajp[jp_ind])
    
    htvs <- c(htvs,mean(seg_df$htv,na.rm=T))
    delt_ds <- c(delt_ds,sum(seg_df$delt_d))
    stimes <- c(stimes,seg_df$stime[1])
    ads <- c(ads,seg_df$ad[1])
    logWLows <- c(logWLows, mean(seg_df$logWLow))
    logWHighs <- c(logWHighs, mean(seg_df$logWHigh))
    divenums <- c(divenums,divenum)
  }
}

Data_less_fine_final <- data.frame(knownState = knownStates,
                                   jp_normed = jps_normed,
                                   htv = htvs,
                                   rajp = rajps,
                                   delt_d = delt_ds,
                                   divenum = divenums,
                                   stime = stimes,
                                   ad = ads,
                                   logWLow = logWLows,
                                   logWHigh = logWHighs)

Data_less_fine_final$ID <- Data_less_fine_final$divenum

Data_less_fine_final <- prepData(Data_less_fine_final,coordNames=NULL)

checkPar0(data=Data_less_fine_final,
          nbStates=N,
          dist=dist,
          DM=DM,
          beta0=beta0,
          delta0=delta0,
          Par0=Par0,
          fixPar=fixPar
)

# get knownStates for lambda purposes
knownStates <- Data_less_fine_final$knownState[Data_less_fine_final$ID %in% train_dives]
knownStates[knownStates %in% 7] <- NA

# add category for knownStates Data_fine_final
Data_less_fine_final$knownState[is.na(Data_less_fine_final$knownState)] <- N+1

# fit the HMM
if(lambda == -1){
  hmm <- fitHMM(data=Data_less_fine_final[Data_less_fine_final$ID %in% train_dives,],
                nbStates=N,
                dist=dist,
                DM=DM,
                beta0=beta0,
                delta0=delta0,
                Par0=Par0,
                fixPar=fixPar,
                lambda=1.0,
                stateNames = bhavs,
                nlmPar = list('iterlim'=1000,
                              'print.level'=2))
} else {
  hmm <- fitHMM(data=Data_less_fine_final[Data_less_fine_final$ID %in% train_dives,],
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
                nlmPar = list('iterlim'=1000,
                              'print.level'=2))
}

print(hmm)

model_name <- paste0("hmm_",
                     paste0(setdiff(names(dist),"knownState"),collapse = "_"),"_",
                     log10(lambda),"_",
                     k,".rds")

saveRDS(hmm,paste0(directory,"/params/",model_name))