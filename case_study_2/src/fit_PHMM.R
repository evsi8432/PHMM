# set random seed
set.seed(rand_seed)

# define states
bhavs <- c("descent","bottom",
           "chase","crunch",
           "ascent1","ascent2")
N <- length(bhavs)

# set distributions
dist <- list()
dist[["delt_d"]] <- "norm"
dist[["htv"]] <- "gamma"
dist[["jp_normed"]] <- "gamma"

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
                     rbind(diag((N-1)),c(rep(0,(N-2)),1)))
  } else if (dist[[feature]] == "vm"){
    DM0 <- rbind(diag((N-1)),c(rep(0,(N-2)),1))
  } 
  DM[[feature]] <- DM0
}

# Set fixed Parameters
fixPar <- list()
eps <- 1e-50
fixPar$delt_d <- c(NA, 0, 0, 0,NA, # means each state
                   NA,NA,NA,NA,NA) # sds each state

# fix beta
fixPar$beta <- c(       NA,-1e2,-1e2,  NA,-1e2,  # descent
                 -1e2,       NA,-1e2,  NA,-1e2,  # bottom
                 -1e2,  NA,       NA,  NA,-1e2,  # chase
                 -1e2,-1e2,-1e2,     -1e2,  NA,  # crunch
                 -1e2,-1e2,-1e2,-1e2,     -1e2,  # ascent 1
                 -1e2,-1e2,-1e2,-1e2,-1e2      ) # ascent 2

fixPar$delta <- c(1.0-5*eps,eps,eps,eps,eps,eps)

# adjust standard deviation if rand_seed = 1
if(rand_seed == 1){
  sig <- 0
} else {
  sig <- 1
}

# set initial parameters
Par0 <- list()
Par0[["delt_d"]] <- c(c(5,  0,  0,  0, -5) + rnorm(5, sd = c(2*sig,0,0,0,2*sig)), # means
                  log(c(5,2.0,2.0,2.0,5)) + rnorm(5, sd = 0.5)) # sds

Par0[["htv"]] <- c(log(c(0.2,0.2,0.5,1.0,0.2)) + rnorm(5, sd = 0.5*sig), # means
                   log(c(0.1,0.1,0.4,1.0,0.1)) + rnorm(5, sd = 0.5*sig)) # sds

Par0[["jp_normed"]] <- c(log(c(10,5,10,15,10)) + rnorm(5, sd = 0.5*sig), # means
                         log(c(10,5,10,15,10)) + rnorm(5, sd = 0.5*sig)) # sds

# pick initial beta
beta0  <- c(       -2 + 1*sig*rnorm(1),-1e2,-1e2,  -2 + 1*sig*rnorm(1),-1e2,  # descent
            -1e2,       -2 + 1*sig*rnorm(1),-1e2,  -2 + 1*sig*rnorm(1),-1e2,  # bottom
            -1e2,  -2 + 1*sig*rnorm(1),       -2 + 1*sig*rnorm(1),  -2 + 1*sig*rnorm(1),-1e2,  # chase
            -1e2,-1e2,-1e2,     -1e2, 1*sig*rnorm(1),  # capture
            -1e2,-1e2,-1e2,-1e2,     -1e2,  # ascent 1
            -1e2,-1e2,-1e2,-1e2,-1e2      ) # ascent 2

beta0 <- matrix(beta0,nrow=1)

# pick initial delta
delta0 <- matrix(c(1.0-5*eps,eps,eps,eps,eps,eps),nrow=1)

# get knownStates for lambda purposes
knownStates <- df$knownState[df$ID %in% train_dives]
knownStates[knownStates %in% 7] <- NA

# fit the HMM
model_name <- paste0("hmm_",
                     paste0(setdiff(names(dist),"knownState"),collapse = "_"),"_",
                     log10(lambda),"_",
                     k,"_",K,"_",
                     rand_seed,".rds")

if(file.exists(paste0(directory,"/params/",model_name))){
  hmm <- readRDS(paste0(directory,"/params/",model_name))
} else {
  df <- prepData(df,coordNames=NULL)
  hmm <- fitHMM(data=df[df$ID %in% train_dives,],
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
                nlmPar = list('iterlim'=1000))
  
  saveRDS(hmm,paste0(directory,"/params/",model_name))
}

print(hmm)
