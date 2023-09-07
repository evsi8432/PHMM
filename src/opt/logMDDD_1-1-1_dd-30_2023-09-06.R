#setwd("/Users/evsi8432/Documents/Research/PHMM/src/bash")

# directory to save things to
directory <- "../../exp/logMDDD_1-1-1_dd-30_2023-09-06"
dir.create(directory, showWarnings = FALSE)
dir.create(paste0(directory,"/params"), showWarnings = FALSE)
dir.create(paste0(directory,"/plt"), showWarnings = FALSE)

# select number of initializations
n_retries <- 10
iterlim <- 1000

# Select Sex
sex <- c("Male","Female")

# set if we are doing hierarchical
hier <- FALSE
share_fine <- FALSE
span <- 5 # minutes

# Select number of States
statesPerBehaviour <- c(1,1,1) # rest, trav, forg
N <- sum(statesPerBehaviour)
N0 <- statesPerBehaviour[1]

# Select distribution
dist <- list()
#dist[["logMaxDepth"]] <- "norm"
#dist[["logDiveDuration"]] <- "norm"
#dist[["logWLow"]] <- "norm"
#dist[["logWHigh"]] <- "norm"
#dist[["logWTotal"]] <- "norm"
#dist[["logW"]] <- "mvnorm2"
dist[["logMDDD"]] <- "mvnorm2"
#dist[["maxDepthCat"]] <- "cat5"

# hold on to features that will be need in the future
features1 <- names(dist)
features2 <- c()
for(feature in names(dist)){
  if(dist[[feature]] == "mvnorm2"){
    features2 <- c(features2,
                   paste0(feature,".x"),
                   paste0(feature,".y"))
  } else {
    features2 <- c(features2,feature)
  }
}

# Select dive duration threshold
dd_thresh <- 30 # seconds
md_thresh <- 0.5 # meters
md_threshs <- c(5,10,30,50) # meters

# Set user and work bounds
userBounds <- list()
workBounds <- list()

for(feature in features1){
  if(dist[[feature]] == "mvnorm2"){
    userBounds[[feature]] <- matrix(c(rep(c(-Inf,-Inf,0.01,0.01,0.01),each=N),
                                      rep(c(Inf,Inf,Inf,Inf,Inf),each=N)),
                                    nrow=5*N,ncol=2)
    workBounds[[feature]] <- matrix(c(rep(c(-Inf,-Inf,-Inf,-Inf,-Inf),each=N),
                                      rep(c(Inf,Inf,Inf,-0.01,Inf),each=N)),
                                    nrow=5*N,ncol=2)
  } else if (dist[[feature]] == "norm"){
    userBounds[[feature]] <- matrix(c(rep(c(-Inf,0.01),each=N),
                                      rep(c(Inf,Inf),each=N)),
                                    nrow=2*N,ncol=2)
    workBounds[[feature]] <- matrix(c(rep(c(-Inf,-Inf),each=N),
                                      rep(c(Inf,Inf),each=N)),
                                    nrow=2*N,ncol=2)
  } else if (substring(dist[[feature]], 1,3) == "cat") {
    ncats <- as.integer(substring(dist[[feature]], 4))

    userBounds[[feature]] <- matrix(c(rep(0.00,each=(ncats-1)*N),
                                      rep(1.00,each=(ncats-1)*N)),
                                    nrow=(ncats-1)*N,ncol=2)
    workBounds[[feature]] <- matrix(c(rep(-Inf,each=(ncats-1)*N),
                                      rep(Inf,each=(ncats-1)*N)),
                                    nrow=(ncats-1)*N,ncol=2)
  }
}
