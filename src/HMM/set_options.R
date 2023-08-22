# first delete workspace
rm(list = ls())

# select number of initializations
n_retries <- 10

# Select Sex
sex <- c("Male","Female")

# set if we are doing hierarchical
hier <- TRUE
span <- 10 # minutes

# Select number of States
statesPerBehaviour <- c(3,3,3) # rest, trav, forg
N <- sum(statesPerBehaviour)

# Select distribution
dist <- list()
#dist[["logMaxDepth"]] <- "norm"
#dist[["logDiveDuration"]] <- "norm"
#dist[["logWLow"]] <- "norm"
#dist[["logWHigh"]] <- "norm"
#dist[["logWTotal"]] <- "norm"
#dist[["logW"]] <- "mvnorm2"
dist[["logMDDD"]] <- "mvnorm2"

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
dd_thresh <- 2 # seconds
md_thresh <- 0.5 # meters

# directory
directory <- "../../exp/hier_logMDDD_3_share"

# Set user and work bounds
userBounds <- list()
workBounds <- list()

for(feature in features1){
  if(dist[feature] == "mvnorm2"){
    userBounds[[feature]] <- matrix(c(rep(c(-Inf,-Inf,0.01,0.01,0.01),each=N),
                                      rep(c(Inf,Inf,Inf,Inf,Inf),each=N)),
                                    nrow=5*N,ncol=2)
    workBounds[[feature]] <- matrix(c(rep(c(-Inf,-Inf,-Inf,-Inf,-Inf),each=N),
                                      rep(c(Inf,Inf,Inf,-0.01,Inf),each=N)),
                                    nrow=5*N,ncol=2)
  } else if(dist[feature] == "norm"){
    userBounds[[feature]] <- matrix(c(rep(c(-Inf,0.01),each=N),
                                      rep(c(Inf,Inf),each=N)),
                                    nrow=2*N,ncol=2)
    workBounds[[feature]] <- matrix(c(rep(c(-Inf,-Inf),each=N),
                                      rep(c(Inf,Inf),each=N)),
                                    nrow=2*N,ncol=2)
  }
}

save.image(file = "options.RData")
