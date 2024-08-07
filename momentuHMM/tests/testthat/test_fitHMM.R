
context("fitHMM")

test_that("Exceptions are thrown",{
  data <- example$m$data
  simPar <- example$simPar
  par0 <- example$par0

  expect_error(fitHMM(data=data,nbStates=simPar$nbStates,Par=par0$Par,
                     beta0=par0$beta0,delta0=par0$delta0,formula=par0$formula,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean), NA)

  # if nbStates<1
  expect_error(fitHMM(data=data,nbStates=0,Par=par0$Par,
                     beta0=par0$beta0,delta0=par0$delta0,formula=par0$formula,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean))

  # if data empty
  expect_error(fitHMM(data=data.frame(),nbStates=simPar$nbStates,Par=par0$Par,
                     beta0=par0$beta0,delta0=par0$delta0,formula=par0$formula,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean))

  # if Par empty
  expect_error(fitHMM(data=data,nbStates=simPar$nbStates,Par=list(),
                     beta0=par0$beta0,delta0=par0$delta0,formula=par0$formula,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean))

  # if stepDist not in list
  expect_error(fitHMM(data=data,nbStates=simPar$nbStates,Par=par0$Par,
                     beta0=par0$beta0,delta0=par0$delta0,formula=par0$formula,
                     dist=list(step="unif",angle=simPar$dist$angle),estAngleMean=example$m$conditions$estAngleMean))

  # if angleDist not in list
  expect_error(fitHMM(data=data,nbStates=simPar$nbStates,Par=par0$Par,
                     beta0=par0$beta0,delta0=par0$delta0,formula=par0$formula,
                     dist=list(step=simPar$dist$step,angle="unif"),estAngleMean=example$m$conditions$estAngleMean))

  # if stepPar not within bounds
  expect_error(fitHMM(data=data,nbStates=simPar$nbStates,Par=list(step=-par0$Par$step,angle=par0$Par$angle),
                     beta0=par0$beta0,delta0=par0$delta0,formula=par0$formula,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean))

  # if wrong number of initial parameters
  expect_error(fitHMM(data=data,nbStates=simPar$nbStates,Par=list(step=par0$stepPar0[-1],angle=par0$anglePar0),
                     beta0=par0$beta0,delta0=par0$delta0,formula=par0$formula,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean))

  # if stepPar are missing
  expect_error(fitHMM(data=data,nbStates=simPar$nbStates,Par=par0$Par[-1],
                     beta0=par0$beta0,delta0=par0$delta0,formula=par0$formula,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean))

  # if beta0 has the wrong dimensions
  expect_error(fitHMM(data=data,nbStates=simPar$nbStates,Par=par0$Par,
                     beta0=par0$beta0[1,],delta0=par0$delta0,formula=par0$formula,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean))

  # if delta0 has the wrong length
  expect_error(fitHMM(data=data,nbStates=simPar$nbStates,Par=par0$Par,
                     beta0=par0$beta0,delta0=par0$delta0[-1],formula=par0$formula,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean))
  
  # invalid userBounds
  expect_error(fitHMM(data=data,nbStates=simPar$nbStates,Par=par0$Par,
                     beta0=par0$beta0,delta0=par0$delta0,formula=par0$formula,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean,
                     userBounds=list(step=matrix(c(-Inf,-Inf,0,0,Inf,Inf,Inf,Inf),simPar$nbStates*2,2))))
  
  # invalid DM
  expect_error(fitHMM(data=data,nbStates=simPar$nbStates,Par=list(step=log(par0$Par$step),angle=par0$Par$angle),
                     beta0=par0$beta0,delta0=par0$delta0,formula=par0$formula,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean,
                     DM=list(step=diag(3))))
  
  # invalid workBounds
  expect_error(fitHMM(data=data,nbStates=simPar$nbStates,Par=list(step=log(par0$Par$step),angle=par0$Par$angle),
                     beta0=par0$beta0,delta0=par0$delta0,formula=par0$formula,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean,
                     DM=list(step=diag(4)),workBounds=list(step=c(1,1,1))))
  
  expect_error(fitHMM(data=data,nbStates=simPar$nbStates,Par=list(step=log(par0$Par$step),angle=par0$Par$angle),
                     beta0=par0$beta0,delta0=par0$delta0,formula=par0$formula,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean,
                     DM=list(step=diag(4)),workBounds=list(step=matrix(c(0,0,0,0,0,0),3,2))))

})

test_that("The output has the right class",{
  data <- example$m$data
  simPar <- example$simPar
  par0 <- example$par0

  m <- fitHMM(data=data,nbStates=simPar$nbStates,Par=par0$Par,
                beta0=par0$beta0,delta0=par0$delta0,formula=par0$formula,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean)

  expect_equal(length(which(class(m)=="momentuHMM")),1)
})

test_that("Step length only + zero-inflation works",{
  
  oldRNG<-setRNG::setRNG()
  
  setRNG::setRNG(kind="Mersenne-Twister",normal.kind="Inversion",seed=1)
  
  #set.seed(1)
  nbAnimals <- 2
  nbStates <- 2
  nbCovs <- 2
  mu <- c(100,600)
  sigma <- c(10,40)
  zeromass <- c(0.1,0.05)
  stepPar <- c(mu,sigma,zeromass)
  anglePar <- NULL
  stepDist <- "gamma"
  angleDist <- "none"
  zeroInflation <- TRUE
  nbAnim <- c(50,100)

  data <- simData(nbAnimals=nbAnimals,nbStates=nbStates,dist=list(step=stepDist),
                  Par=list(step=stepPar),nbCovs=nbCovs,zeroInflation=list(step=zeroInflation),
                  obsPerAnimal=nbAnim)

  mu0 <- c(100,600)
  sigma0 <- c(10,40)
  zeromass0 <- c(0.1,0.05)
  stepPar0 <- c(mu0,sigma0,zeromass0)
  anglePar0 <- NULL
  angleMean <- NULL
  formula <- ~cov1+cov2

  expect_error(fitHMM(data=data,nbStates=nbStates,Par=list(step=c(log(stepPar0[1:(2*nbStates)]),stats::qlogis(zeromass0))),DM=list(step=diag(3*nbStates)),dist=list(step=stepDist),formula=formula,
              nlmPar=list(print.level=0)), NA)
  
  setRNG::setRNG(oldRNG)
})

test_that("equivalent momentuHMM and moveHMM models match",{

  oldRNG<-setRNG::setRNG()
  
  simPar <- example$simPar
  par0 <- example$par0
  nbStates<-simPar$nbStates
  
  setRNG::setRNG(kind="Mersenne-Twister",normal.kind="Inversion",seed=4)
  
  data<-simData(nbAnimals=2,model=example$m)
  momentuHMM_fit<-fitHMM(data=data,nbStates=nbStates,Par=list(step=log(par0$Par$step),angle=par0$Par$angle),stationary=TRUE,
                         beta0=par0$beta0[1,,drop=FALSE],delta0=par0$delta0,DM=list(step=diag(4)),dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean)
  moveHMM_fit<-moveHMM::fitHMM(moveData(data),nbStates=nbStates,stepPar=par0$Par$step,anglePar=par0$Par$angle,stepDist=simPar$dist$step,angleDist=simPar$dist$angle,stationary=TRUE,
                               beta0=par0$beta0[1,,drop=FALSE],delta0=par0$delta0)
  expect_equal(abs(momentuHMM_fit$mod$estimate-moveHMM_fit$mod$estimate)<1.e-6,rep(TRUE,length(momentuHMM_fit$mod$estimate)))
  expect_equal(abs(momentuHMM_fit$mod$minimum-moveHMM_fit$mod$minimum)<1.e-6,TRUE)
  
  #moveHMM called from within momentuHMM
  momentuHMM_fit2<-fitHMM(data=data,nbStates=nbStates,Par=list(step=par0$Par$step,angle=par0$Par$angle),stationary=TRUE,
                         beta0=par0$beta0[1,,drop=FALSE],delta0=par0$delta0,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean)
  expect_equal(abs(momentuHMM_fit2$mod$estimate-moveHMM_fit$mod$estimate)<1.e-6,rep(TRUE,length(momentuHMM_fit2$mod$estimate)))
  expect_equal(abs(momentuHMM_fit2$mod$minimum-moveHMM_fit$mod$minimum)<1.e-6,TRUE)
  
  #zeroInflation
  nbAnimals <- 2
  nbStates <- 2
  nbCovs <- 2
  mu <- c(100,600)
  sigma <- c(10,40)
  zeromass <- c(0.1,0.05)
  stepPar <- c(mu,sigma,zeromass)
  anglePar <- c(0,0,0.25,0.75)
  stepDist <- "gamma"
  angleDist <- "wrpcauchy"
  zeroInflation <- TRUE
  nbAnim <- c(50,100)
  beta0 <- matrix(-1.5,1,nbStates)
  
  data <- simData(nbAnimals=nbAnimals,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
                  Par=list(step=stepPar,angle=anglePar),nbCovs=nbCovs,zeroInflation=list(step=zeroInflation),
                  obsPerAnimal=nbAnim)
  
  mu0 <- c(100,600)
  sigma0 <- c(10,40)
  zeromass0 <- c(0.1,0.05)
  stepPar0 <- c(mu0,sigma0,zeromass0)
  anglePar0 <- anglePar
  stepDM <- diag(6)
  Par0<-getParDM(data,nbStates,dist=list(step=stepDist,angle=angleDist),Par=list(step=stepPar0,angle=anglePar0),zeroInflation=list(step=TRUE),estAngleMean=list(angle=TRUE),
                 DM=list(step=stepDM))
  
  momentuHMM_fit<-fitHMM(data=data,nbStates=nbStates,Par0=Par0,stationary=TRUE,
                         beta0=beta0,delta0=example$par0$delta0,DM=list(step=stepDM),dist=list(step=stepDist,angle=angleDist),estAngleMean=list(angle=TRUE))
  moveHMM_fit<-moveHMM::fitHMM(moveData(data),nbStates=nbStates,stepPar=stepPar0,anglePar=anglePar0,stationary=TRUE,
                               beta0=beta0,delta0=example$par0$delta0,stepDist=stepDist,angleDist=angleDist)
  expect_equal(abs(momentuHMM_fit$mod$estimate-moveHMM_fit$mod$estimate)<1.e-6,rep(TRUE,length(momentuHMM_fit$mod$estimate)))
  expect_equal(abs(momentuHMM_fit$mod$minimum-moveHMM_fit$mod$minimum)<1.e-6,TRUE)
  
  #moveHMM called from within momentuHMM
  momentuHMM_fit2<-fitHMM(data=data,nbStates=nbStates,Par0=list(step=stepPar0,angle=anglePar0),stationary=TRUE,
                         beta0=beta0,delta0=example$par0$delta0,dist=list(step=stepDist,angle=angleDist),estAngleMean=list(angle=TRUE))
  expect_equal(abs(momentuHMM_fit2$mod$estimate-moveHMM_fit$mod$estimate)<1.e-6,rep(TRUE,length(momentuHMM_fit2$mod$estimate)))
  expect_equal(abs(momentuHMM_fit2$mod$minimum-moveHMM_fit$mod$minimum)<1.e-6,TRUE)
  
  setRNG::setRNG(oldRNG)
  
})

test_that("equivalent betaRef models match",{
  
  oldRNG<-setRNG::setRNG()
  
  par0 <- getPar0(example$m)
  nbStates<-length(example$m$stateNames)
  
  setRNG::setRNG(kind="Mersenne-Twister",normal.kind="Inversion",seed=1)
  
  data<-simData(nbAnimals=2,model=example$m)
  momentuHMM_fit1<-fitHMM(data=data,nbStates=nbStates,Par=list(step=log(par0$Par$step),angle=par0$Par$angle),stationary=TRUE,
                          DM=list(step=diag(4)),dist=example$m$conditions$dist,estAngleMean=example$m$conditions$estAngleMean)
  momentuHMM_fit2<-fitHMM(data=data,nbStates=nbStates,Par=list(step=log(par0$Par$step),angle=par0$Par$angle),stationary=TRUE,
                          DM=list(step=diag(4)),dist=example$m$conditions$dist,estAngleMean=example$m$conditions$estAngleMean,betaRef=c(2,2))
  expect_equal(all(abs(momentuHMM_fit2$mod$estimate[-9]-momentuHMM_fit1$mod$estimate[-9])<1.e-6),TRUE)
  expect_equal(abs(momentuHMM_fit2$mod$minimum-momentuHMM_fit1$mod$minimum)<1.e-6,TRUE)
  
  setRNG::setRNG(oldRNG)
  
})

test_that("equivalent models with and without dummy covariate match",{
  
  oldRNG<-setRNG::setRNG()
  
  simPar <- example$simPar
  par0 <- example$par0
  nbStates<-simPar$nbStates
  
  setRNG::setRNG(kind="Mersenne-Twister",normal.kind="Inversion",seed=1)
  
  data<-simData(nbAnimals=2,model=example$m,obsPerAnimal = 30)
  data$cov1<-rep(1,nrow(data)) # dummy covariate
  momentuHMM_nocov<-fitHMM(data=data,nbStates=nbStates,Par=list(step=log(par0$Par$step),angle=par0$Par$angle),
                         beta0=par0$beta0[1,,drop=FALSE],delta0=par0$delta0,DM=list(step=diag(4)),dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean)
  DM1<-list(step=list(mean=~cov1,sd=~cov1),angle=list(mean=~cov1,concentration=~cov1))
  Par0_1<-getParDM(data.frame(cov1=1),nbStates,simPar$dist,par0$Par,DM=DM1,estAngleMean=example$m$conditions$estAngleMean)
  momentuHMM_cov1<-fitHMM(data=data,nbStates=nbStates,Par=Par0_1,formula=~cov1,
                           beta0=rbind(par0$beta0[1,,drop=FALSE]/2,par0$beta0[1,,drop=FALSE]/2),delta0=par0$delta0,DM=DM1,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean)
  DM2<-list(step=matrix(c(1,"cov1",0,0,0,0,0,0,
                          0,0,1,"cov1",0,0,0,0,
                          0,0,0,0,1,"cov1",0,0,
                          0,0,0,0,0,0,1,"cov1"),nrow=4,byrow=TRUE),
            angle=matrix(c(1,"cov1",0,0,0,0,0,0,
                           0,0,1,"cov1",0,0,0,0,
                           0,0,0,0,1,"cov1",0,0,
                           0,0,0,0,0,0,1,"cov1"),nrow=4,byrow=TRUE))
  Par0_2<-getParDM(data.frame(cov1=1),nbStates,simPar$dist,par0$Par,DM=DM2,estAngleMean=example$m$conditions$estAngleMean)
  momentuHMM_cov2<-fitHMM(data=data,nbStates=nbStates,Par=Par0_2,formula=~cov1,
                          beta0=rbind(par0$beta0[1,,drop=FALSE]/2,par0$beta0[1,,drop=FALSE]/2),delta0=par0$delta0,DM=DM2,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean)
  
  expect_equal(abs(unlist(lapply(momentuHMM_nocov$CIreal,function(x) x$est),use.names=FALSE)-unlist(lapply(momentuHMM_cov1$CIreal,function(x) x$est),use.names=FALSE))<1.e-3,rep(TRUE,length(unlist(lapply(momentuHMM_nocov$CIreal,function(x) x$est),use.names=FALSE))))
  expect_equal(abs(momentuHMM_nocov$mod$minimum-momentuHMM_cov1$mod$minimum)<1.e-6,TRUE)
  expect_equal(abs(unlist(lapply(momentuHMM_cov2$CIreal,function(x) x$est),use.names=FALSE)-unlist(lapply(momentuHMM_cov1$CIreal,function(x) x$est),use.names=FALSE))<1.e-3,rep(TRUE,length(unlist(lapply(momentuHMM_nocov$CIreal,function(x) x$est),use.names=FALSE))))
  expect_equal(abs(momentuHMM_cov2$mod$minimum-momentuHMM_cov1$mod$minimum)<1.e-6,TRUE)
  
  # circular-circular regression
  data$cov2<-rep(0,nrow(data)) # dummy covariate
  DM=list(angle=list(mean=~1,concentration=~1))
  Par0<-getParDM(data,nbStates,simPar$dist,par0$Par,DM=DM,estAngleMean=example$m$conditions$estAngleMean,circularAngleMean = list(angle=TRUE))
  momentuHMM_nocov<-fitHMM(data=data,nbStates=nbStates,Par=Par0,
                           beta0=par0$beta0[1,,drop=FALSE],delta0=par0$delta0,DM=DM,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean,circularAngleMean = list(angle=TRUE))
  DM1<-list(step=list(mean=~cov2,sd=~cov2),angle=list(mean=~cov2,concentration=~cov2))
  Par0_1<-getParDM(data.frame(cov2=0),nbStates,simPar$dist,par0$Par,DM=DM1,estAngleMean=example$m$conditions$estAngleMean,circularAngleMean = list(angle=TRUE))
  momentuHMM_cov1<-fitHMM(data=data,nbStates=nbStates,Par=Par0_1,formula=~cov2,
                          beta0=rbind(par0$beta0[1,,drop=FALSE],c(0,0)),delta0=par0$delta0,DM=DM1,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean,circularAngleMean = list(angle=TRUE))
  DM2<-list(step=matrix(c(1,"cov2",0,0,0,0,0,0,
                          0,0,1,"cov2",0,0,0,0,
                          0,0,0,0,1,"cov2",0,0,
                          0,0,0,0,0,0,1,"cov2"),nrow=4,byrow=TRUE),
            angle=matrix(c("cov2",0,0,0,0,0,
                           0,"cov2",0,0,0,0,
                           0,0,1,"cov2",0,0,
                           0,0,0,0,1,"cov2"),nrow=4,byrow=TRUE,dimnames=list(c(paste0("mean_",1:nbStates),paste0("concentration_",1:nbStates)),
                                                                             c(paste0("mean_",1:nbStates,":(cov2)"),paste0("concentration_",rep(1:nbStates,each=2),rep(c(":(Intercept)",":cov2"),2))))))
  Par0_2<-getParDM(data.frame(cov2=0),nbStates,simPar$dist,par0$Par,DM=DM2,estAngleMean=example$m$conditions$estAngleMean,circularAngleMean = list(angle=TRUE))
  momentuHMM_cov2<-fitHMM(data=data,nbStates=nbStates,Par=Par0_2,formula=~cov2,
                          beta0=rbind(par0$beta0[1,,drop=FALSE],c(0,0)),delta0=par0$delta0,DM=DM2,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean,circularAngleMean = list(angle=TRUE))
  
  DM3<-list(step=matrix(c(1,"cov2",0,0,0,0,0,0,
                          0,0,1,"cov2",0,0,0,0,
                          0,0,0,0,1,"cov2",0,0,
                          0,0,0,0,0,0,1,"cov2"),nrow=4,byrow=TRUE),
            angle=matrix(c("angleFormula(cov2,cov1)",0,0,0,0,0,
                           0,"angleFormula(cov2,cov1)",0,0,0,0,
                           0,0,1,"cov2",0,0,
                           0,0,0,0,1,"cov2"),nrow=4,byrow=TRUE,dimnames=list(c(paste0("mean_",1:nbStates),paste0("concentration_",1:nbStates)),
                                                                             c(paste0("mean_",1:nbStates,":cov1:(cov2)"),paste0("concentration_",rep(1:nbStates,each=2),rep(c(":(Intercept)",":cov2"),2))))))
  momentuHMM_cov3<-fitHMM(data=data,nbStates=nbStates,Par=Par0_2,formula=~cov2,
                          beta0=rbind(par0$beta0[1,,drop=FALSE],c(0,0)),delta0=par0$delta0,DM=DM3,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean,circularAngleMean = list(angle=TRUE))

  DM4<-list(step=matrix(c(1,"cov2",0,0,0,0,0,0,
                          0,0,1,"cov2",0,0,0,0,
                          0,0,0,0,1,"cov2",0,0,
                          0,0,0,0,0,0,1,"cov2"),nrow=4,byrow=TRUE),
            angle=list(mean=~angleFormula(cov2,cov1),concentration=~cov2))
  
  momentuHMM_cov4<-fitHMM(data=data,nbStates=nbStates,Par=Par0_2,formula=~cov2,
                          beta0=rbind(par0$beta0[1,,drop=FALSE],c(0,0)),delta0=par0$delta0,DM=DM4,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean,circularAngleMean = list(angle=TRUE))

  DM5<-list(step=matrix(c(1,"cov2",0,0,0,0,0,0,
                        0,0,1,"cov2",0,0,0,0,
                        0,0,0,0,1,"cov2",0,0,
                        0,0,0,0,0,0,1,"cov2"),nrow=4,byrow=TRUE),
          angle=list(mean=~angleFormula(cov2,cov1,by=ID),concentration=~cov2))
    
  Par0_5 <- Par0_2
  Par0_5$angle <- c(1,1,1,1,0,0,0,0)
  
  momentuHMM_cov5<-fitHMM(data=data,nbStates=nbStates,Par=Par0_5,formula=~cov2,
                          beta0=rbind(par0$beta0[1,,drop=FALSE],c(0,0)),delta0=par0$delta0,DM=DM5,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean,circularAngleMean = list(angle=TRUE))
    
  DM6<-list(step=matrix(c(1,"cov2",0,0,0,0,0,0,
                            0,0,1,"cov2",0,0,0,0,
                            0,0,0,0,1,"cov2",0,0,
                            0,0,0,0,0,0,1,"cov2"),nrow=4,byrow=TRUE),
              angle=matrix(c("angleFormula(cov2,cov1,by=ID1)","angleFormula(cov2,cov1,by=ID2)",0,0,0,0,0,0,
                             0,0,"angleFormula(cov2,cov1,by=ID1)","angleFormula(cov2,cov1,by=ID2)",0,0,0,0,
                             0,0,0,0,1,"cov2",0,0,
                             0,0,0,0,0,0,1,"cov2"),nrow=4,byrow=TRUE,dimnames=list(c(paste0("mean_",1:nbStates),paste0("concentration_",1:nbStates)),
                                                                               c(paste0("mean_",rep(1:nbStates,each=2),":ID",rep(1:2,2),":cov1:(cov2)"),paste0("concentration_",rep(1:nbStates,each=2),rep(c(":(Intercept)",":cov2"),2))))))
  
  momentuHMM_cov6<-fitHMM(data=data,nbStates=nbStates,Par=Par0_5,formula=~cov2,
                          beta0=rbind(par0$beta0[1,,drop=FALSE],c(0,0)),delta0=par0$delta0,DM=DM6,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean,circularAngleMean = list(angle=TRUE))

  expect_equal(abs(unlist(lapply(momentuHMM_nocov$CIreal,function(x) x$est),use.names=FALSE)-unlist(lapply(momentuHMM_cov1$CIreal,function(x) x$est),use.names=FALSE))<1.e-3,rep(TRUE,length(unlist(lapply(momentuHMM_nocov$CIreal,function(x) x$est),use.names=FALSE))))
  expect_equal(abs(momentuHMM_nocov$mod$minimum-momentuHMM_cov1$mod$minimum)<1.e-6,TRUE)
  expect_equal(abs(unlist(lapply(momentuHMM_cov2$CIreal,function(x) x$est),use.names=FALSE)-unlist(lapply(momentuHMM_cov1$CIreal,function(x) x$est),use.names=FALSE))<1.e-3,rep(TRUE,length(unlist(lapply(momentuHMM_nocov$CIreal,function(x) x$est),use.names=FALSE))))
  expect_equal(abs(momentuHMM_cov2$mod$minimum-momentuHMM_cov1$mod$minimum)<1.e-6,TRUE)
  expect_equal(abs(unlist(lapply(momentuHMM_cov3$CIreal,function(x) x$est),use.names=FALSE)-unlist(lapply(momentuHMM_cov1$CIreal,function(x) x$est),use.names=FALSE))<1.e-3,rep(TRUE,length(unlist(lapply(momentuHMM_nocov$CIreal,function(x) x$est),use.names=FALSE))))
  expect_equal(abs(momentuHMM_cov3$mod$minimum-momentuHMM_cov1$mod$minimum)<1.e-6,TRUE)
  expect_equal(abs(unlist(lapply(momentuHMM_cov4$CIreal,function(x) x$est),use.names=FALSE)-unlist(lapply(momentuHMM_cov1$CIreal,function(x) x$est),use.names=FALSE))<1.e-3,rep(TRUE,length(unlist(lapply(momentuHMM_nocov$CIreal,function(x) x$est),use.names=FALSE))))
  expect_equal(abs(momentuHMM_cov4$mod$minimum-momentuHMM_cov1$mod$minimum)<1.e-6,TRUE)
  expect_equal(abs(unlist(lapply(momentuHMM_cov5$CIreal,function(x) x$est),use.names=FALSE)-unlist(lapply(momentuHMM_cov1$CIreal,function(x) x$est),use.names=FALSE))<1.e-3,rep(TRUE,length(unlist(lapply(momentuHMM_nocov$CIreal,function(x) x$est),use.names=FALSE))))
  expect_equal(abs(momentuHMM_cov5$mod$minimum-momentuHMM_cov1$mod$minimum)<1.e-6,TRUE)
  expect_equal(abs(unlist(lapply(momentuHMM_cov6$CIreal,function(x) x$est),use.names=FALSE)-unlist(lapply(momentuHMM_cov1$CIreal,function(x) x$est),use.names=FALSE))<1.e-3,rep(TRUE,length(unlist(lapply(momentuHMM_nocov$CIreal,function(x) x$est),use.names=FALSE))))
  expect_equal(abs(momentuHMM_cov6$mod$minimum-momentuHMM_cov1$mod$minimum)<1.e-6,TRUE)
  
  setRNG::setRNG(oldRNG)
  
})

test_that("retryFits works",{
  data <- example$m$data
  simPar <- example$simPar
  par0 <- example$par0
  par0$Par$step <- c(60.39849, 2316.51195,   70.13456,   18.51266)
  par0$Par$angle <- c(-0.1811843,  1.0969285,  1.1738157,  2.5099845)
  par0$beta0 <- matrix(c(-50.397203, -11.051239,
                        4.751330,-13.734779,
                        -2.345075,  -1.895257),3,2,byrow=TRUE)
  par0$delta0 <- c(1-3.774067e-14,3.774067e-14)
  
  setRNG::setRNG(kind="Mersenne-Twister",normal.kind="Inversion",seed=3)
  mod1 <- fitHMM(data=data,nbStates=simPar$nbStates,Par=par0$Par,
                  beta0=par0$beta0,delta0=par0$delta0,formula=par0$formula,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean)
  
  Par0 <- getPar(mod1)
  mod2 <- fitHMM(data=data,nbStates=simPar$nbStates,Par=Par0$Par,
                 beta0=Par0$beta,delta0=Par0$delta,formula=par0$formula,dist=simPar$dist,estAngleMean=example$m$conditions$estAngleMean,retryFits=5,retrySD=list(step=5,beta=20,delta=20))
  expect_equal(isTRUE(all.equal(mod1$mod$estimate,mod2$mod$estimate)),FALSE)
  
})