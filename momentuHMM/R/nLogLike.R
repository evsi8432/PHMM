
#' Negative log-likelihood function
#'
#' @param optPar Vector of working parameters.
#' @param nbStates Number of states of the HMM.
#' @param formula Regression formula for the transition probability covariates.
#' @param bounds Named list of 2-column matrices specifying bounds on the natural (i.e, real) scale of the probability
#' distribution parameters for each data stream.
#' @param parSize Named list indicating the number of natural parameters of the data stream probability distributions
#' @param data An object \code{momentuHMMData}.
#' @param dist Named list indicating the probability distributions of the data streams.
#' @param covs data frame containing the beta model covariates (if any)
#' @param estAngleMean Named list indicating whether or not to estimate the angle mean for data streams with angular
#' distributions ('vm' and 'wrpcauchy').
#' @param circularAngleMean Named list indicating whether to use circular-linear or circular-circular
#' regression on the mean of circular distributions ('vm' and 'wrpcauchy') for turning angles. See \code{\link{fitHMM}}.
#' @param consensus Named list indicating whether to use the circular-circular regression consensus model
#' @param zeroInflation Named list of logicals indicating whether the probability distributions of the data streams are zero-inflated.
#' @param oneInflation Named list of logicals indicating whether the probability distributions of the data streams are one-inflated.
#' @param stationary \code{FALSE} if there are time-varying covariates in \code{formula} or any covariates in \code{formulaDelta}. If \code{TRUE}, the initial distribution is considered
#' equal to the stationary distribution. Default: \code{FALSE}.
#' @param fullDM Named list containing the full (i.e. not shorthand) design matrix for each data stream.
#' @param DMind Named list indicating whether \code{fullDM} includes individual- and/or temporal-covariates for each data stream
#' specifies (-1,1) bounds for the concentration parameters instead of the default [0,1) bounds.
#' @param Bndind Named list indicating whether \code{DM} is NULL with default parameter bounds for each data stream.
#' @param knownStates Vector of values of the state process which are known prior to fitting the
#' model (if any).
#' @param lambda is a scalar that weights un-labelled examples.
#' @param fixPar Vector of working parameters which are assumed known prior to fitting the model (NA indicates parameters is to be estimated).
#' @param wparIndex Vector of indices for the elements of \code{fixPar} that are not NA.
#' @param nc indicator for zeros in fullDM
#' @param meanind index for circular-circular regression mean angles with at least one non-zero entry in fullDM
#' @param covsDelta data frame containing the delta model covariates (if any)
#' @param workBounds named list of 2-column matrices specifying bounds on the working scale of the probability distribution, transition probability, and initial distribution parameters
#' @param prior A function that returns the log-density of the working scale parameter prior distribution(s)
#' @param betaCons Matrix of the same dimension as \code{beta0} composed of integers identifying any equality constraints among the t.p.m. parameters.
#' @param betaRef Indices of reference elements for t.p.m. multinomial logit link.
#' @param deltaCons Matrix of the same dimension as \code{delta0} composed of integers identifying any equality constraints among the initial distribution working scale parameters.
#' @param optInd indices of constrained parameters
#' @param recovs data frame containing the recharge model theta covariates (if any)
#' @param g0covs data frame containing the recharge model g0 covariates (if any)
#' @param mixtures Number of mixtures for the state transition probabilities
#' @param covsPi data frame containing the pi model covariates
#' @param recharge recharge model specification (only used for hierarchical models)
#' @param aInd vector of indices of first observation for each animal
#'
#' @return The negative log-likelihood of the parameters given the data.
#'
#' @examples
#' \dontrun{
#' # data is a momentuHMMData object (as returned by prepData), automatically loaded with the package
#' data <- example$m$data
#' m<-example$m
#' Par <- getPar(m)
#' nbStates <- length(m$stateNames)
#'
#' inputs <- momentuHMM:::checkInputs(nbStates,m$conditions$dist,Par$Par,m$conditions$estAngleMean,
#'           m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,
#'           m$conditions$DM,m$conditions$userBounds,
#'           m$stateNames)
#'
#' wpar <- momentuHMM:::n2w(Par$Par,m$conditions$bounds,list(beta=Par$beta),
#'         log(Par$delta[-1]/Par$delta[1]),nbStates,m$conditions$estAngleMean,
#'         m$conditions$DM,m$conditions$Bndind,
#'         m$conditions$dist)
#'
#' l <- momentuHMM:::nLogLike(wpar,nbStates,m$conditions$formula,m$conditions$bounds,
#'      inputs$p$parSize,data,inputs$dist,model.matrix(m$conditions$formula,data),
#'      m$conditions$estAngleMean,m$conditions$circularAngleMean,inputs$consensus,
#'      m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$stationary,
#'      m$conditions$fullDM,m$conditions$DMind,
#'      m$conditions$Bndind,m$knownStates,m$lambda,unlist(m$conditions$fixPar),
#'      m$conditions$wparIndex,covsDelta=m$covsDelta,workBounds=m$conditions$workBounds,
#'      betaRef=m$conditions$betaRef,covsPi=m$covsPi)
#' }
#'

nLogLike <- function(optPar,nbStates,formula,bounds,parSize,data,dist,covs,
                     estAngleMean,circularAngleMean,consensus,zeroInflation,oneInflation,
                     stationary=FALSE,fullDM,DMind,Bndind,knownStates,lambda=1.0,fixPar,wparIndex,nc,meanind,covsDelta,workBounds,prior=NULL,betaCons=NULL,betaRef,deltaCons=NULL,optInd=NULL,recovs=NULL,g0covs=NULL,mixtures=1,covsPi,recharge=NULL,aInd)
{

  # check arguments
  distnames<-names(dist)

  #covs <- stats::model.matrix(formula,data)
  nbCovs <- ncol(covs)-1
  if(!is.null(recovs)) {
    nbRecovs <- ncol(recovs)-1
    nbG0covs <- ncol(g0covs)-1
  } else nbRecovs <- nbG0covs <- 0

  # convert the parameters back to their natural scale
  wpar <- expandPar(optPar,optInd,fixPar,wparIndex,betaCons,deltaCons,nbStates,ncol(covsDelta)-1,stationary,nbCovs,nbRecovs+nbG0covs,mixtures,ncol(covsPi)-1)
  par <- w2n(wpar,bounds,parSize,nbStates,nbCovs,estAngleMean,circularAngleMean,consensus,stationary,fullDM,DMind,nrow(data),dist,Bndind,nc,meanind,covsDelta,workBounds,covsPi)

  if(nbRecovs){
    for(i in 1:length(unique(data$ID))){
      idInd <- which(data$ID==unique(data$ID)[i])
      if(inherits(data,"hierarchical")) {
        recLevels <- length(recharge)
        recLevelNames <- names(recharge)
        rechargeNames <- paste0("recharge",gsub("level","",recLevelNames))
        colInd <- lapply(recLevelNames,function(x) which(grepl(paste0("I((level == \"",gsub("level","",x),"\")"),colnames(recovs),fixed=TRUE)))
      } else {
        recLevels <- 1
        rechargeNames <- "recharge"
        colInd <- list(1:ncol(recovs))
      }
      for(iLevel in 1:recLevels){
        g0 <- par$g0 %*% t(g0covs[(i-1)*recLevels+iLevel,,drop=FALSE])
        theta <- par$theta
        covs[idInd,grepl(rechargeNames[iLevel],colnames(covs))] <- cumsum(c(g0,theta[colInd[[iLevel]]]%*%t(recovs[idInd[-length(idInd)],colInd[[iLevel]]])))
      }
    }
  }

  if(is.null(knownStates)) knownStates <- -1
  else knownStates[which(is.na(knownStates))] <- 0

  # NULL arguments don't suit C++
  if(any(unlist(lapply(dist,is.null)))){
    par[which(unlist(lapply(dist,is.null)))]<-matrix(NA)
  }

  if(stationary)
    par$delta <- c(NA)
  if(nbStates==1) {
    par$beta <- matrix(NA)
    #par[["pi"]] <- c(NA)
    par$delta <- c(NA)
    par[distnames] <- lapply(par[distnames],as.matrix)
  }

  nllk <- tryCatch(nLogLike_rcpp(nbStates,as.matrix(covs),data,names(dist),dist,
                        par,
                        aInd,zeroInflation,oneInflation,stationary,
                        knownStates,lambda,
                        betaRef,mixtures),error=function(e) e)

  if(inherits(nllk,"error")) nllk <- NaN

  if(!is.null(prior)) nllk <- nllk - prior(wpar)

  return(nllk)
}
