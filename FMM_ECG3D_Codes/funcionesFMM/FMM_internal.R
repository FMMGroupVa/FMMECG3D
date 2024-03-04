################################################################################
# Auxiliary internal functions
# Functions:
#   step1FMM:        M, A and beta initial parameter estimations.
#   bestStep1:       to find the optimal initial parameters estimation.
#   step2FMM:        second step of FMM fitting process.
#   refineFMM:       fitFMM from a previous objectFMM.
#   PV:              percentage of variability explained.
#   PVj:             percentage of variability explained by each component of
#                    FMM model.
#   seqTimes:        to build a sequence of equally time points spaced in range
#                    [0,2*pi).
#   calculateCosPhi: to calculate components' cos(phi(t)).
#   getApply:        returns the parallelized apply function depending on the OS.
################################################################################


################################################################################
# Internal function: to estimate M, A and beta initial parameters
# also returns residual sum of squared (RSS).
# Arguments:
#    alphaOmegaParameters: vector of the parameters alpha and omega
#    vData: data to be fitted an FMM model.
#    timePoints: one single period time points.
# Returns a 6-length numerical vector: M, A, alpha, beta, omega and RSS
################################################################################
step1FMM <- function(optBase, vData) {
  
  pars <- optBase[["base"]] %*% vData
  mobiusRegression <- pars[1] + pars[2]*optBase[["cost"]] + pars[3]*optBase[["sint"]]
  
  residualSS <- sum((vData - mobiusRegression)^2)/length(vData)
  
  return(c(optBase[["alpha"]], optBase[["omega"]], residualSS))
}

################################################################################
# Internal function: to find the optimal initial parameter estimation
# Arguments:
#    vData: data to be fitted an FMM model.
#    step1: a data.frame with estimates of
#           M, A, alpha, beta, omega, RSS as columns.
# Returns the optimal row of step1 argument.
# optimum: minimum RSS with several stability conditions.
################################################################################
bestStep1 <- function(vData, step1){

  # step1 in decreasing order by RSS
  orderedModelParameters <- order(step1[,"RSS"])

  maxVData <- max(vData)
  minVData <- min(vData)
  nObs <- length(vData)

  # iterative search: go through rows ordered step 1
  #    until the first one that verifies the stability conditions
  bestModelFound <- FALSE
  i <- 1
  while(!bestModelFound){
    # parameters
    mParameter <- step1[orderedModelParameters[i], "M"]
    aParameter <- step1[orderedModelParameters[i], "A"]
    sigma <- sqrt(step1[orderedModelParameters[i], "RSS"]*nObs/(nObs-5))

    # stability conditions
    amplitudeUpperBound <- mParameter + aParameter
    amplitudeLowerBound <- mParameter - aParameter
    rest1 <- amplitudeUpperBound <= maxVData + 1.96*sigma
    rest2 <- amplitudeLowerBound >= minVData - 1.96*sigma

    # it is necessary to check that there are no NA,
    # because it can be an extreme solution
    if(is.na(rest1)) rest1 <- FALSE
    if(is.na(rest2)) rest2 <- FALSE

    if(rest1 & rest2){
      bestModelFound <- TRUE
    } else {
      i <- i+1
    }

    if(i > nrow(step1))
      return(NULL)
  }

  return(step1[orderedModelParameters[i],])
}



################################################################################
# Internal function: second step of FMM fitting process
# Arguments:
#   parameters: M, A, alpha, beta, omega initial parameter estimations
#   vData: data to be fitted an FMM model.
#   timePoints: one single period time points.
#   omegaMax: max value for omega.
################################################################################
step2FMM <- function(parameters, vData, timePoints, omegaMax){

  nObs <- length(timePoints)

  # FMM model and residual sum of squares
  modelFMM <- parameters[1] + parameters[2] *
    cos(parameters[4]+2*atan2(parameters[5]*sin((timePoints - parameters[3])/2),
                                cos((timePoints - parameters[3])/2)))
  residualSS <- sum((modelFMM - vData)^2)/nObs
  sigma <- sqrt(residualSS*nObs/(nObs - 5))

  # When amplitude condition is valid, it returns RSS
  # else it returns infinite.
  amplitudeUpperBound <- parameters[1] + parameters[2]
  amplitudeLowerBound <- parameters[1] - parameters[2]
  rest1 <- amplitudeUpperBound <= max(vData) + 1.96*sigma
  rest2 <- amplitudeLowerBound >= min(vData) - 1.96*sigma

  # Other integrity conditions that must be met
  rest3 <- parameters[2] > 0  # A > 0
  rest4 <- parameters[5] > 0  &  parameters[5] <= omegaMax # omega in (0, omegaMax]
  if(rest1 & rest2 & rest3 & rest4)
    return(residualSS)
  else
    return(Inf)
}

#EXPERIMENTAL ---------------------
step2FMMExp <- function(par, vData, varErr = 1){
  return(-logLikFMMWave(alphaOmegaParameters = par, vData = vData, varErr))
}

################################################################################
# Internal function: to calculate the percentage of variability explained by
#   the FMM model
# Arguments:
#   vData: data to be fitted an FMM model.
#   pred: fitted values.
################################################################################
PV <- function(vData, pred){
  return(1 - sum((vData - pred)^2)/sum((vData - mean(vData))^2))
}

################################################################################
# Internal function: to calculate the percentage of variability explained by
#   each component of FMM model
# Arguments:
#   vData: data to be fitted an FMM model.
#   timePoints: one single period time points.
#   alpha, beta, omega: vectors of corresponding parameter estimates.
################################################################################
PVj <- function(vData, timePoints, alpha, beta, omega){
  # Fitted values of each wave
  nComponents <- length(alpha)
  waves <- calculateCosPhi(alpha = alpha, beta = beta, omega = omega,
                           timePoints = timePoints)

  # The percentage of variability explained up to wave i is determined
  cumulativePV <- sapply(1:nComponents, function(x){PV(vData, predict(lm(vData ~ waves[,1:x])))})

  # individual percentage of variability is the part that adds to the whole
  return(c(cumulativePV[1], diff(cumulativePV)))
}

################################################################################
# Internal function: to build a sequence of equally time points spaced
#                    in range [0,2*pi).
# Arguments:
#   nObs: secuence length.
################################################################################
seqTimes <- function(nObs){
  return(seq(0, 2*pi, length.out = nObs+1)[1:nObs])
}

################################################################################
# Internal function: to calculate components' cos(phi(t)).
# Arguments:
#   alpha, beta, omega: parameters.
#   timePoints: time points in which the FMM model is computed.
# Returns a matrix of each component's cos(phi(t)) as columns.
################################################################################
calculateCosPhi <- function(alpha, beta, omega, timePoints){
  calculateSingleCosPhi <- function(alpha, beta, omega){
    return(cos(beta + 2*atan(omega*tan((timePoints - alpha)/2))))
  }
  return(mapply(FUN = calculateSingleCosPhi, alpha = alpha, beta = beta, omega = omega))
}

################################################################################
# Internal function: returns the parallelized apply function depending on the OS.
# Returns the apply function to be used.
################################################################################
getApply <- function(parallelize = FALSE){

  getApply_Rbase <- function(){
    usedApply <- function(FUN, X, ...) t(apply(X = X, MARGIN = 1, FUN = FUN, ...))
  }

  getParallelApply_Windows <- function(parallelCluster){
    usedApply <- function(FUN, X, ...) t(parallel::parApply(parallelCluster, FUN = FUN,
                                                            X = X, MARGIN = 1, ...))
    return(usedApply)
  }

  parallelFunction_Unix<-function(nCores){
    # A parallelized apply function does not exist, so it must be translated to a lapply
    usedApply <- function(FUN, X, ...){
      matrix(unlist(parallel::mclapply(X = asplit(X, 1), FUN = FUN, mc.cores = nCores, ...)),
             nrow = nrow(X), byrow = T)
    }
    return(usedApply)
  }

  nCores <- min(12, parallel::detectCores() - 1)

  if(parallelize){
    # different ways to implement parallelization depending on OS:
    if(.Platform$OS.type == "windows"){
      parallelCluster <- parallel::makePSOCKcluster(nCores)
      doParallel::registerDoParallel(parallelCluster)
      usedApply <- getParallelApply_Windows(parallelCluster)
    }else{
      usedApply <- parallelFunction_Unix(nCores)
      parallelCluster <- NULL
    }
  }else{
    # R base apply:
    usedApply <- getApply_Rbase()
    parallelCluster <- NULL
  }

  return(list(usedApply, parallelCluster))
}
