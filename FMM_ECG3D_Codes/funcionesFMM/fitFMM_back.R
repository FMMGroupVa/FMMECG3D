###########################################################################################
# Internal function: fit multicomponent FMM model
# Arguments:
#   vData: data to be fitted an FMM model.
#   nback: number of FMM components to be fitted.
#   timePoints: one single period time points.
#   maxiter: maximum number of iterations for the backfitting algorithm.
#   stopFunction: function to check the criterion convergence for the backfitting algorithm.
#   lengthAlphaGrid, lengthOmegaGrid: precision of the grid of alpha and omega parameters.
#   alphaGrid, omegaGrid: grids of alpha and omega parameters.
#   omegaMax: max value for omega.
#   numReps: number of times the alpha-omega grid search is repeated.
#   showProgress: TRUE to display a progress indicator on the console.
# Returns an object of class FMM.
###########################################################################################
fitFMM_back<-function(vData, nback, timePoints = seqTimes(length(vData)),
                      maxiter = nback, stopFunction = alwaysFalse,
                      lengthAlphaGrid = 48, lengthOmegaGrid = 24,
                      alphaGrid = seq(0, 2*pi, length.out = lengthAlphaGrid),
                      omegaMin = 0.0001, omegaMax = 1,
                      omegaGrid = exp(seq(log(omegaMin), log(omegaMax),
                                          length.out = lengthOmegaGrid)),
                      numReps = 3, showProgress = TRUE, usedApply = getApply(FALSE)[[1]]){

  nObs <- length(vData)

  if(showProgress){
    totalMarks <- 50
    partialMarkLength <- 2
    progressHeader<-paste(c("|", rep("-", totalMarks), "|\n|"), collapse = "")
    cat(progressHeader)
    completedPercentage <- 0.00001
    previousPercentage <- completedPercentage
  }

  # Object initialization.
  fittedValuesPerComponent <- matrix(rep(0, nObs*nback), ncol = nback)
  fittedFMMPerComponent <- list()
  prevFittedFMMvalues <- NULL
  stopCriteria <- "Stopped by reaching maximum iterations ("

  # Backfitting algorithm: iteration
  for(i in 1:maxiter){
    # Backfitting algorithm: component
    for(j in 1:nback){
      # data for component j: difference between vData and all other components fitted values
      backFittingData <- vData - apply(as.matrix(fittedValuesPerComponent[,-j]), 1, sum)

      # component j fitting using fitFMM_unit function
      fittedFMMPerComponent[[j]] <- fitFMM_unit(backFittingData, timePoints = timePoints, lengthAlphaGrid = lengthAlphaGrid,
                                            lengthOmegaGrid = lengthOmegaGrid, alphaGrid = alphaGrid, omegaMin = omegaMin,
                                            omegaMax = omegaMax, omegaGrid = omegaGrid, numReps = numReps, usedApply)
      fittedValuesPerComponent[,j] <- getFittedValues(fittedFMMPerComponent[[j]])
      # showProgress
      if(showProgress){
        completedPercentage <- completedPercentage + 100/(nback*maxiter)
        if(ceiling(previousPercentage) < floor(completedPercentage)){
          progressDone <- paste(rep("=",sum((seq(ceiling(previousPercentage), floor(completedPercentage), by = 1)
                                           %% partialMarkLength == 0))), collapse = "")
          cat(progressDone)
          previousPercentage <- completedPercentage
        }
      }
    }

    # Check stop criterion
    # Fitted values as sum of all components
    fittedFMMvalues <- apply(fittedValuesPerComponent, 1, sum)

    if(!is.null(prevFittedFMMvalues)){
      if(PV(vData, prevFittedFMMvalues) > PV(vData, fittedFMMvalues)){
        fittedFMMPerComponent <- previousFittedFMMPerComponent
        fittedFMMvalues <- prevFittedFMMvalues
        stopCriteria <- "Stopped by reaching maximum R2 ("
        break
      }
      if(stopFunction(vData, fittedFMMvalues, prevFittedFMMvalues)){
        stopCriteria <- "Stopped by the stopFunction ("
        break
      }
    }
    prevFittedFMMvalues <- fittedFMMvalues
    previousFittedFMMPerComponent <- fittedFMMPerComponent
  }
  nIter <- i

  # showProgress
  if(showProgress){
    if(completedPercentage < 100){
      completedPercentage <- 100
      nMarks <- ifelse(ceiling(previousPercentage) < floor(completedPercentage),
                       sum((seq(ceiling(previousPercentage),
                                floor(completedPercentage), by = 1) %% partialMarkLength == 0)), 0)
      if (nMarks > 0) {
        cat(paste(rep("=",nMarks), collapse = ""))
        previousPercentage <- completedPercentage
      }
    }
    cat("|\n", paste(stopCriteria, nIter, sep = ""),"iteration(s))","\n")
  }

  alpha <- unlist(lapply(fittedFMMPerComponent, getAlpha))
  beta <- unlist(lapply(fittedFMMPerComponent, getBeta))
  omega <- unlist(lapply(fittedFMMPerComponent, getOmega))

  # A and M estimates are recalculated by linear regression; cosPhi is the design matrix
  cosPhi <- calculateCosPhi(alpha = alpha, beta = beta, omega = omega, timePoints = timePoints)
  linearModel <- lm(vData ~ cosPhi)
  M <- as.vector(linearModel$coefficients[1])
  A <- as.vector(linearModel$coefficients[-1])
  fittedFMMvalues <- predict(linearModel)

  # Residual sum of squares
  SSE <- sum((fittedFMMvalues - vData)^2)

  # Reorder waves by explained variability
  iPV <- sapply(1:nback, function(x){PV(vData, predict(lm(vData ~ cosPhi[,x])))})
  waveOrder<-which.max(iPV)
  while(length(waveOrder)<nback){
     iPV <- sapply((1:nback)[-waveOrder], function(x){PV(vData, predict(lm(vData ~ cosPhi[,c(waveOrder,x)])))})
     waveOrder<-c(waveOrder,(1:nback)[-waveOrder][which.max(iPV)])
  }
  A <- A[waveOrder]
  alpha <- alpha[waveOrder]
  beta <- beta[waveOrder]
  omega <- omega[waveOrder]

  # Returns an object of class FMM.
  return(FMM(
    M = M,
    A = A,
    alpha = alpha,
    beta = beta,
    omega = omega,
    timePoints = timePoints,
    summarizedData = vData,
    fittedValues= fittedFMMvalues,
    SSE = SSE,
    R2 = PVj(vData, timePoints, alpha, beta, omega),
    nIter = nIter
  ))
}
