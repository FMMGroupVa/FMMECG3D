#### Dependencies ####
require("RColorBrewer")
require("FMM")
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(tibble))
suppressPackageStartupMessages(require(tidyr))
suppressWarnings(require(writexl))

#### Fit of multiFMM model functions ####
fitMultiFMM <- function(vDataMatrix, nBack = 5, maxIter = 10, weightError = TRUE, 
                        lengthAlphaGrid = 48, lengthOmegaGrid = 24, 
                        alphaGrid = seq(0, 2*pi, length.out = lengthAlphaGrid),
                        omegaMin = 0.0001, omegaMax = 1,
                        omegaGrid = exp(seq(log(max(omegaMin, omegaMin)), log(1), length.out = lengthOmegaGrid)),
                        parallelize = TRUE){
  
  nSignals <- ncol(vDataMatrix)
  
  #### Preparations for the fit ####
  
  # Preparation of the paralellization
  usedApply_Cluster <- getApply(parallelize = parallelize)
  usedApply <- usedApply_Cluster[[1]]
  
  # The results are stored in:
  paramsPerWave <- replicate(nSignals, simplify = FALSE, 
                             data.frame(M=double(), A=double(), Alpha=double(), Beta=double(), Omega=double()))
  fittedWaves<-lapply(1:nSignals,function(x){
    vData<-vDataMatrix[,x]; vData<-vData[!is.na(vData)];
    return(array(0,c(length(vData),nBack)))})
  
  # The MSE must be balanced in the fit of multiple signals
  errorWeights<-rep(1,nSignals); totalMSE<-rep(Inf, maxIter)
  
  # Backfitting algorithm: iteration
  currentBack<-1; continueBackfitting<-TRUE
  cat(paste("Backfitting: ",sep=""))
  
  while(continueBackfitting) {
    cat(paste(currentBack," ",sep=""))
    fmmObjectList<-list()
    paramsPerWave <- replicate(nSignals, simplify = FALSE,
                               data.frame(M=double(), A=double(),
                                          Alpha=double(), Beta=double(), Omega=double()))
    # Backfitting algorithm: component
    for(j in 1:nBack){
      #### First Step: determine optimal common parameters (alpha and omega) ####
      optimalParams <- optimizeAlphaOmega(vDataMatrix = vDataMatrix, fittedWaves = fittedWaves,
                                          currentComp = j, alphaGrid = alphaGrid,
                                          omegaGrid = omegaGrid, omegaMax = omegaMax,
                                          errorWeights = errorWeights, usedApply = usedApply)
      
      #### Second Step: fit single FMM wave in each signal with common parameters ####
      for(signalIndex in 1:nSignals){
        vData<-vDataMatrix[,signalIndex]; vData<-vData[!is.na(vData)]; nObs<-length(vData)
        vData<-vData - apply(as.matrix(fittedWaves[[signalIndex]][,-j]), 1, sum)
        fmmObjectList[[signalIndex]]<-optimizeOtherParameters(vData = vData, fixedParams = optimalParams)
        
        fittedWaves[[signalIndex]][,j]<-getFittedValues(fmmObjectList[[signalIndex]])
      }
      
      #### Get FMM parameters per wave ####
      params<-data.frame("M"=sapply(fmmObjectList,getM),"A"=sapply(fmmObjectList,getA),"Alpha"=sapply(fmmObjectList,getAlpha),
                         "Beta"=sapply(fmmObjectList,getBeta),"Omega"=sapply(fmmObjectList,getOmega), "Var"=rep(NA,nSignals))
      paramsPerWave<-lapply(1:nSignals, function(x) rbind(paramsPerWave[[x]], params[x,]))
      paramsPerWave<-lapply(1:nSignals, function(x) recalculateMA(vDatai=vDataMatrix[,x], paramsPerSignal=paramsPerWave[[x]]))
      #### Error is weighted across signals ####
      sigma<-sapply(1:nSignals, function(x){sum((vDataMatrix[,x]-apply(as.matrix(fittedWaves[[x]]), 1, sum))^2)/(nObs-1)})
      if(weightError) errorWeights<-1/sigma
      else errorWeights<-rep(1,nSignals)
    }
    totalMSE[currentBack]<-sum(sapply(1:nSignals, function(x){sum((vDataMatrix[,x]-apply(as.matrix(fittedWaves[[x]]), 1, sum))^2)}))
    
    # Check if backfitting should stop
    stopCondition1<-currentBack>=maxIter
    stopCondition2<-ifelse(currentBack>1, totalMSE[currentBack]>totalMSE[currentBack-1], FALSE)
    if(stopCondition1 | stopCondition2){
      stopCriteria<-ifelse(stopCondition1, "Maximum Iterations", "Minimum MSE")
      cat(paste("  Stop: ", stopCriteria, "\n", sep=""))
      continueBackfitting<-FALSE
    }else{
      currentBack<-currentBack+1
    }
  }

  # Unname waves and stop parallelized cluster
  for(i in 1:nSignals) rownames(paramsPerWave[[i]])<-1:nBack
  cluster <- usedApply_Cluster[[2]]
  if(!is.null(cluster)) parallel::stopCluster(cluster)

  #### Return results ####
  return(paramsPerWave)
}


#### Internal multiFMM functions ####
## MultiFMM, first step: optimize common parameters

step1FMM3D <- function(alphaOmegaParameters, vData, timePoints, sigmas = rep(1, ncol(vData))) {
  alpha <- as.numeric(alphaOmegaParameters[1])
  omega <- as.numeric(alphaOmegaParameters[2])
  tStar <- alpha + 2*atan2(omega*sin((timePoints-alpha)/2), cos((timePoints-alpha)/2))
  
  fittedValues <- apply(vData, 2, function(x){
    dM <- cbind(rep(1, length(x)), cos(tStar), sin(tStar))
    mDeltaGamma <- solve(crossprod(dM), t(dM)%*%x)
    yFit <- mDeltaGamma[1] + mDeltaGamma[2]*cos(tStar) + mDeltaGamma[3]*sin(tStar)
    return(yFit)
  })
  residuales <- (vData - fittedValues)^2
  logL <- -sum(t(sigmas^2*t(residuales))/2)
  return(c(alpha, omega, logL))
}

logLik3DFMM1 <- function(alphaOmegaParameters, vData, timePoints, sigmas = rep(1, ncol(vData))) {
  alpha <- as.numeric(alphaOmegaParameters[1])
  omega <- as.numeric(alphaOmegaParameters[2])
  tStar <- alpha + 2*atan2(omega*sin((timePoints-alpha)/2), cos((timePoints-alpha)/2))
  fittedValues <- apply(vData, 2, function(x){
    dM <- cbind(rep(1, length(x)), cos(tStar), sin(tStar))
    mDeltaGamma <- solve(t(dM)%*%dM)%*%t(dM)%*%x
    yFit <- mDeltaGamma[1] + mDeltaGamma[2]*cos(tStar) + mDeltaGamma[3]*sin(tStar)
    return(yFit)
  })
  residuales <- (vData - fittedValues)^2
  logL <- -sum(t(sigmas^2*t(residuales))/2) # Cada columna por su sigma^2
  return(logL)
}

optimizeAlphaOmega<-function(vDataMatrix, fittedWaves, currentComp,
                             alphaGrid, omegaGrid, errorWeights, usedApply,
                             omegaMax = 0.7){
  nObs<-nrow(vDataMatrix)
  nSignals <- ncol(vDataMatrix)
  timePoints = seqTimes(nObs)
  
  grid <- expand.grid(alphaGrid, omegaGrid)
  residualsMatrix <- vDataMatrix
  
  rssMatrix<-matrix(NA, nrow = nrow(grid), ncol = (nSignals)+2)
  rssMatrix[,1]<-grid[,1]; rssMatrix[,2]<-grid[,2]
  
  for(signalIndex in 1:nSignals){
    vData<-vDataMatrix[,signalIndex]
    residualsMatrix[,signalIndex] <- vData - apply(as.matrix(fittedWaves[[signalIndex]][,-currentComp]), 1, sum)
    
    step1 <- usedApply(FUN = FMM:::step1FMM, X = grid, vData = residualsMatrix[,signalIndex],
                       timePoints = timePoints)
    rssMatrix[,signalIndex+2] <- errorWeights[signalIndex]*step1[,6]
  }
  
  bestParamsIndex <- which.min(apply(rssMatrix[,-c(1:2)],1,mean))
  initialParams <- as.numeric(rssMatrix[bestParamsIndex, c(1:2)])
  
  nelderMead <- optim(par = initialParams, fn = logLik3DFMM1, 
                      vData = residualsMatrix, timePoints = timePoints, sigmas = 1/errorWeights, 
                      method = "L-BFGS-B", control = list(fnscale = -1), lower = c(0, omegaGrid[1]), upper = c(2*pi, 1))
  
  return(nelderMead$par) # Best alpha and omega
  
}

## MultiFMM, second step: determine non-common parameters
optimizeOtherParameters <- function(vData, fixedParams, timePoints = seqTimes(length(vData))){
  
  alpha <- as.numeric(fixedParams[1])
  omega <- as.numeric(fixedParams[2])
  tStar <- alpha + 2*atan2(omega*sin((timePoints-alpha)/2), cos((timePoints-alpha)/2))
  
  dM <- cbind(rep(1, length(tStar)), cos(tStar), sin(tStar))
  mDeltaGamma <- solve(t(dM)%*%dM)%*%t(dM)%*%vData
  M <- mDeltaGamma[1]
  delta <- mDeltaGamma[2]
  gamma <- mDeltaGamma[3]
  
  fittedFMMvalues <- mDeltaGamma[1] + mDeltaGamma[2]*cos(tStar) + mDeltaGamma[3]*sin(tStar)
  SSE <- sum((vData-fittedFMMvalues)^2)
  
  fmmObject<-FMM(
    M = mDeltaGamma[1],
    A = sqrt(delta^2+gamma^2),
    alpha = alpha,
    beta = atan2(-gamma, delta) + alpha,
    omega = omega,
    timePoints = timePoints,
    summarizedData = vData,
    fittedValues = fittedFMMvalues,
    SSE = SSE,
    R2 = PV(vData, fittedFMMvalues),
    nIter = 0
  )
  
  fmmObject@nPeriods<-1
  fmmObject@data<-vData
  
  return(fmmObject)
}

step2FMM_fixedAlphaOmega <- function(parameters, vData, timePoints, fixedAlpha, fixedOmega){
  nObs <- length(timePoints)
  modelFMM <- parameters[1] + parameters[2] *
    cos(parameters[4]+2*atan2(fixedOmega*sin((timePoints - fixedAlpha)/2),
                              cos((timePoints - fixedAlpha)/2)))
  residualSS <- sum((modelFMM - vData)^2)/nObs
  rest3 <- parameters[2] > 0  # A > 0
  if(rest3)
    return(residualSS)
  else
    return(Inf)
}

#### Other multiFMM useful functions ####
recalculateMA<-function(vDatai, paramsPerSignal){
  
  nObs<-length(vDatai)
  
  ## Extract Params. and calculate cosPhi
  alpha<-paramsPerSignal$Alpha; beta<-paramsPerSignal$Beta
  omega<-paramsPerSignal$Omega; maxComp<-length(alpha)
  cosPhi<-t(calculateCosPhi(alpha = alpha, beta = beta,
                            omega = omega, timePoints = seqTimes(nObs)))
  
  ## Recalculate M and As
  currentFormula<-as.formula(paste("vDatai~", paste("cosPhi[",1:maxComp,",]",sep="",collapse="+")))
  regression<-lm(currentFormula)
  recalculatedM<-regression$coefficients[1]
  recalculatedAs<-as.vector(regression$coefficients[-1])
  
  ## Calculate R2, including other waves
  individualR2<-PVj(vData = vDatai, timePoints = seqTimes(nObs),
                    alpha = paramsPerSignal$Alpha,
                    beta = paramsPerSignal$Beta, omega = paramsPerSignal$Omega)
  paramsPerSignal$Var<-individualR2
  
  mResults<-rep(NA,maxComp); aResults<-rep(NA,maxComp)
  mResults<-rep(recalculatedM, maxComp); aResults<-recalculatedAs
  paramsPerSignal$M<-mResults; paramsPerSignal$A<-aResults
  
  ## Negative A => Translation on the beta parameter
  aValidContion <- aResults>0 | is.na(aResults)
  if(any(!aValidContion)){
    paramsPerSignal$Beta[!aValidContion]<-(paramsPerSignal$Beta[!aValidContion]+pi)%%(2*pi)
    paramsPerSignal$A[!aValidContion]<- -paramsPerSignal$A[!aValidContion]
  }
  
  return(paramsPerSignal)
}
