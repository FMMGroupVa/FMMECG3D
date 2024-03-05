
#### Dependencies ####
require("RColorBrewer")
require(R.utils)
sourceDirectory("funcionesFMM/")
# source("auxMultiFMMPlot.R") # incluidas al final

#### Internal multiFMM functions ####



step2_commonAlphaOmega <- function(initialParams, vDataMatrix, omegaMax){
  
  initialAlpha<-initialParams[1]; initialOmega<-initialParams[2]
  nSignals<-ncol(vDataMatrix); residualSS<-0
  rest1<-TRUE; rest2<-TRUE; rest3<-TRUE; rest4<-TRUE
  
  for(i in 1:nSignals){
    parameters<-initialParams[3:5+((i-1)*3)]
    vData<-vDataMatrix[,i]; vData<-vData[!is.na(vData)]; nObs<-length(vData)
    timePoints<-seqTimes(nObs)
    
    # FMM model and residual sum of squares
    modelFMM <- parameters[1] + parameters[2] *
      cos(parameters[3]+2*atan2(initialOmega*sin((timePoints - initialAlpha)/2),
                                cos((timePoints - initialAlpha)/2)))
    residualSS <- residualSS+sum((modelFMM - vData)^2)/nObs
    sigma <- sqrt(residualSS*nObs/(nObs - 5))
    
    # When amplitude condition is valid, it returns RSS
    # else it returns infinite.
    # amplitudeUpperBound <- parameters[1] + parameters[2]
    # amplitudeLowerBound <- parameters[1] - parameters[2]
    # rest1 <- amplitudeUpperBound <= max(vData) + 1.96*sigma
    # if(!rest1){
    #   warning("step2FMM_fixedAlphaOmega - rest1 not met: amplitudeUpperBound > max(x(t))")
    #   rest1<-!rest1
    # }
    # rest2 <- amplitudeLowerBound >= min(vData) - 1.96*sigma
    # if(!rest2){
    #   warning("step2FMM_fixedAlphaOmega - rest2 not met: amplitudeLowerBound < max(x(t))")
    #   rest2<-!rest2
    # }
    
    # Other integrity conditions that must be met
    if(parameters[2] < 0) rest3 <- FALSE  # A > 0
    if(initialOmega < 0  |  initialOmega > omegaMax) rest4<-FALSE # omega in (0, omegaMax]
  }
  
  #if(rest1 & rest2 & rest3 & rest4)  # Disabled rest1 and rest2
  if(rest3 & rest4)  # Disabled rest1 and rest2
    return(residualSS)
  else
    return(Inf)
}

optimizeAlphaOmega<-function(vDataMatrix, baseGrid, fittedWaves, currentComp,
                             alphaGrid, omegaGrid, errorWeights, usedApply,
                             omegaMax = 0.7){
  
  nSignals<-ncol(vDataMatrix)
  step1OutputNames <- c("alpha","omega","RSS")
  
  grid <- expand.grid(alphaGrid, omegaGrid)
  
  rssMatrix<-matrix(NA, nrow = nrow(grid), ncol = (nSignals)+2,
                    dimnames = list(NULL,c("alpha", "omega", paste("RSS", 1:nSignals, sep = ""))))
  rssMatrix[,1]<-grid[,1]; rssMatrix[,2]<-grid[,2]
  
  #### Initial estimation based on step1 ####
  residualsMatrix<-vDataMatrix
  
  for(signalIndex in 1:nSignals){
    vData<-vDataMatrix[,signalIndex]; vData<-vData[!is.na(vData)]; nObs<-length(vData)
    residualsMatrix[,signalIndex] <- vData - apply(as.matrix(fittedWaves[[signalIndex]][,-currentComp]), 1, sum)
    
    step1 <- lapply(FUN = step1FMM, X = baseGrid, vData = residualsMatrix[,signalIndex])
    
    step1 <- matrix(unlist(step1), ncol=3, byrow=T)
    colnames(step1) <- step1OutputNames
    rssMatrix[,signalIndex+2]<-errorWeights[signalIndex]*step1[,"RSS"]
  }
  
  # Solutions with big omegas are not permitted
  notValidIndex<-rssMatrix[,"omega"]>omegaMax
  if(all(notValidIndex)){
    notValidIndex<-rep(FALSE,nrow(grid))
  }
  rssMatrix[notValidIndex,"RSS1"]<-Inf

  # The best alpha and omega is the index that fulfills:
  bestParamsIndex <- which.min(apply(rssMatrix[,-c(1:2)],1,mean))

  # OPTIMIZACION ANTIGUA
  # 1ro: calcular parametros M, A, beta para todas las leads 
  
  alpha <- rssMatrix[bestParamsIndex, 1]
  omega <- rssMatrix[bestParamsIndex, 2]
  
  fullParameters <- rep(0,2+3*nSignals)
  fullParameters[1] <- alpha
  fullParameters[2] <- omega
  timePoints <- seqTimes(nObs=nrow(vDataMatrix))
  tStar <- 2*atan2(omega*sin((timePoints-alpha)/2), cos((timePoints-alpha)/2))
  dM <- cbind(1+0*tStar, cos(tStar), sin(tStar))
  dMP <- solve(t(dM)%*%dM)%*%t(dM)

  for(signalIndex in 1:nSignals){
    mDeltaGamma <- dMP %*% residualsMatrix[,signalIndex]
    delta <- mDeltaGamma[2]
    gamma <- mDeltaGamma[3]
    fullParameters[3+((signalIndex-1)*3)] <- c(mDeltaGamma[1])
    fullParameters[4+((signalIndex-1)*3)] <- c(sqrt(delta^2+gamma^2))
    fullParameters[5+((signalIndex-1)*3)] <- (c(atan2(-gamma, delta)) + 2*pi) %% (2*pi)
  }
  nelderMead <- optim(par = round(fullParameters, 10), fn = step2_commonAlphaOmega,
                      vDataMatrix = residualsMatrix, omegaMax = 0.7)
  return(nelderMead$par[1:2])
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

fitMultiFMM<-function(vDataMatrix, commonOmega=TRUE, nBack=5, maxIter=10,
                      weightError=TRUE, lengthAlphaGrid = 48, lengthOmegaGrid = 24,
                      alphaGrid = seq(0, 2*pi, length.out = lengthAlphaGrid),
                      omegaMin = 0.0001, omegaMax = 1,
                      omegaGrid = exp(seq(log(omegaMin), log(omegaMax), length.out = lengthOmegaGrid)),
                      numReps = 3, parallelize = TRUE, plotWithPlotECG=FALSE,
                      plotJustLast=FALSE, filename=NA){

  nSignals<-ncol(vDataMatrix)

  # Preprocessing
  #preprocessingWave<-list()

  # Is there trend? Subtract it if its the case. Furthermore, normalize signal
  #detrendPreprocesses<-multiDetrend(vDataMatrix, filename=NA)
  #fittedTrends<-sapply(1:nSignals, function(x) detrendPreprocesses[,x]$fittedTrend)
  #fittedTrends<-sapply(1:nSignals, function(x) detrendPreprocesses[,x])

  ## Is there trend? Subtract it if its the case. Furthermore, normalize signal
  #for(i in 1:nSignals) {
  #  preprocessingWave[[i]]<-list()
  #  preprocessingWave[[i]]$trend<-fittedTrends[,i]
  #  preprocessingWave[[i]]$scaling<- minimax(vDataMatrix[,i]-preprocessingWave[[i]]$trend)
  #  vDataMatrix[,i] <- preprocessingWave[[i]]$scaling$minimax
  #}

  #### Preparations for the fit ####

  # Preparation of the paralellization
  if(typeof(parallelize)=="logical"){
    usedApply_Cluster <- getApply(parallelize = parallelize)
    usedApply <- usedApply_Cluster[[1]]

    ## The paralellized cluster can also be previously instantiated and used
  }else{
    usedApply <- parallelize
  }

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
      #### First Step: determine optimal common parameters (just alpha or alpha and omega) ####
      optimalParams <- optimizeCommonParameters(commonOmega = commonOmega, vDataMatrix = vDataMatrix,
                                                fittedWaves = fittedWaves, currentComp = j,
                                                alphaGrid = alphaGrid, omegaGrid = omegaGrid, omegaMax = omegaMax,
                                                errorWeights = errorWeights, usedApply = usedApply)

      #### Second Step: fit single FMM wave in each signal with common parameters ####
      for(signalIndex in 1:nSignals){
        vData<-vDataMatrix[,signalIndex]; vData<-vData[!is.na(vData)]; nObs<-length(vData)
        vData<-vData - apply(as.matrix(fittedWaves[[signalIndex]][,-j]), 1, sum)
        fmmObjectList[[signalIndex]]<-optimizeOtherParameters(commonOmega = commonOmega,
                                                              vData = vData, fixedParams = optimalParams,
                                                              usedApply = usedApply)
        fittedWaves[[signalIndex]][,j]<-getFittedValues(fmmObjectList[[signalIndex]])
      }

      #### Get FMM parameters per wave ####
      params<-data.frame("M"=sapply(fmmObjectList,getM),"A"=sapply(fmmObjectList,getA),"Alpha"=sapply(fmmObjectList,getAlpha),
                         "Beta"=sapply(fmmObjectList,getBeta),"Omega"=sapply(fmmObjectList,getOmega), "Var"=rep(NA,nSignals))
      paramsPerWave<-lapply(1:nSignals, function(x) rbind(paramsPerWave[[x]], params[x,]))
      ## Recalculate As and Ms
      paramsPerWave<-lapply(1:nSignals, function(x) recalculateMA(vDatai=vDataMatrix[,x], currentBackResults=paramsPerWave[[x]]))

      #### Error is weighted across signals ####
      sigma<-sapply(1:nSignals, function(x){sum((vDataMatrix[,x]-apply(as.matrix(fittedWaves[[x]]),
                                                                       1, sum))^2)/(nObs-1)})
      if(weightError) errorWeights<-1/sigma
      else errorWeights<-rep(1,nSignals)
    }
    totalMSE[currentBack]<-sum(sigma)

    #### Plot ####
    if(!plotJustLast){
      if(!plotWithPlotECG){
        plotMultiFMM(vDatai=vDataMatrix, fittedWaves = fittedWaves, currentBack = currentBack,
                     leadNames=colnames(vDataMatrix), currentBackResults = paramsPerWave,
                     plotToFile = !is.na(filename), unassigned=TRUE, filename = filename)
      }else{
        plotMultiFMM_ECG(vDataMatrix=vDataMatrix, fittedWaves = fittedWaves, leadNames = colnames(vDataMatrix), #leadNames = gsub("^.{0,4}", "", colnames(vDataMatrix)),
                         currentBack = currentBack, paramsPerLead = paramsPerWave, unassigned = TRUE, filename = filename, plotToFile = TRUE)
      }
    }

    # Check if backfitting should stop
    stopCondition1<-currentBack>=maxIter
    stopCondition2<-ifelse(currentBack>1,totalMSE[currentBack]>totalMSE[currentBack-1], FALSE)
    if(stopCondition1 | stopCondition2){
      stopCriteria<-ifelse(stopCondition1, "Maximum Iterations",
                           "Minimum MSE")
      cat(paste("  Stop: ",stopCriteria,"\n",sep=""))
      continueBackfitting<-FALSE

      #### Plot, if plotJustLast ####
      if(plotJustLast){
        if(!plotWithPlotECG){
          plotMultiFMM(vDatai=vDataMatrix, fittedWaves = fittedWaves, currentBack = currentBack,
                       leadNames=colnames(vDataMatrix), currentBackResults = paramsPerWave,
                       plotToFile = !is.na(filename), unassigned=TRUE, filename = filename)
        }else{
          plotMultiFMM_ECG(vDataMatrix=vDataMatrix, fittedWaves = fittedWaves, leadNames = colnames(vDataMatrix), #leadNames = gsub("^.{0,4}", "", colnames(vDataMatrix)),
                           currentBack = currentBack, paramsPerLead = paramsPerWave, unassigned = TRUE, filename = filename, plotToFile = TRUE)
        }
      }

    }else{
      currentBack<-currentBack+1
    }
  }

  # Unname waves
  for(i in 1:nSignals) rownames(paramsPerWave[[i]])<-1:nBack

  # If the cluster has been instantiated in the function, it must be stopped
  if(typeof(parallelize)=="logical"){
    cluster <- usedApply_Cluster[[2]]
    if(!is.null(cluster)) parallel::stopCluster(cluster)
  }

  #### Return results ####
  return(paramsPerWave)
}

minimax<-function(y, minY=NULL, maxY=NULL){
  if(is.null(minY) | is.null(maxY))  minY<-min(y, na.rm = TRUE); maxY<-max(y, na.rm = TRUE)
  minimaxList<-list("minY"=minY, "maxY"=maxY,
                    "minimax"=-1+2*(y-minY)/(maxY-minY))
  return(minimaxList)
}

revMinimax<-function(values, minY, maxY){
  return(((maxY-minY)*(values+1)/2)+minY)
}

#### Other multiFMM useful functions ####

getFMMPeaksParams <- function(currentBackResults, nObs) {
  
  M<-currentBackResults["M"][1,1]; A<-currentBackResults["A"]
  assignedIndex<-!is.na(A)
  A<-A[assignedIndex]; alpha<-currentBackResults["Alpha"][assignedIndex]
  beta<-currentBackResults["Beta"][assignedIndex]; omega<-currentBackResults["Omega"][assignedIndex]
  nComp <- length(alpha)
  
  # Fiducial points estimation
  peakU <- rep(NA,nComp); peakL <- rep(NA,nComp); tpeakU <- rep(NA,nComp)
  tpeakL <- rep(NA,nComp); ZU <- rep(NA,nComp); ZL <- rep(NA,nComp)
  
  peakU[assignedIndex]<-(alpha+2*atan2(1/omega*sin(-beta/2),cos(-beta/2)))%%(2*pi)
  peakL[assignedIndex]<-(alpha+2*atan2(1/omega*sin((pi-beta)/2),cos((pi-beta)/2)))%%(2*pi)
  
  tpeakU[assignedIndex]<-peakU[assignedIndex]*nObs/(2*pi) + 1
  tpeakL[assignedIndex]<-peakL[assignedIndex]*nObs/(2*pi) + 1
  
  # Signal estimation
  phU <- lapply(1:nComp,function(k) (peakU[k]-alpha)/2)
  phL <- lapply(1:nComp,function(k) (peakL[k]-alpha)/2)
  ZU[assignedIndex] <- sapply(1:nComp, function(k) M+sum(A*cos(beta + 2*atan(omega*tan(phU[[k]]))), na.rm = TRUE))
  ZL[assignedIndex] <- sapply(1:nComp, function(k) M+sum(A*cos(beta + 2*atan(omega*tan(phL[[k]]))), na.rm = TRUE))
  names(ZU)<-NULL; names(ZL)<-NULL
  
  fiducialPoints<-cbind(peakU, peakL, tpeakU, tpeakL, ZU, ZL)
  colnames(fiducialPoints)<-c("PeakU", "PeakL", "tPeakU", "tPeakL", "ZU", "ZL")
  rownames(fiducialPoints)<-rownames(currentBackResults)
  
  return(as.data.frame(fiducialPoints))
}


#### Plot multiFMM functions ####
plotMultiFMM<-function(vDatai, fittedWaves, currentBack, currentBackResults,
                       filename=NA, leadNames=1:length(currentBackResults), path="./",
                       plotToFile=TRUE, unassigned=FALSE, extra=FALSE){
  
  nSignals<-length(currentBackResults)
  if(plotToFile){
    if(!is.na(filename)){
      png(filename=paste(path,filename,"_Back_",currentBack,".png",sep=""), type = "cairo-png",
          width=ifelse(nSignals<=12 | nSignals%%3!=0,
                       300*nSignals,50*nSignals),
          height=ifelse(nSignals<=12, 600,
                        ifelse(nSignals%%3!=0, 600*round(nSignals/12), 700*round(nSignals/12))))
    }else{
      stop("Filename must be provided if plotToFile=TRUE")
    }
  }
  if(nSignals>12 & nSignals%%3==0){
    sixthSignals<-nSignals/3; sumPlotOrder<-1:nSignals; comPlotOrder<-(nSignals+1):(2*nSignals)
    plotOrder<-c(rbind(sumPlotOrder, comPlotOrder))
    plotLayout <- matrix(plotOrder, nrow = sixthSignals, ncol = (2*nSignals)/sixthSignals)
    plotLayout <- rbind(plotLayout, rep((2*nSignals)+1, (2*nSignals)/sixthSignals))
    layout(mat = plotLayout, heights = c(rep(0.95/sixthSignals,sixthSignals),0.05))
    par(mar=c(0,0,0,0))
  }else if(nSignals>4){
    plotOrder<-1:(2*nSignals); halfSignals<-ceiling(nSignals/2)
    while(length(plotOrder)%%4!=0 | length(plotOrder)%%halfSignals!=0){
      plotOrder<-c(plotOrder, tail(plotOrder, 1))
    }
    plotLayout <- matrix(plotOrder, nrow = 4, ncol = halfSignals, byrow = TRUE)
    plotLayout <- rbind(plotLayout, rep((2*nSignals)+1, halfSignals))
    layout(mat = plotLayout, heights = c(0.1,0.1,0.1,0.1,0.1))
    par(mar=c(0.1,0.1,0.1,0.1))
  }else{
    plotOrder<-1:(2*nSignals)
    plotLayout <- matrix(plotOrder, nrow = 2, ncol = nSignals, byrow = TRUE)
    plotLayout <- rbind(plotLayout, rep((2*nSignals)+1, nSignals))
    layout(mat = plotLayout,heights = c(0.4,0.4,0.2))
    par(mar=c(0.1,0.1,0.1,0.1))
  }
  
  limitsSumPlot<-NA
  sapply(1:nSignals,function(x)
    plotMultiFMM_Sum(vDatai=vDatai[!is.na(vDatai[x,]),x], fittedWaves = fittedWaves[[x]],
                     currentBackResults = currentBackResults[[x]], currentBack=currentBack,
                     filename=filename, leadName=leadNames[x], unassigned=unassigned,
                     extra=extra, yLimits = limitsSumPlot))
  
  sapply(1:nSignals,function(y)
    plotMultiFMM_Comps(vDatai=vDatai[!is.na(vDatai[y,]),y], fittedWaves = fittedWaves[[y]],
                       currentBackResults = currentBackResults[[y]],
                       currentBack=currentBack, filename=filename, leadName=leadNames[y],
                       plotLegend=FALSE, unassigned=unassigned))
  
  # Add legend to the plot
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  if(unassigned){
    usedColors<-c(brewer.pal(n = 9, name = "Set1"), "aquamarine", "saddlebrown", "darkolivegreen1")
    usedColors<-usedColors[1:nrow(currentBackResults[[1]])]
    names(usedColors)<-1:nrow(currentBackResults[[1]])
  }else{
    usedColors<-getUsedColors(currentBackResults[[1]])
  }
  
  legend(x = "top",inset = 0, legend = names(usedColors), col=usedColors,
         lwd=5, cex=1.2, horiz = TRUE)
  
  if(plotToFile) dev.off()
}

plotMultiFMM_Sum<-function(vDatai, fittedWaves, currentBackResults, currentBack, leadName,
                           filename=NA, path="./", plotToFile=FALSE, unassigned=FALSE,
                           extra=FALSE, yLimits=NA){
  par(mar = c(2,2,3,1))
  
  if(!is.na(filename)){
    ecgId_beatId<-strsplit(filename, "_")[[1]]
    ecgId<-substr(ecgId_beatId[1], 3, nchar(ecgId_beatId[1]))
    beatId<-ecgId_beatId[2]
  }
  
  nObs<-length(vDatai)
  if(unassigned){
    assignedCondition<-rep(TRUE,nrow(currentBackResults))
  }else{
    assignedCondition<-!substr(rownames(currentBackResults),1,1) %in% c("X","O")
  }
  
  assignedResults<-currentBackResults[assignedCondition,]
  assignedWavesSum<-generateFMM(M = assignedResults$M[1], A = assignedResults$A,
                                alpha = assignedResults$Alpha, beta = assignedResults$Beta,
                                omega = assignedResults$Omega, length.out = nObs, plot = F)$y
  
  ## Colors definition: color linked to assigned wave
  totalR2<-sum(assignedResults$Var)
  
  if(plotToFile){
    png(filename=paste(path,"/02 Results/Plots/",filename,"_Back_",currentBack,".png",sep=""),
        width=1200, height=900)
  }
  
  ## Fit plot
  # if(!is.na(filename)){
  #   mainText<-ifelse(extra,
  #                    paste("Lead ",leadName,": Patient ",ecgId,", beat ",beatId,", extraback. ",currentBack,sep=""),
  #                    paste("Lead ",leadName,": Patient ",ecgId,", beat ",beatId,", back. ",currentBack,sep=""))
  # }else{
  mainText<-paste("Signal ",leadName,", back. ",currentBack,sep="")
  #}
  
  assignedWavesSumPlot<-assignedWavesSum
  
  if(any(is.na(yLimits))){
    yLimits<-c(min(assignedWavesSumPlot,vDatai),
               max(assignedWavesSumPlot,vDatai))
  }
  
  plot(1:nObs,vDatai,type="l",ylim=yLimits,cex.main=1.6,
       main=mainText)
  lines(1:nObs,assignedWavesSumPlot,col=4,lwd=2)
  legend("bottomright",legend=sprintf(totalR2*100, fmt = "R2=%#.1f %%"), col=NA, lty=1,
         lwd=NA, pch=NA, text.width=40, x.intersp=0.25, bty = "n", cex=1.5,
         inset=c(0.2,0))
  
  if(plotToFile) dev.off()
}

plotMultiFMM_Comps<-function(vDatai, fittedWaves, currentBackResults,
                             currentBack, leadName, filename=NA, path="./",
                             plotLegend=TRUE, plotToFile=FALSE, unassigned=FALSE,
                             yLimits=NA){
  par(mar = c(2,2,1,1))
  
  maxComp<-nrow(currentBackResults)
  if(unassigned){
    assignedCondition<-rep(TRUE,nrow(currentBackResults))
    
    ## Colors definition: color linked to assigned wave
    assignedWavesColors<-c(brewer.pal(n = 9, name = "Set1"), "aquamarine", "saddlebrown", "darkolivegreen1")
    assignedWavesColors<-assignedWavesColors[1:nrow(currentBackResults)]
  }else{
    assignedCondition<-!substr(rownames(currentBackResults),1,1) %in% c("X","O")
    
    ## Colors definition: color linked to assigned wave
    assignedWavesColors<-brewer.pal(n = 5, name = "Set1")
    names(assignedWavesColors)<-c("R","T","S","P","Q")
  }
  otherWavesColors<-rev(brewer.pal(n = 7, name = "Set2"))
  
  if(plotToFile){
    png(filename=paste(path,"/02 Results/Plots/",filename,"_Back_",currentBack,".png",sep=""),
        width=1200, height=900)
  }
  
  ## Components plot
  if(any(is.na(yLimits))){
    yLimits<-c(min(sapply(1:maxComp, function(x) fittedWaves[,x]-median(fittedWaves[,x]))),
               max(sapply(1:maxComp, function(x) fittedWaves[,x]-median(fittedWaves[,x]))))
  }
  
  currentColor<-ifelse(assignedCondition[1],
                       ifelse(unassigned,assignedWavesColors[1],
                              assignedWavesColors[names(assignedWavesColors)==rownames(currentBackResults)[1]]),
                       otherWavesColors[1])
  
  plot(fittedWaves[,1]-median(fittedWaves[,1]), lty=ifelse(assignedCondition[1],1,2),
       col=currentColor, type="l", lwd=3, ylim=yLimits)
  if(assignedCondition[1]) usedColors<-c(currentColor)
  
  if(maxComp>1){
    notAssignedColorIndex<-ifelse(!assignedCondition[1],2,1)
    assignedColorIndex<-ifelse(assignedCondition[1],2,1)
    for(i in 2:maxComp){
      currentColor<-ifelse(assignedCondition[i],
                           ifelse(unassigned,assignedWavesColors[assignedColorIndex],
                                  assignedWavesColors[names(assignedWavesColors)==rownames(currentBackResults)[i]]),
                           otherWavesColors[notAssignedColorIndex])
      if(assignedCondition[i]){
        if(!exists("usedColors")){usedColors<-c(currentColor)
        }else{usedColors<-c(usedColors, currentColor)}
        assignedColorIndex<-assignedColorIndex+1
      }else{
        notAssignedColorIndex<-notAssignedColorIndex+1
      }
      lines(fittedWaves[,i]-median(fittedWaves[,i]),col=currentColor,
            lwd=3, lty=ifelse(assignedCondition[i],1,2))
    }
  }
  
  ## Legend of the components plot
  if(plotLegend){
    legend("bottomright",legend=paste("Wave",1:maxComp),col=usedColors,lty=rep(1,maxComp),
           lwd=rep(2,maxComp),pch=rep(NA,maxComp), cex = 1, text.width=25, x.intersp=0.25,
           y.intersp = 0.1, bty = "n")
  }
  
  if(plotToFile) dev.off()
}

getUsedColors<-function(currentBackResults){
  maxComp<-nrow(currentBackResults)
  assignedCondition<-!substr(rownames(currentBackResults),1,1) %in% c("X","O")
  
  ## Colors definition: color linked to assigned wave
  assignedWavesColors<-brewer.pal(n = 5, name = "Set1")
  names(assignedWavesColors)<-c("R","T","S","P","Q")
  otherWavesColors<-rev(brewer.pal(n = 7, name = "Set2"))
  
  currentColor<-ifelse(assignedCondition[1],
                       assignedWavesColors[names(assignedWavesColors)==rownames(currentBackResults)[1]],
                       otherWavesColors[1])
  usedColors<-c(currentColor)
  
  notAssignedColorIndex<-ifelse(!assignedCondition[1],2,1)
  assignedColorIndex<-ifelse(assignedCondition[1],2,1)
  for(i in 2:maxComp){
    currentColor<-ifelse(assignedCondition[i],
                         assignedWavesColors[names(assignedWavesColors)==rownames(currentBackResults)[i]],
                         otherWavesColors[notAssignedColorIndex])
    if(assignedCondition[i]){
      if(!exists("usedColors")){usedColors<-c(currentColor)
      }else{usedColors<-c(usedColors, currentColor)}
      assignedColorIndex<-assignedColorIndex+1
    }else{
      usedColors<-c(usedColors, currentColor)
      notAssignedColorIndex<-notAssignedColorIndex+1
    }
  }
  names(usedColors)<-rownames(currentBackResults)
  
  usedColors<-usedColors[order(names(usedColors))]
  
  return(usedColors)
}
