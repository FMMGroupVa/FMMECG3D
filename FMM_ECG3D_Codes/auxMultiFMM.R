
#### Dependencies ####
require("RColorBrewer")

#### Internal multiFMM functions ####
## MultiFMM, first step: optimize common parameters

step1FMM3D <- function(optBase, vDataMatrix, 
                       weights = rep(1, nrow(as.matrix(vDataMatrix)))) {
  
  pars <- optBase[["base"]] %*% vDataMatrix
  mobiusRegression <- apply(X = pars, MARGIN = 2, FUN = function(x){
    mobiusRegression <- x[1] + x[2]*optBase[["cost"]] + x[3]*optBase[["sint"]]
  }, simplify = TRUE)
  
  residuals <- (vDataMatrix - mobiusRegression)^2
  
  RSS <- sum(t(weights*t(residuals)))
  return(c(optBase[["alpha"]], optBase[["omega"]], RSS))
}

# Function for optimization: only depends on (alpha, omega).
step2_commonAlphaOmega <- function(initialParams, vDataMatrix, weights){
  
  initialAlpha<-as.numeric(initialParams[1])
  initialOmega<-as.numeric(initialParams[2])
  timePoints <- seqTimes(length(vDataMatrix[,1]))

  tStar <- 2*atan(initialOmega*tan((timePoints-initialAlpha)/2))
  
  # Calculates estimated component per lead
  fittedValues <- apply(vDataMatrix, 2, function(x, tStar){
    dM <- cbind(rep(1, length(x)), cos(tStar), sin(tStar))
    mDeltaGamma <- solve(t(dM)%*%dM)%*%t(dM)%*%x
    yFit <- mDeltaGamma[1] + mDeltaGamma[2]*cos(tStar) + mDeltaGamma[3]*sin(tStar)
    return(yFit)
  }, tStar = tStar)
  
  # Weighted residuals
  residuales <- (vDataMatrix - fittedValues)^2
  RSS <- sum(t(weights*t(residuales)))
  
  return(RSS)
}

optimizeAlphaOmega<-function(vDataMatrix, baseGrid, fittedWaves, currentComp,
                             errorWeights, omegaMax = 0.7){
  
  residualsMatrix<-as.matrix(vDataMatrix)
  for(signalIndex in 1:ncol(vDataMatrix)){
    vData<-vDataMatrix[,signalIndex]; vData<-vData[!is.na(vData)]; nObs<-length(vData)
    residualsMatrix[,signalIndex] <- vData - apply(as.matrix(fittedWaves[[signalIndex]][,-currentComp]), 1, sum)
  }
  # Grid step. RSS is a weighted mean of the RSS(i), i in cols(vDataMatrix)
  step1 <- lapply(FUN = step1FMM3D, X = baseGrid, vDataMatrix = residualsMatrix, 
                  weights = errorWeights)
  
  step1 <- matrix(unlist(step1), ncol=3, byrow=T)
  
  bestParamsIndex <- which.min(step1[,3])
  
  alpha <- step1[bestParamsIndex, 1]
  omega <- step1[bestParamsIndex, 2]
  
  # Post-optimization. Depends on (alpha, omega)
  nelderMead <- optim(par = c(alpha, omega), fn = step2_commonAlphaOmega,
                      vDataMatrix = residualsMatrix, method = "L-BFGS-B", 
                      lower = c(-2*pi, 0.0001), upper = c(4*pi, omegaMax), weights = errorWeights)

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
  
  fmmObject<-FMM:::FMM(
    M = mDeltaGamma[1],
    A = sqrt(delta^2+gamma^2),
    alpha = (alpha+2*pi)%%(2*pi),
    beta = (atan2(-gamma, delta) + alpha + 2*pi)%%(2*pi),
    omega = omega,
    timePoints = timePoints,
    summarizedData = vData,
    fittedValues = fittedFMMvalues,
    SSE = SSE,
    R2 = FMM:::PV(vData, fittedFMMvalues),
    nIter = 0
  )
  
  fmmObject@nPeriods<-1
  fmmObject@data<-vData
  
  return(fmmObject)
}

recalculateMA<-function(vDatai, paramsPerSignal){
  
  nObs<-length(vDatai)
  
  ## Extract Params. and calculate cosPhi
  alpha<-paramsPerSignal$Alpha; maxComp<-length(alpha)
  
  ## Have the waves been assigned?
  if("R" %in% rownames(paramsPerSignal)){
    assignedWaves<-TRUE
    useWave<-!substr(rownames(paramsPerSignal),1,1) %in% c("X","O")
  }else{
    assignedWaves<-FALSE; useWave<-rep(TRUE, maxComp)
  }
  
  alpha<-alpha[useWave]; beta<-paramsPerSignal$Beta[useWave]
  omega<-paramsPerSignal$Omega[useWave]
  cosPhi<-t(calculateCosPhi(alpha = alpha, beta = beta,
                            omega = omega, timePoints = seqTimes(nObs)))
  
  ## Recalculate M and As
  currentFormula<-as.formula(paste("vDatai~", paste("cosPhi[",1:sum(useWave),",]",
                                                    sep="",collapse="+")))
  regression<-lm(currentFormula)
  recalculatedM<-regression$coefficients[1]
  recalculatedAs<-as.vector(regression$coefficients[-1])
  
  ## Calculate R2, including other waves
  individualR2<-PVj(vData = vDatai, timePoints = seqTimes(nObs),
                    alpha = paramsPerSignal$Alpha,
                    beta = paramsPerSignal$Beta, omega = paramsPerSignal$Omega)
  paramsPerSignal$Var<-individualR2
  
  # Arrays for the saving of the results
  mResults<-rep(NA,maxComp); aResults<-rep(NA,maxComp)
  
  ## If there are discarded waves, the R2 are calculated with just the assigned waves
  if(assignedWaves & any(substr(rownames(paramsPerSignal),1,1)=="X")){
    paramsPerSignal$assignedR2<-rep(NA,nrow(paramsPerSignal))
    paramsPerSignal$assignedR2[useWave]<-PVj(vData = vDatai, timePoints = seqTimes(nObs),
                                             alpha = alpha, beta = beta, omega = omega)
    
    ## The A of the X waves is calculated with the other waves included
    xIndex<-which(substr(rownames(paramsPerSignal),1,1)=="X")
    xParams<-paramsPerSignal[xIndex,]
    cosPhi<-rbind(cosPhi,t(calculateCosPhi(alpha = xParams$Alpha, beta = xParams$Beta,
                                           omega = xParams$Omega, timePoints = seqTimes(nObs))))
    currentFormula<-as.formula(paste("vDatai~", paste("cosPhi[",1:(sum(useWave)+1),",]",
                                                      sep="",collapse="+")))
    aResults[xIndex]<-lm(currentFormula)$coefficients[sum(useWave)+2]
    paramsPerSignal$Var<-paramsPerSignal$assignedR2; paramsPerSignal$assignedR2<-NULL
  }
  
  # Results are saved
  mResults[useWave]<-rep(recalculatedM, sum(useWave)); aResults[useWave]<-recalculatedAs
  paramsPerSignal$M<-mResults; paramsPerSignal$A<-aResults
  
  ## If there is any negative A, they become 0.001
  # aValidContion<-aResults>0 | is.na(aResults)
  # if(any(!aValidContion)) aResults[!aValidContion]<-0.001
  #
  
  ## Negative A => Translation on the beta parameter
  aValidContion<-aResults>0 | is.na(aResults)
  if(any(!aValidContion)){
    # Change beta parameters
    paramsPerSignal$Beta[!aValidContion]<-(paramsPerSignal$Beta[!aValidContion]+pi)%%(2*pi)
    
    # Change A parameters
    paramsPerSignal$A[!aValidContion]<- -paramsPerSignal$A[!aValidContion]
  }
  
  return(paramsPerSignal)
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
  # print(fittedWaves)
  
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
