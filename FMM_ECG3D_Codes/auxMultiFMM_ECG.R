#### Dependencies ####
source("auxMultiFMM.R")
source("thresholdsMultiFMM_ECG.R")

#### Fit of multiFMM_ECG model functions ####
fitMultiFMM_ECG<-function(vDataMatrix, annotation,
                          commonOmega=TRUE, weightError=TRUE, showProgress=TRUE, parallelize = FALSE,
                          nBack=5, extraBack=2, extraBackP=3, maxIter=10, numReps = 3, excellentR2=0.98,
                          lengthAlphaGrid = 48, lengthOmegaGrid = 24, alphaGrid = seq(0, 2*pi, length.out = lengthAlphaGrid),
                          omegaMin = 0.0001, omegaMax = 1, omegaGrid = exp(seq(log(omegaMin), log(omegaMax), length.out = lengthOmegaGrid)),

 				  ecgId=NA, beatId=NA, patientClass=NA, patientSuperclass=NA,
                          resultsFile=NA, plotToFile=FALSE, saveResults=FALSE,
                          plotBasename=NA, plotPath="./03 Results", moveBackplotsFolder=TRUE){


dos<-FALSE;tres<-FALSE;cuatro<-FALSE;cinco<-FALSE
  ## Define leads sets and extract data
  eightLeads<-c("I", "II", "V1", "V2", "V3", "V4", "V5", "V6")
  vDataMatrix<-vDataMatrix[,eightLeads] # In the model fitting process only eight leads are used

  ## Check data integrity
  checkedData<-checkData(vDataMatrix)
  vDataMatrix<-checkedData$vDataMatrix; denseDataMatrix<-checkedData$denseDataMatrix

  ## Get measures related to available data
  leadNames<-colnames(denseDataMatrix); nObs<-nrow(denseDataMatrix); nSignals<-ncol(denseDataMatrix)
  availableRelevantLeads<-relevantLeads[relevantLeads %in% colnames(denseDataMatrix)]
  alphaN<-getAlphaN(annoActB=annotation, nObs=nObs)

  #### Preparations for the fit ####
  if(!is.na(ecgId) & !is.na(beatId)){
    singleFilename<-paste("HR",ecgId,"_", beatId,collapse = "", sep="")
  }else{
    singleFilename<-plotBasename
  }

  # Preparation of the paralellization
  if(typeof(parallelize)=="logical"){
    usedApply_Cluster <- getApply(parallelize = parallelize)
    usedApply <- usedApply_Cluster[[1]]

  ## The paralellized cluster can also be previously instantiated and used
  }else{
    usedApply <- parallelize
  }

  # The results are stored in:
  paramsPerLead <- replicate(nSignals, simplify = FALSE,
                             data.frame(M=double(), A=double(), Alpha=double(), Beta=double(), Omega=double()))
  fittedWaves<-lapply(1:nSignals,function(x){
    vData<-denseDataMatrix[,x]; vData<-vData[!is.na(vData)];
    return(array(0,c(length(vData),nBack+extraBack+extraBackP)))})
  paramsPerLeadPerBack<-list(); rChanges<<-""; uChanges<<-"" # Keep track of R and U changes

  # The MSE must be balanced in the fit of multiple signals
  errorWeights<-rep(1,nSignals);  signalWeights<-rep(1,nSignals); names(signalWeights)<-leadNames
  totalMSE<-rep(Inf, maxIter)
  # Frontal and Horizontal Plane Leads should have different weights
  if(weightError){
    nFrontalLeads<-sum(grepl("I",leadNames,fixed = TRUE))
    nHorizontalLeads<-length(leadNames)-nFrontalLeads
    signalWeights<-ifelse(grepl("I",leadNames,fixed = TRUE), 1/nFrontalLeads, 1/nHorizontalLeads)
    names(signalWeights)<-leadNames; errorWeights<-signalWeights
  }

  # Backfitting algorithm: iteration
  currentBack<-1; continueBackfitting<-TRUE
  if(showProgress){
    if(!is.na(ecgId) & !is.na(beatId)){
      cat(paste("Analysing ",paste("PTBXL ",ecgId,"-", beatId,", leads ",paste(leadNames,collapse = ", "),
                                   collapse = "", sep=""),"\n   Backfitting: ",sep=""))
    }else{
      cat(paste("Analysing leads ",paste(leadNames,collapse = ", "),"\n   Backfitting: ",sep=""))
    }
  }
  
  grid <- expand.grid(alphaGrid, omegaGrid)
  time = seqTimes(n)
  
  mobiusBase = function(par, times){
    par <- as.numeric(par)
    nonlinearMob = 2*atan(par[2]*tan((times-par[1])/2))
    M <- cbind(time*0+1, cos(nonlinearMob), sin(nonlinearMob))
    return(list(base = solve(t(M)%*%M)%*%t(M), 
                alpha = par[1], omega = par[2], 
                cost = cos(nonlinearMob), 
                sint = sin(nonlinearMob))) 
  }
  
  optBase <- apply(grid, 1, FUN = mobiusBase, times = time, simplify = FALSE)
  
  while(continueBackfitting) {
    if(showProgress) cat(paste(currentBack," ",sep=""))
    fmmObjectList<-list()
    paramsPerLead <- replicate(nSignals, simplify = FALSE,
                               data.frame(M=double(), A=double(),
                                          Alpha=double(), Beta=double(), Omega=double()))

    # Backfitting algorithm: component
    for(j in 1:nBack){
      #### First Step: determine optimal common parameters (just alpha or alpha and omega) ####
      optimalParams <- optimizeAlphaOmega(vDataMatrix = denseDataMatrix, baseGrid = optBase, fittedWaves = fittedWaves, currentComp = j,
                                          omegaMax = omegaMax, errorWeights = errorWeights, usedApply = usedApply)

      #### Second Step: fit single FMM wave in each signal with common parameters ####
      for(signalIndex in 1:nSignals){
        vData<-denseDataMatrix[,signalIndex]; vData<-vData[!is.na(vData)];
        vData<-vData - apply(as.matrix(fittedWaves[[signalIndex]][,-j]), 1, sum)
        fmmObjectList[[signalIndex]]<-optimizeOtherParameters(vData = vData, fixedParams = optimalParams)
        fittedWaves[[signalIndex]][,j]<-getFittedValues(fmmObjectList[[signalIndex]])
      }

      #### Get FMM parameters per wave ####
      params<-data.frame("M"=sapply(fmmObjectList,getM),"A"=sapply(fmmObjectList,getA),"Alpha"=sapply(fmmObjectList,getAlpha),
                         "Beta"=sapply(fmmObjectList,getBeta),"Omega"=sapply(fmmObjectList,getOmega), "Var"=rep(NA,nSignals))
      paramsPerLead<-lapply(1:nSignals, function(x) rbind(paramsPerLead[[x]], params[x,], make.row.names=FALSE))

      ## Recalculate As and Ms
      paramsPerLead<-lapply(1:nSignals, function(x) recalculateMA(vDatai=denseDataMatrix[,x], paramsPerSignal=paramsPerLead[[x]]))

      ## Order by RelVar
      paramsPerLead<-setNames(paramsPerLead,leadNames)
      paramsPerLead_fittedWaves<-wavesInRelValOrder(vDataMatrix=denseDataMatrix, paramsPerLead=paramsPerLead, fittedWaves=fittedWaves)
      paramsPerLead<-paramsPerLead_fittedWaves$paramsPerLead; fittedWaves<-paramsPerLead_fittedWaves$fittedWaves

      #### Error is weighted across signals ####
      sigma<-getAssignedSigma(vDataMatrix=denseDataMatrix, fittedWaves=fittedWaves, paramsPerLead=paramsPerLead,
                              weightError=weightError, signalWeights=signalWeights, assigned=FALSE)

	}
  
    #### Wave assignment ####
    paramsPerLead<-setNames(paramsPerLead,leadNames)
    paramsPerLead<-waveAssignation(paramsPerLead=paramsPerLead, nObs = nObs, alphaN=alphaN)
  storeSinMA<- paramsPerLead
    ## Recalculate As and Ms
    paramsPerLead<-lapply(1:nSignals, function(x) recalculateMA(vDatai=denseDataMatrix[,x], paramsPerSignal=paramsPerLead[[x]]))
    # print(paramsPerLead)
    
    ## Order by RelVar
    paramsPerLead<-setNames(paramsPerLead,leadNames)
    paramsPerLead_fittedWaves<-wavesInRelValOrder(vDataMatrix=denseDataMatrix, paramsPerLead=paramsPerLead, fittedWaves=fittedWaves)
    paramsPerLead<-paramsPerLead_fittedWaves$paramsPerLead; fittedWaves<-paramsPerLead_fittedWaves$fittedWaves
  
    ## Final error weights and MSE are calculated with just assigned waves
    assignedWaves<-substr(rownames(paramsPerLead[[1]]),1,1)!="X"
    sigma<-getAssignedSigma(vDataMatrix=denseDataMatrix, fittedWaves=fittedWaves, paramsPerLead=paramsPerLead,
                            weightError=weightError, signalWeights=signalWeights, assigned=TRUE)
    totalMSE[currentBack]<-sum(sigma)

    # If there is a wave that is not assigned, search for more
    unusedWave<-!assignedWaves; additionalComponents<-NA; extraComponent<-NA; paramsPerLead<-setNames(paramsPerLead,leadNames)
    if(any(unusedWave)){
      additionalComponents<-searchAdditionalComponents(vDataMatrix=denseDataMatrix, baseGrid = optBase, fittedWaves=fittedWaves, alphaN=alphaN, leadNames=leadNames,
                                                       weightError=weightError, errorWeights=errorWeights, signalWeights=signalWeights,
                                                       extraBack=(nBack+1):(nBack+extraBack), paramsPerLead=paramsPerLead, usedApply=usedApply,
                                                       singleFilename=singleFilename, plotBacks=FALSE)
      # The extra components find the other waves
      if(sum(assignedWaves)<sum(substr(rownames(additionalComponents$paramsPerLead[[1]]),1,1)!="X")){
        paramsPerLead<-additionalComponents$paramsPerLead
  dos<-TRUE
  varX<-paramsPerLead

        fittedWaves<-additionalComponents$fittedWaves
        errorWeights<-additionalComponents$errorWeights
        extraComponent<-additionalComponents$extraComponent
        assignedWaves<-substr(rownames(paramsPerLead[[1]]),1,1)!="X"
        sigma<-getAssignedSigma(vDataMatrix=denseDataMatrix, fittedWaves=fittedWaves, paramsPerLead=paramsPerLead,
                                weightError=weightError, signalWeights=signalWeights, assigned = TRUE)
        totalMSE[currentBack]<-sum(sigma)
      }
    }

    # Check P: if P is doubtful, there may be other better P candidates in waves 6-7
    newP<-FALSE
    if(doubtfulP(paramsPerLead, alphaN=alphaN, inSetPlusPlusP=TRUE)){

      # If it has not been done before, additional components are searched
      # Also, only 6 components may have been searched
      if(any(is.na(additionalComponents)) | (!is.na(extraComponent) & extraComponent==6)){
        if(is.na(extraComponent)){
          extraBacksNeeded<-(nBack+1):(nBack+extraBack)
        }else{
          extraBacksNeeded<-(nBack+extraBack)
        }
        paramsPerLead<-setNames(paramsPerLead,leadNames)
        additionalComponents<-searchAdditionalComponents(vDataMatrix=denseDataMatrix, baseGrid = optBase, fittedWaves=fittedWaves, alphaN=alphaN, leadNames=leadNames,
                                                         weightError=weightError, errorWeights=errorWeights, signalWeights=signalWeights,
                                                         extraBack = extraBacksNeeded, paramsPerLead=paramsPerLead, usedApply=usedApply,
                                                         singleFilename=singleFilename, plotBacks=FALSE, searchAllExtraBack=TRUE)

        oldWavesNames<-rownames(paramsPerLead[[1]])
        oldWavesNames<-c(oldWavesNames, paste("X",c(extraBacksNeeded-7)+20,sep=""))
        paramsPerLead<-additionalComponents$paramsPerLead
tres<-TRUE
varX<-paramsPerLead

        for(i in 1:nSignals) rownames(paramsPerLead[[i]])<-oldWavesNames
        fittedWaves<-additionalComponents$fittedWaves
      }

      # Calculate variance explained without assigned P
      selectedWaves<-substr(rownames(paramsPerLead[[1]]),1,1)!="X" & !substr(rownames(paramsPerLead[[1]]),1,1)=="P"
      alreadyExplainedVarWithoutP<-apply(sapply(availableRelevantLeads, function(x) ifelse(selectedWaves, paramsPerLead[[x]]$Var,0)), 2, sum)

      # Is there another candidate P?
      oldIndexP<-which(rownames(paramsPerLead[[1]])=="P")
      setPlusP<-getSetPlusP(paramsPerLead=paramsPerLead, alphaN=alphaN, setPlusPlusP=TRUE, freeWaves = !selectedWaves)
      setPlusP[oldIndexP]<-FALSE

      if(any(setPlusP, na.rm = TRUE)){

        # Should J' discard P?
        discardJP<-sapply(1:length(setPlusP), function(x){
          ifelse(setPlusP[x], shouldDiscardPJ(paramsPerLead, jIndexInSetPlusPlusP=x,
                                              alreadyExplainedVarWithoutP=alreadyExplainedVarWithoutP), NA)})

        # If a new wave is apt, it substitutes the current P
        if(any(discardJP, na.rm=TRUE)){
          pNewIndex<-which(discardJP)[1]; newP<-TRUE
          oldIndexQ<-ifelse("Q" %in% rownames(paramsPerLead[[1]]), which(rownames(paramsPerLead[[1]])=="Q"), NA)

          # Should old P discard Q?
          discardQP<-shouldDiscardQP(paramsPerLead, pNewIndex=pNewIndex)

          if(discardQP){ # Old P becomes Q
            oldLabelQ<-"XQ2"; oldLabelP<-"Q"
          }else{ # Q is left as it is
            oldLabelQ<-"Q"; oldLabelP<-"XP"
          }

          # Relabel waves
          if(!is.na(oldIndexQ)){
            for(i in 1:nSignals) rownames(paramsPerLead[[i]])[c(oldIndexP,oldIndexQ,pNewIndex)]<-c(oldLabelP,oldLabelQ,"P")
          }else{
            for(i in 1:nSignals) rownames(paramsPerLead[[i]])[c(oldIndexP,pNewIndex)]<-c(oldLabelP,"P")
          }

          # Get fitted waves
          fittedWaves<-additionalComponents$fittedWaves; errorWeights<-additionalComponents$errorWeights

          # Recalculate MSE
          assignedWaves<-substr(rownames(paramsPerLead[[1]]),1,1)!="X"
          sigma<-getAssignedSigma(vDataMatrix=denseDataMatrix, fittedWaves=fittedWaves, paramsPerLead=paramsPerLead,
                                  weightError=weightError, signalWeights=signalWeights)
          totalMSE[currentBack]<-sum(sigma)
        }
      }
    }

    # P still not assigned?
    if(!"P" %in% rownames(paramsPerLead[[1]])){

      # P could be within the first five waves but be more distant than expected
      assignedWaves<-substr(rownames(paramsPerLead[[1]]),1,1)!="X"
      alreadyExplainedVar<-apply(sapply(availableRelevantLeads,function(x) ifelse(assignedWaves,
                                                                                 paramsPerLead[[x]]$Var,0)), 2, sum)
      lastChanceP_paramsPerLead<-lastChanceP(paramsPerLead, alreadyExplainedVar=alreadyExplainedVar,
                                             alphaN=alphaN)

      if("P" %in% rownames(lastChanceP_paramsPerLead[[1]])){
        paramsPerLead<-lastChanceP_paramsPerLead; assignedWaves<-substr(rownames(paramsPerLead[[1]]),1,1)!="X"
        sigma<-getAssignedSigma(vDataMatrix=denseDataMatrix, fittedWaves=fittedWaves, paramsPerLead=paramsPerLead,
                                weightError=weightError, signalWeights=signalWeights, assigned=TRUE)
        totalMSE[currentBack]<-sum(sigma)
      }else{
        previouslyFoundWaves<-rownames(paramsPerLead[[1]])[substr(rownames(paramsPerLead[[1]]),1,1)!="X"]
        # P may be in the waves 8-10
        paramsPerLead<-setNames(paramsPerLead,leadNames)
        firstBackToSearch<-ifelse(is.na(extraComponent),nBack+1,nBack+extraBack+1)
        additionalComponentsP<-searchAdditionalComponents(vDataMatrix=denseDataMatrix, baseGrid = optBase, fittedWaves=fittedWaves, alphaN=alphaN, leadNames=leadNames,
                                                          weightError=weightError, errorWeights=errorWeights, signalWeights=signalWeights,
                                                          extraBack = firstBackToSearch:(nBack+extraBack+extraBackP), paramsPerLead=paramsPerLead,
                                                          usedApply=usedApply, singleFilename=singleFilename, plotBacks=FALSE, searchAllExtraBack=FALSE)

        # This search is more strict than the habitual P search
        if("P" %in% rownames(additionalComponentsP$paramsPerLead[[1]])){
          pFarIndex<-which(rownames(additionalComponentsP$paramsPerLead[[1]])=="P")
          for(i in nSignals) rownames(additionalComponentsP$paramsPerLead[[1]])[pFarIndex]<-"XFarP"
        }

        secondLastChanceP_paramsPerLead<-lastChanceP(additionalComponentsP$paramsPerLead, alreadyExplainedVar=alreadyExplainedVar,
                                                     alphaN = alphaN)

        if("P" %in% rownames(secondLastChanceP_paramsPerLead[[1]]) &
           sum(assignedWaves)<sum(substr(rownames(secondLastChanceP_paramsPerLead[[1]]),1,1)!="X")){
          paramsPerLead<-secondLastChanceP_paramsPerLead
cuatro<-TRUE
varX<-paramsPerLead

          # Get fitted waves
          fittedWaves<-additionalComponentsP$fittedWaves; errorWeights<-additionalComponentsP$errorWeights

          # Recalculate MSE
          assignedWaves<-substr(rownames(paramsPerLead[[1]]),1,1)!="X"
          sigma<-getAssignedSigma(vDataMatrix=denseDataMatrix, fittedWaves=fittedWaves, paramsPerLead=paramsPerLead,
                                  weightError=weightError, signalWeights=signalWeights, assigned=TRUE)
          totalMSE[currentBack]<-sum(sigma)
        }
      }
    }
    
    ## Recalculate As and Ms with just assigned waves
    unassignedWaves<-substr(rownames(paramsPerLead[[1]]),1,1)=="X"
    if(any(substr(rownames(paramsPerLead[[1]]),1,1)=="X")){
      paramsPerLead<-lapply(1:nSignals, function(x) recalculateMA(vDatai=denseDataMatrix[,x],
                                                                  paramsPerSignal=paramsPerLead[[x]]))
    }
    if(showProgress) plotMultiFMM_ECG(vDataMatrix=denseDataMatrix, fittedWaves = fittedWaves, leadNames = leadNames,
                                      currentBack = currentBack, paramsPerLead=paramsPerLead,
                                      filename = singleFilename, plotToFile = plotToFile, path=plotPath)
    
    ## Save results adequately
    paramsPerLeadPerBack[[currentBack]]<-setNames(paramsPerLead,leadNames)
    paramsPerLead<-setNames(paramsPerLead,leadNames)
    if(currentBack==1){ lastParamsPerLead<-NA
    }else{ lastParamsPerLead<-paramsPerLeadPerBack[[currentBack-1]] }

    # Check if backfitting should stop
    nAssignedWaves<-sum(!unassignedWaves); last_nAssignedWaves<-sum(substr(rownames(lastParamsPerLead[[1]]),1,1)!="X")
    currentMeanR2<-mean(sapply(1:nSignals, function(x) sum(paramsPerLead[[x]][!unassignedWaves,]$Var)))

    stopCondition1<-currentBack>=maxIter
    stopCondition2<-nAssignedWaves==5 & currentMeanR2>excellentR2
    stopCondition3<-nAssignedWaves<5 & last_nAssignedWaves==5 # N? of assigned waves is reduced
    stopCondition4<-ifelse(currentBack>1,totalMSE[currentBack]>totalMSE[currentBack-1], FALSE) &
      !(isImprovedP(paramsPerLead, lastParamsPerLead, alphaN=alphaN) | newP)

    if(((currentBack==1 & maxIter==1)| currentBack>2) &
       (stopCondition1 | stopCondition2 | stopCondition3 | stopCondition4)){
      continueBackfitting<-FALSE
      stopCriteria<-ifelse(stopCondition1, "maxIter",
                           ifelse(stopCondition2, "excellentR2",
                                  ifelse(stopCondition3, "Waves Unassignment",
                                         "MSE increase + Unimproved waves")))
      nSpacesLeft<-ifelse(maxIter==10,
                          ifelse(currentBack!=maxIter, (maxIter-currentBack)*2+1, 0), 0)
      if(showProgress) cat(paste(paste(rep(" ", nSpacesLeft), collapse = ""),
                "  Stop: ",stopCriteria,"\n",sep=""))

      # The backfitting with the most assigned waves and the highest R2 is returned
      bestBackIndex<-getBestBack(paramsPerLeadPerBack=paramsPerLeadPerBack,
                                 alphaN=alphaN, totalMSE=totalMSE)

      unusedBackfittings<-NA
      # The best fit is returned
      if(currentBack!=bestBackIndex){
        paramsPerLead<-paramsPerLeadPerBack[[bestBackIndex]]
        unusedBackfittings<-(bestBackIndex+1):(currentBack)
        currentBack<-bestBackIndex

        # Unused backfitting plots are renamed
        if(plotToFile){
          renameUnusedBackplots(filename=singleFilename, totalBack=currentBack,
                                unusedBackfittings=unusedBackfittings)
        }
      }

    }else{
      unassignedWaves<-substr(rownames(paramsPerLead[[1]]),1,1)=="X"
      currentBack<-currentBack+1

      ## Assigned Waves are moved to the first five positions
      if(any(which(!unassignedWaves)>5)){
        nWavesToMove<-length(which(!unassignedWaves)>5)
        for(i in 1:nWavesToMove){
          unassignedIndex<-which(unassignedWaves)[1]
          assignedIndex<-tail(which(!unassignedWaves),1)

          # Exchange contents
          for(i in 1:nSignals){
            fittedWaves[[i]][,unassignedIndex]<-fittedWaves[[i]][,assignedIndex]
          }
          unassignedWaves[assignedIndex]<-TRUE
          unassignedWaves[unassignedIndex]<-FALSE
        }
      }
      fittedWaves<-lapply(1:nSignals, function(x){
        fittedWaves[[x]][,unassignedWaves]<-0
        return(fittedWaves[[x]])})
    }
  }

  # If the cluster has been instantiated in the function, it must be stopped
  if(typeof(parallelize)=="logical"){
    cluster <- usedApply_Cluster[[2]]
    if(!is.null(cluster)) parallel::stopCluster(cluster)
  }

  #### Save results ####
   if(plotToFile & moveBackplotsFolder){
     moveBackplots(filename=singleFilename, totalBack=currentBack, unusedBackfittings=unusedBackfittings)
   }

   if(saveResults){
     if(is.na(resultsFile)){
       resultsFile<-"multiLeadResults.csv"
     }
     saveResultsCSV(ecgId, beatId, leadNames, paramsPerLead, superclass=patientSuperclass,
                    class=patientClass, totalBack=currentBack, stopCriteria=stopCriteria, resultsFile=resultsFile,
                    annotation=annotation, alphaN=alphaN, nObs=nObs, rChanges=paste(unique(rChanges[rChanges!=""]), collapse = "; "),
                    uChanges=paste(unique(uChanges[uChanges!=""]), collapse = "; "))
   }
if(!dos & !tres & !cuatro)varX<-storeSinMA
  return(paramsPerLead)
}

# Used to fit a MultiFMM of m1 components and save just plots + R2
fitMultiFMM_justR2<-function(vDataMatrix, commonOmega=TRUE, weightError=TRUE, showProgress=TRUE, parallelize = FALSE,
                             nBack=12, maxIter=5, numReps = 3, excellentR2=0.98,
                             lengthAlphaGrid = 48, lengthOmegaGrid = 24, alphaGrid = seq(0, 2*pi, length.out = lengthAlphaGrid),
                             omegaMin = 0.0001, omegaMax = 1, omegaGrid = exp(seq(log(omegaMin), log(omegaMax), length.out = lengthOmegaGrid)),

                             ## Arguments to be suppressed in app version
                             ecgId=NA, beatId=NA, patientClass=NA, patientSuperclass=NA,
                             resultsFile=NA, plotToFile=FALSE, saveResults=FALSE,
                             plotBasename=NA, plotPath="./03 Results"){

  ## Define leads sets and extract data
  eightLeads<-c("I", "II", "V1", "V2", "V3", "V4", "V5", "V6")
  vDataMatrix<-vDataMatrix[,eightLeads] # In the model fitting process only eight leads are used

  ## Check data integrity
  checkedData<-checkData(vDataMatrix)
  vDataMatrix<-checkedData$vDataMatrix; denseDataMatrix<-checkedData$denseDataMatrix

  ## Get measures related to available data
  leadNames<-colnames(denseDataMatrix); nObs<-nrow(denseDataMatrix); nSignals<-ncol(denseDataMatrix)
  availableRelevantLeads<-relevantLeads[relevantLeads %in% colnames(denseDataMatrix)]

  #### Preparations for the fit ####
  singleFilename<-paste0(paste("HR",ecgId,"_", beatId,collapse = "", sep=""),"_m=",nBack)

  #### Fit models ####
  multiFmm<-fitMultiFMM(vDataMatrix=denseDataMatrix, commonOmega=commonOmega, nBack=nBack, maxIter=maxIter, numReps = numReps,
                         weightError=weightError, lengthAlphaGrid = lengthAlphaGrid, alphaGrid = alphaGrid,
                         lengthOmegaGrid = lengthOmegaGrid, omegaMin = omegaMin, omegaMax = omegaMax, omegaGrid = omegaGrid,
                         parallelize = parallelize, plotWithPlotECG=TRUE, filename=singleFilename, plotJustLast = TRUE)

  #### Save R2 per lead ####
  # Some leads may not have not been used
  if(nSignals<8){
    unusedLeads<-eightLeads[!eightLeads %in% leadNames]
    for(i in unusedLeads){
      nSignals<-nSignals+1
      multiFmm[[nSignals]]<-as.data.frame(matrix(NA, nrow = nrow(multiFmm[[1]]), ncol = ncol(multiFmm[[1]]),
                                                  dimnames = dimnames(multiFmm[[1]])))
      leadNames<-c(leadNames, i)
    }
  }
  names(multiFmm)<-leadNames
  multiFmm<-multiFmm[eightLeads]
  r2Fmm<-data.frame(t(sapply(1:nBack, function(x) sapply(multiFmm, function(y) sum(y$Var[1:x])))))
  rownames(r2Fmm)<-paste0("R2_",1:nBack)
  resultWide<-r2Fmm %>%rownames_to_column() %>%
    pivot_longer(cols = -rowname) %>%
    unite(name, rowname, name, sep = "_") %>%
    pivot_wider()
  resultWide$EcgId<-ecgId; resultWide$BeatId<-beatId
  resultWide<-resultWide[,c(tail(1:ncol(resultWide),2), head(1:ncol(resultWide),ncol(resultWide)-2))]

  #### Save PRD per lead ####
  PRD <- function(vData, pred) return(sqrt((sum((vData - pred)^2))/(sum(vData^2))))
  PRDj <- function(vData, alpha, beta, omega){
    nComponents <- length(alpha); timePoints = FMM::seqTimes(length(vData))
    waves <-  calculateCosPhi(alpha = alpha, beta = beta, omega = omega, timePoints = timePoints)

    # The percentage of variability explained up to wave i is determined
    cumulativePRD <- sapply(1:nComponents, function(x){PRD(vData, predict(lm(vData ~ waves[,1:x])))})

    # individual percentage of variability is the part that adds to the whole
    return(cumulativePRD)
  }

  #cumulativePRD<-PRDj(vData = vDataMatrix[,1], alpha = multiFmm[[1]]$Alpha, beta = multiFmm[[1]]$Beta, omega = multiFmm[[1]]$Omega)
  prdFmm<-data.frame(sapply(eightLeads, function(x){
    if(!x %in% colnames(denseDataMatrix)){
      return(rep(NA, nBack-4))
    }else{
      return(PRDj(vData = vDataMatrix[,x], alpha = multiFmm[[x]]$Alpha,
                  beta = multiFmm[[x]]$Beta, omega = multiFmm[[x]]$Omega)[5:nBack])}}))
  rownames(prdFmm)<-paste0("PRD_",5:nBack)
  resultWide<-cbind(resultWide, prdFmm %>%rownames_to_column() %>%
    pivot_longer(cols = -rowname) %>%
    unite(name, rowname, name, sep = "_") %>%
    pivot_wider())
  resultWide$errorMessage<-NA

  #### Save results ####
  if(is.na(resultsFile)){
    resultsFile<-"multiLeadResults.csv"
  }
  destFile<-paste("./03 Results", resultsFile, sep="/")
  if(file.exists(destFile)){
    write.table(resultWide, file = destFile,sep = ",",
                append = TRUE, quote = FALSE,col.names = FALSE,
                row.names = FALSE)
  } else{
    write.table(resultWide, file = destFile,sep = ",",
                append = TRUE, quote = FALSE,col.names = TRUE,
                row.names = FALSE)

  }
}

#### Functions related to the MultiFMM_ECG fit ####
wavesInRelValOrder<-function(vDataMatrix, paramsPerLead, fittedWaves, unusedComp=NA){
  availableRelevantLeads<-relevantLeads[relevantLeads %in% names(paramsPerLead)]
  nSignals<-length(paramsPerLead); nWaves<-nrow(paramsPerLead[[1]])

  # With just a single wave, they are correctly ordered
  if(nWaves==1){
    return(list("fittedWaves"=fittedWaves, "paramsPerLead"=paramsPerLead, "unusedComp"=unusedComp))
  }else{
    relevantVar<-apply(sapply(availableRelevantLeads,function(x) paramsPerLead[[x]]$Var), 1, mean)
    orderByRelVar<-order(relevantVar, decreasing = TRUE)
    currentOrder<-1:nWaves

    # The waves order must be changed and M and A must be recalculated
    if(!all(orderByRelVar==currentOrder)){

      # Reorder waves parameters
      leadNames<-names(paramsPerLead)
      paramsPerLead<-lapply(1:nSignals, function(x) paramsPerLead[[x]][orderByRelVar,])
      paramsPerLead<-lapply(1:nSignals, function(x) recalculateMA(vDatai=vDataMatrix[,x], paramsPerSignal=paramsPerLead[[x]]))
      paramsPerLead<-setNames(paramsPerLead,leadNames)

      # Reorder waves fitted values
      maxWaves<-ncol(fittedWaves[[1]])
      if(nWaves<maxWaves){
        entireOrder<-c(orderByRelVar, (nWaves+1):(maxWaves))
      }else{
        entireOrder<-orderByRelVar
      }
      fittedWaves<-lapply(1:nSignals, function(x) fittedWaves[[x]][,entireOrder])

      # Reorder unused components
      unusedComp<-sapply(1:length(unusedComp), function(x) ifelse(orderByRelVar[unusedComp[x]]!=unusedComp[x],
                                                                  which(unusedComp[x]==orderByRelVar), unusedComp[x]))

      return(list("fittedWaves"=fittedWaves, "paramsPerLead"=paramsPerLead, "unusedComp"=unusedComp))
    }else{
      return(list("fittedWaves"=fittedWaves, "paramsPerLead"=paramsPerLead, "unusedComp"=unusedComp))
    }
  }
}

searchAdditionalComponents<-function(vDataMatrix, baseGrid , fittedWaves, leadNames, alphaN,
                                     weightError, errorWeights, signalWeights, paramsPerLead,
                                     usedApply, singleFilename, plotBacks=FALSE, commonOmega=TRUE,
                                     extraBack=6:7, lengthAlphaGrid = 48, lengthOmegaGrid = 24,
                                     alphaGrid = seq(0, 2*pi, length.out = lengthAlphaGrid),
                                     omegaMin = 0.0001, omegaMax = 1,
                                     omegaGrid = exp(seq(log(omegaMin), log(omegaMax), length.out = lengthOmegaGrid)),
                                     numReps = 3, searchAllExtraBack=FALSE){

  validComponentFound<-FALSE; currentIndex<-1; nSignals<-ncol(vDataMatrix); fmmObjectList<-list()
  unusedComp<-which(substr(rownames(paramsPerLead[[1]]),1,1)=="X")

  previousParamsPerLead<-paramsPerLead
  while(!validComponentFound){
    j=extraBack[currentIndex]

    #### First Step: determine optimal common parameters (just alpha or alpha and omega) ####
    optimalParams <- optimizeAlphaOmega(vDataMatrix = vDataMatrix, baseGrid = baseGrid, fittedWaves = fittedWaves, currentComp=j,
                                        omegaMax = omegaMax, errorWeights = errorWeights, usedApply = usedApply)
    
    #### Second Step: fit single FMM wave in each signal with common parameters ####
    for(signalIndex in 1:nSignals){
      vData<-vDataMatrix[,signalIndex]; vData<-vData[!is.na(vData)]; nObs<-length(vData)
      vData<-vData - apply(as.matrix(fittedWaves[[signalIndex]]), 1, sum)
      fmmObjectList[[signalIndex]]<-optimizeOtherParameters(vData = vData, fixedParams = optimalParams)
      
      fittedWaves[[signalIndex]][,j]<-getFittedValues(fmmObjectList[[signalIndex]])
    }

    #### Get FMM parameters per wave ####
    params<-data.frame("M"=sapply(fmmObjectList,getM),"A"=sapply(fmmObjectList,getA),"Alpha"=sapply(fmmObjectList,getAlpha),
                       "Beta"=sapply(fmmObjectList,getBeta),"Omega"=sapply(fmmObjectList,getOmega), "Var"=rep(NA,nSignals))
    paramsPerLead<-lapply(1:nSignals, function(x) rbind(paramsPerLead[[x]], params[x,], make.row.names=FALSE))
salida<-paramsPerLead
    ## Recalculate As and Ms
    paramsPerLead<-lapply(1:nSignals, function(x) recalculateMA(vDatai=vDataMatrix[,x], paramsPerSignal=paramsPerLead[[x]]))

    # If no assignation is done, it is unnecessary to order by variance
    if(!searchAllExtraBack){
      ## Order by RelVar
      paramsPerLead<-setNames(paramsPerLead,leadNames)
      paramsPerLead_fittedWaves<-wavesInRelValOrder(vDataMatrix=vDataMatrix, paramsPerLead=paramsPerLead, fittedWaves=fittedWaves,
                                                    unusedComp=unusedComp)
      paramsPerLead<-paramsPerLead_fittedWaves$paramsPerLead; fittedWaves<-paramsPerLead_fittedWaves$fittedWaves; unusedComp<-paramsPerLead_fittedWaves$unusedComp
    }

    #### Error is weighted across signals ####
    sigma<-getAssignedSigma(vDataMatrix=vDataMatrix, fittedWaves=fittedWaves, paramsPerLead=paramsPerLead,
                            weightError=weightError, signalWeights=signalWeights, assigned=FALSE)

    #### Wave preassignment ####
    paramsPerLead<-setNames(paramsPerLead,leadNames)
    # If just other waves are searched, no assignation is needed
    if(!searchAllExtraBack){
      paramsPerLead<-waveAssignation(paramsPerLead = paramsPerLead, alphaN = alphaN,
                                     nObs = nObs, unusedComp = unusedComp, inSearchAdditionalComponents=TRUE)
      ## Recalculate As and Ms
      paramsPerLead<-lapply(1:nSignals, function(x) recalculateMA(vDatai=vDataMatrix[,x], paramsPerSignal=paramsPerLead[[x]]))
      paramsPerLead<-setNames(paramsPerLead,leadNames)
    }

    if(plotBacks)
      plotMultiFMM(vDatai=vDataMatrix, fittedWaves = fittedWaves, leadNames = leadNames,
                   currentBack = j, currentBackResults = paramsPerLead, filename = singleFilename,
                   plotToFile=FALSE, unassigned=TRUE, extra=TRUE)

    currentIndex<-currentIndex+1

    if(sum(substr(rownames(paramsPerLead[[1]]),1,1)!="X")==5 |
       currentIndex>(diff(range(extraBack))+1)){
      validComponentFound<-TRUE
    }else{
      if(!substr(rownames(paramsPerLead[[1]]),1,1)[j]!="X"){ # The founded wave has not been assigned
        unusedComp<-c(unusedComp, j)
      }
    }

  }

  #### Return final error weights and total MSE ####
  sigma<-getAssignedSigma(vDataMatrix=vDataMatrix, fittedWaves=fittedWaves, paramsPerLead=paramsPerLead,
                          weightError=weightError, signalWeights=signalWeights, assigned=TRUE)
#Cambio
sigmaAll<-getAssignedSigma(vDataMatrix=vDataMatrix, fittedWaves=fittedWaves, paramsPerLead=paramsPerLead,
                            weightError=weightError, signalWeights=signalWeights, assigned=FALSE)
  totalMSE<-sum(sigma)
  
  return(list("paramsPerLead"=paramsPerLead, "extraComponent"=currentIndex+5,
              "fittedWaves"=fittedWaves, "errorWeights"=errorWeights, "totalMSE"=totalMSE, salida))
}

# Funtion to get the index of the optimal backfitting
getBestBack<-function(paramsPerLeadPerBack, alphaN, totalMSE){

  # Calculate number of assigned waves and MSE per backfitting
  nBacks<-length(paramsPerLeadPerBack)
  nAssignedWavesPerBack<-sapply(1:nBacks, function(back){
    return(sum(substr(rownames(paramsPerLeadPerBack[[back]][[1]]),1,1)!="X"))
  })
  maxAssignedWaves<-max(nAssignedWavesPerBack)
  msePerBack<-sapply(1:nBacks, function(back){
    assignedWaves<-substr(rownames(paramsPerLeadPerBack[[back]][[1]]),1,1)!="X"
    return(ifelse(sum(assignedWaves)<maxAssignedWaves, Inf, totalMSE[back]))
  })

  if(maxAssignedWaves==5){

    doubtfulPs<-sapply(1:nBacks, function(back){
      return(doubtfulP(paramsPerLeadPerBack[[back]], alphaN = alphaN))
    })

    laxDoubtfulPs<-sapply(1:nBacks, function(back){
      return(doubtfulP(paramsPerLeadPerBack[[back]], alphaN = alphaN, strict=FALSE))
    })

    # A backfitting with max assigned waves and a valid P exists
    if(any(!doubtfulPs & nAssignedWavesPerBack==max(nAssignedWavesPerBack))){
      selectedBacks<-which(nAssignedWavesPerBack==max(nAssignedWavesPerBack) & !doubtfulPs)

    # A backfitting with max assigned waves and a valid lax P exists
    }else if(any(!laxDoubtfulPs & nAssignedWavesPerBack==max(nAssignedWavesPerBack))){
      selectedBacks<-which(nAssignedWavesPerBack==max(nAssignedWavesPerBack) & !laxDoubtfulPs)

    # Backfittings with just max assigned waves exists
    }else{
      selectedBacks<-which(nAssignedWavesPerBack==max(nAssignedWavesPerBack))
    }

    # Return, among selected backfittings, the one with the lowest MSE
    return(which.min(abs(msePerBack-min(msePerBack[selectedBacks]))))

  # If not all waves have been assigned, return the one with the lowest MSE
  }else{
    return(which.min(msePerBack))
  }
}

getAssignedSigma<-function(vDataMatrix, fittedWaves, paramsPerLead, weightError,
                           signalWeights, assigned=TRUE, scaledSigma=FALSE){
  nSignals<-ncol(vDataMatrix); nObs<-nrow(vDataMatrix)

  if(assigned){
    assignedWaves<-substr(rownames(paramsPerLead[[1]]),1,1)!="X"
  }else{
    assignedWaves<-rep(TRUE, nrow(paramsPerLead[[1]]))
  }

  # if(scaledSigma){
  #   sigma<-sapply(1:nSignals, function(x){
  #     descaling<-minimax(vDataMatrix[,x])
  #     vDataMatrixDescaled<-descaling$minimax
  #     return(sum((vDataMatrixDescaled-minimax(apply(as.matrix(fittedWaves[[x]][,assignedWaves]), 1, sum),
  #                                             minY = descaling$minY, maxY = descaling$maxY)$minimax)^2)/(nObs-1))})
  # }else{
  sigma<-sapply(1:nSignals, function(x){
    sum((vDataMatrix[,x]-apply(as.matrix(fittedWaves[[x]][,assignedWaves]),1, sum))^2)/(nObs-1)})
  # }

  if(weightError){
    errorWeights<-(1/sigma)*signalWeights
  }else{
    errorWeights<-rep(1,nSignals)
  }
  return(sigma)
}



# Function to get the betas, Ms and As of the unused leads in the fitting process
getOtherBetas_MAs<-function(vDataMatrix, fragmentResults, otherLeads){
  derivedLeads<-c("III", "aVR", "aVL", "aVF"); wavesNames<-c("P","Q","R","S","T")
  nLeads<-length(otherLeads); timePoints<- seqTimes(nrow(vDataMatrix))

  # Already calculated results extraction
  alphas<-fragmentResults[1,paste0("Alpha",wavesNames)]; omegas<-fragmentResults[1,paste0("Omega",wavesNames)]
  betasI<-fragmentResults[1,paste0("Beta",wavesNames,"_I")]; betasII<-fragmentResults[1,paste0("Beta",wavesNames,"_II")]
  asI<-fragmentResults[1,paste0("A",wavesNames,"_I")]; asII<-fragmentResults[1,paste0("A",wavesNames,"_II")]
  mI<-fragmentResults[1,"M_I"]; mII<-fragmentResults[1,"M_II"]

  # Function for parameter estimation in derived leads ("III" "aVR" "aVL" "aVF")
  derivedM<-function(mI, mII, coefI, coefII) return((mI*coefI)+(mII*coefII))
  deriveBetaA<-function(betaI, betaII, aI, aII, coefI, coefII){
    lambda<-aI*coefI; mu<-aII*coefII
    aC<-lambda*cos(betaI)+mu*cos(betaII)
    bC<- -lambda*sin(betaI)-mu*sin(betaII)
    return(list("Beta"=((3*pi/2)+atan2(aC, bC))%%(2*pi), "A"=sqrt((aC^2)+(bC^2))))
  }
  derivedLeadsCoefs<-data.frame("Lead"=derivedLeads, "CoefI"=c(-1,-0.5,1,-0.5), "CoefII"=c(1,-0.5,-0.5,1))

  # Variables for other results
  mMatrix<-data.frame(matrix(NA, ncol = nLeads, nrow = 1, dimnames = list("M",otherLeads)))
  r2Matrix<-data.frame(matrix(NA, ncol = nLeads, nrow = 1, dimnames = list("R2",otherLeads)))
  aMatrix<-data.frame(matrix(NA, ncol = nLeads, nrow = 5, dimnames = list(paste0("A",wavesNames),otherLeads)))
  betaMatrix<-data.frame(matrix(NA, ncol = nLeads, nrow = 5, dimnames = list(paste0("Beta",wavesNames),otherLeads)))

  #### Parameters calculus ####
  unavailableLeads<-data.frame("Patient"=c(12722, 12722, 12722, 19299),
                               "Lead"=c("V4", "V5", "V6", "V6"))

  for(i in otherLeads){
    vData<-vDataMatrix[,i]

    if(nrow(unavailableLeads[unavailableLeads$Patient==fragmentResults$EcgId & unavailableLeads$Lead==i,])==0){
      # Leads that can be derived analytically
      if(i %in% derivedLeads & any(!is.na(betasI), na.rm = TRUE) & any(!is.na(betasII), na.rm = TRUE)){
        currentCoefs<-derivedLeadsCoefs[derivedLeadsCoefs$Lead==i,]
        for(j in 1:5){
          betaI=as.numeric(betasI[j]); betaII=as.numeric(betasII[j])
          if(!is.na(betaI)){
            derivedBetaA<-deriveBetaA(betaI=betaI, betaII=betaII, aI=as.numeric(asI[j]), aII=as.numeric(asII[j]),
                                      coefI=currentCoefs$CoefI, coefII=currentCoefs$CoefII)
            betaMatrix[j,i]<-derivedBetaA$Beta; aMatrix[j,i]<-derivedBetaA$A
          }else{
            betaMatrix[j,i]<-NA; aMatrix[j,i]<-NA
          }
        }
        mMatrix[1,i]<-derivedM(mI=mI, mII=mII, coefI=currentCoefs$CoefI, coefII=currentCoefs$CoefII)
        cosPhi <-  calculateCosPhi(alpha = alphas, beta = betaMatrix[,i], omega = omegas, timePoints = timePoints)
        mod_aMatrix<-aMatrix
        if(any(is.na(alphas))){ cosPhi[,is.na(alphas)]<-0; mod_aMatrix[is.na(alphas),i]<-0}
        r2Matrix[1,i] <-  PV(vData, pred = (cosPhi%*%mod_aMatrix[,i])+as.vector(mMatrix[1,i]))

        # Otherwise, M, A and beta estimation by minimum error
      }else{

        # Betas calculus
        betaMatrix[,i]<-sapply(1:5, function(j){
          alphaVal<-as.numeric(alphas[j]); omegaVal<-as.numeric(omegas[j])
          if(!is.na(alphaVal)){
            mobiusTerm <- 2*atan(omegaVal*tan((timePoints - alphaVal)/2))
            costStar <- cos(alphaVal + mobiusTerm); sentstar <- sin(alphaVal + mobiusTerm)
            covMatrix <- stats::cov(cbind(vData, costStar, sentstar))
            denominator <- covMatrix[2,2]*covMatrix[3,3] - covMatrix[2,3]^2
            cosCoeff <- (covMatrix[1,2]*covMatrix[3,3] -
                           covMatrix[1,3]*covMatrix[2,3])/denominator
            sinCoeff <- (covMatrix[1,3]*covMatrix[2,2] -
                           covMatrix[1,2]*covMatrix[2,3])/denominator
            return((atan2(-sinCoeff, cosCoeff)+alphaVal)%%(2*pi))
          }else{
            return(NA)
          }
        })

        # M and As calculus
        cosPhi <-  calculateCosPhi(alpha = alphas, beta = betaMatrix[,i],
                                        omega = omegas, timePoints = timePoints)
        if(any(is.na(alphas))) cosPhi[,is.na(alphas)]<-0
        linearModel <- lm(vData ~ cosPhi)

        ## Wrong beta estimation can lead to negative As => Translate beta parameter
        aValues<-as.vector(linearModel$coefficients[-1])
        aValidContion<-aValues>0 | is.na(aValues)
        if(any(!aValidContion)){
          # Change beta parameter
          betaMatrix[!aValidContion,i]<-(betaMatrix[!aValidContion,i]+pi)%%(2*pi)

          cosPhi <-  calculateCosPhi(alpha = alphas, beta = betaMatrix[,i],
                                          omega = omegas, timePoints = timePoints)
          if(any(is.na(alphas))) cosPhi[,is.na(alphas)]<-0
          linearModel <- lm(vData ~ cosPhi)
        }

        aMatrix[,i] <- as.vector(linearModel$coefficients[-1])
        mMatrix[1,i] <- as.vector(linearModel$coefficients[1])
        r2Matrix[1,i] <-  PV(vData, pred = linearModel$fitted.values)
      }

    ## Some patients have not all leads recorded
    }else{
      betaMatrix[,i]<-NA; aMatrix[,i] <- NA; mMatrix[1,i] <- NA; r2Matrix[1,i] <- NA
    }
  }

  #### Format output ####
  toWide<-function(dataframe){
    return(dataframe %>% rownames_to_column() %>% pivot_longer(-rowname) %>%
             unite(name, rowname, name, sep = "_") %>% pivot_wider())}
  wideBetas<-toWide(betaMatrix); wideMs<-toWide(mMatrix); wideAs<-toWide(aMatrix); wideR2<-toWide(r2Matrix)

  #### Output to Fragment Result ####
  newDerivedOutputs<-function(fragmentResults, wideResult){
    # Substitute missing leads values, enter derived leads values
    derivedCols<-grepl(derivedLeads[1], colnames(wideResult), fixed = TRUE) | grepl(derivedLeads[2], colnames(wideResult), fixed = TRUE) |
      grepl(derivedLeads[3], colnames(wideResult), fixed = TRUE) | grepl(derivedLeads[4], colnames(wideResult), fixed = TRUE)
    fragmentResults[colnames(wideResult)[!derivedCols]]<-wideResult[colnames(wideResult)[!derivedCols]]
    fragmentResults<-cbind(fragmentResults,wideResult[colnames(wideResult)[derivedCols]])
    return(fragmentResults)
  }
  fragmentResults<-newDerivedOutputs(fragmentResults=fragmentResults, wideResult=wideBetas)
  fragmentResults<-newDerivedOutputs(fragmentResults=fragmentResults, wideResult=wideAs)
  fragmentResults<-newDerivedOutputs(fragmentResults=fragmentResults, wideResult=wideMs)
  fragmentResults<-newDerivedOutputs(fragmentResults=fragmentResults, wideResult=wideR2)

  return(fragmentResults)
}

calculateDerivedLeads<-function(ecgId, beatId, loadedData, loadedResults){
  twelveLeads<-c("I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6")
  eightLeads<-c("I", "II", "V1", "V2", "V3", "V4", "V5", "V6"); derivedLeads<-c("III", "aVR", "aVL", "aVF")
  wavesNames<-c("P","Q","R","S","T")

  # Get estimated model parameters
  fragmentResults<-loadedResults[loadedResults$EcgId==ecgId & loadedResults$BeatId==beatId,]
  if(!is.na(fragmentResults$Superclass)){
    otherLeads<-is.na(fragmentResults[1,paste("BetaR",eightLeads,sep="_")])
    otherLeads<-twelveLeads[c(otherLeads[1:2],rep(TRUE,4),otherLeads[-c(1:2)])]

    # Load voltage traces
    readedData<-readData(loadedData=loadedData, ecgId=ecgId, beatId=beatId,
                         usedLeads = twelveLeads, otherLeadsCalculus=TRUE)
    checkedData<-checkData(readedData$vDataMatrix, otherLeadsCalculus=TRUE)
    vDataMatrix<-checkedData$vDataMatrix

    #### Calculate other parameters ####
    ## Betas, M, As, Var
    fragmentResults<-getOtherBetas_MAs(vDataMatrix, fragmentResults, otherLeads)

    # Mark derived Leads
    if(any(c("I","II") %in% otherLeads)){
      fragmentResults[paste0("der",twelveLeads)]<- matrix(twelveLeads %in% otherLeads, nrow=1)
    }else{
      fragmentResults[paste0("der",twelveLeads)]<- matrix(twelveLeads %in% otherLeads & !twelveLeads %in% derivedLeads, nrow=1)
    }
  }else{
    addedColumnsNames<-c(paste0("Beta", rep(wavesNames, 4), "_", rep(derivedLeads,each=5)),
                         paste0("M", "_", derivedLeads),
                         paste0("A", rep(wavesNames, 4), "_", rep(derivedLeads,each=5)),
                         #paste0("Var", rep(wavesNames, 4), "_", rep(derivedLeads,each=5)),
                         paste0("R2", "_", derivedLeads), paste0("der",twelveLeads))
    fragmentResults[addedColumnsNames]<- NA
  }

  return(fragmentResults)
}

calculateDerivedLeads_App<-function(fragmentResults, vDataMatrix){
  twelveLeads<-c("I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6")
  eightLeads<-c("I", "II", "V1", "V2", "V3", "V4", "V5", "V6");

  # Get estimated model parameters
  otherLeads<-is.na(fragmentResults[1,paste("BetaR",eightLeads,sep="_")])
  otherLeads<-twelveLeads[c(otherLeads[1:2],rep(TRUE,4),otherLeads[-c(1:2)])]

  #### Calculate other parameters ####
  ## Betas, M, As, Var
  fragmentResults<-getOtherBetas_MAs(vDataMatrix, fragmentResults, otherLeads)

  return(fragmentResults)
}


recalculateR2<-function(loadedData, fragmentResults){
  twelveLeads<-c("I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6")
  waveNames<-c("P","Q","R","S","T"); derivedLeads<-c("III", "aVR", "aVL", "aVF")
  ecgId<-fragmentResults$EcgId; beatId<-fragmentResults$BeatId

  if(!is.na(fragmentResults$Superclass)){
    # Load voltage traces
    readedData<-readData(loadedData=loadedData, ecgId=ecgId, beatId=beatId,
                         usedLeads = twelveLeads, otherLeadsCalculus=TRUE)
    checkedData<-checkData(readedData$vDataMatrix, otherLeadsCalculus=TRUE)
    vDataMatrix<-checkedData$vDataMatrix
    meanMatrix<-matrix(rep(apply(vDataMatrix, 2, mean), each=nrow(vDataMatrix)), nrow = nrow(vDataMatrix), ncol = ncol(vDataMatrix))

    # Generate fitted values
    unavailableLeads<-data.frame("Patient"=c(12722, 12722, 12722, 19299),
                                 "Lead"=c("V4", "V5", "V6", "V6"))
    fittedMatrix<-sapply(twelveLeads, function(y){
      if(nrow(unavailableLeads[unavailableLeads$Patient==fragmentResults$EcgId & unavailableLeads$Lead==y,])==0){
        arguments<-list("data"=vDataMatrix[,y], "M"=as.numeric(fragmentResults[paste0("M_",y)]),
                        "A"=as.numeric(fragmentResults[paste0(rep("A",5),waveNames,"_",y)]),
                        "alpha"=as.numeric(fragmentResults[paste0(rep("Alpha",5),waveNames)]),
                        "beta"=as.numeric(fragmentResults[paste0(rep("Beta",5),waveNames,"_",y)]),
                        "omega"=as.numeric(fragmentResults[paste0(rep("Omega",5),waveNames)]), "R2"=0)

        if(any(is.na(arguments$alpha))){
          arguments$A<-arguments$A[!is.na(arguments$alpha)]; arguments$beta<-arguments$beta[!is.na(arguments$alpha)]
          arguments$omega<-arguments$omega[!is.na(arguments$alpha)]; arguments$alpha<-arguments$alpha[!is.na(arguments$alpha)]
        }
        fmmConstructor(arguments)@fittedValues
      }else{
        rep(NA, length(vDataMatrix[,y]))
      }
    })

    #### R2 calculus ####
    numValue<-apply((vDataMatrix-fittedMatrix)^2, 2, sum)
    r2Measures<-1-(numValue/(apply((vDataMatrix-meanMatrix)^2, 2, sum)))
    r2PrimeMeasures<-1-(numValue/(apply((vDataMatrix)^2, 2, sum)))
    fitMeasures<-data.frame("R2"=r2Measures, "RPrime"=r2PrimeMeasures)

  }else{
    fitMeasures<-data.frame("R2"=rep(NA,12), "RPrime"=rep(NA,12))
    rownames(fitMeasures)<-twelveLeads
  }

  #### Save in fragmentResults ####
  toWide<-function(dataframe){
    return(dataframe %>% rownames_to_column() %>% pivot_longer(-rowname) %>%
             unite(name, rowname, name, sep = "_") %>% pivot_wider())}
  wideFits<-toWide(as.data.frame(t(fitMeasures)))

  newDerivedOutputs<-function(fragmentResults, wideResult){
    # Substitute missing leads values, enter derived leads values
    alreadyInColumns<-colnames(wideResult) %in% colnames(fragmentResults)
    fragmentResults[colnames(wideResult)[alreadyInColumns]]<-wideResult[1,alreadyInColumns]
    fragmentResults<-cbind(fragmentResults,wideResult[1,!alreadyInColumns])
    return(fragmentResults)
  }

  return(newDerivedOutputs(fragmentResults=fragmentResults, wideResult=wideFits))
}

#### Useful functions related to FMM parameters ####
isAlpha1LeftOfAlpha2<-function(alpha1, alpha2){
  if(all(!is.na(alpha1)) & all(!is.na(alpha2))){
    return((alpha2>pi & (alpha1<alpha2 & alpha1>pi)) |
             (alpha2<pi & (alpha1>pi | alpha1<alpha2)))
  }else{
    return(TRUE)
  }
}

isBetaUpDown_Positive<-function(beta){
  return(2*pi/3 <= beta & beta <= 5*pi/3)
}

isBetaNegative<-function(beta){
  return(as.numeric(5*pi/3 < beta | beta < pi/3))
}

betaCharacteristics<-function(betaMatrix){
  nullResult<-rep(0, nrow(betaMatrix))

  # Beta positive is calculated for I, II and V5
  otherLeads<-colnames(betaMatrix)!="V2"; nOtherLeads<-sum(otherLeads)

  if(nOtherLeads==0){
    betaPositive_UpDown<-nullResult
  }else if(nOtherLeads==1){
    betaPositive_UpDown<-as.numeric(isBetaUpDown_Positive(betaMatrix[,otherLeads]))
  }else{
    betaPositive_UpDown<-apply(isBetaUpDown_Positive(betaMatrix[,otherLeads]), 1, sum)
  }

  # Beta negative is calculated for V2
  v2Lead<-colnames(betaMatrix)=="V2"

  if(any(v2Lead, na.rm = TRUE)){
    betaNegative<-isBetaNegative(betaMatrix[,v2Lead])
  }else{
    betaNegative<-nullResult
  }

  return(list("betaPositive_UpDown"=betaPositive_UpDown, "betaNegative"=betaNegative))
}

#### Wave sets definitions functions ####
getNoisySet<-function(paramsPerLead){
  # Parameters extraction
  omega<-paramsPerLead[[1]]$Omega; alpha<-paramsPerLead[[1]]$Alpha

  # Set's conditions
  return(noisyAlpha[1]<alpha & alpha<noisyAlpha[2] | omega>noisyMaxOmega | omega<noisyMinOmega)
}

getFreeSet<-function(paramsPerLead, unusedComp=NA, freeWaves=NA){
  nWaves<-nrow(paramsPerLead[[1]])

  freeComponents<-rep(TRUE, nWaves)

  if(!any(is.na(unusedComp))){
    freeComponents<- freeComponents & (!1:nWaves %in% unusedComp)
  }
  if(!any(is.na(freeWaves))){
    freeComponents<- freeComponents & (1:nWaves %in% freeWaves)
  }

  return(freeComponents)
}

getSuperSetR<-function(paramsPerLead, alphaN, alpha=paramsPerLead[[1]]$Alpha, minusSuperSetR=FALSE, unusedComp=NA){

  # Parameters extraction
  omega<-paramsPerLead[[1]]$Omega; distancesToAlphaN<-1-cos(alpha-alphaN)
  availableRelevantLeads<-relevantLeads[relevantLeads %in% names(paramsPerLead)]
  relevantVar<-apply(sapply(availableRelevantLeads,function(x) paramsPerLead[[x]]$Var), 1, mean)

  # SuperSetR's conditions
  ## Within the top three variance waves
  if(!minusSuperSetR){
    varCondition<-relevantVar %in% sort(relevantVar, decreasing = TRUE)[1:3]
  }else{
    varCondition<-rep(TRUE, nrow(paramsPerLead[[1]]))
  }

  ## Omega and Var condition
  omegaVarCondition<-omega<rMaxOmega[1] |
    (rMaxOmega[1] < omega & omega < rMaxOmega[2] & rMinVarMaxOmega < relevantVar)

  ## Distance conditions
  leftWaves<-isAlpha1LeftOfAlpha2(alpha,alphaN)
  distCondition<-((leftWaves & distancesToAlphaN < rLeftDistanceToAlphaN[3]) |
                    (!leftWaves & distancesToAlphaN < rRightDistanceToAlphaN))

  return(!getNoisySet(paramsPerLead) & getFreeSet(paramsPerLead, unusedComp) &
           varCondition & omegaVarCondition & distCondition)
}

getSetR<-function(paramsPerLead, alphaN, alpha=paramsPerLead[[1]]$Alpha, setPlusR=NA, unusedComp=NA){

  # Parameters extraction
  distancesToAlphaN<-1-cos(alpha-alphaN)
  availableRelevantLeads<-relevantLeads[relevantLeads %in% names(paramsPerLead)]
  relevantBetaMatrix<-sapply(availableRelevantLeads,function(x) paramsPerLead[[x]]$Beta)
  leftWaves<-isAlpha1LeftOfAlpha2(alpha,alphaN)

  # Extract beta characteristics
  betaCharacs<-betaCharacteristics(relevantBetaMatrix)
  betaPositive_UpDown<-betaCharacs$betaPositive_UpDown; betaNegative<-betaCharacs$betaNegative

  # Definition of the R sets
  superSetR<-getSuperSetR(paramsPerLead=paramsPerLead, alpha=alpha, alphaN=alphaN, unusedComp=unusedComp)

  if(is.na(setPlusR)){
    setR<-superSetR & ((leftWaves & distancesToAlphaN < rLeftDistanceToAlphaN[1]) | !leftWaves) &
      ((betaPositive_UpDown+betaNegative)>=2)
    return(setR)
  }else if(setPlusR==1){
    setPlusR<-superSetR & ((leftWaves & distancesToAlphaN < rLeftDistanceToAlphaN[2]) | !leftWaves) &
      ((betaPositive_UpDown+betaNegative)>=2)
    return(setPlusR)
  }else if(setPlusR==2){
    setPlusPlusR<-superSetR & ((betaPositive_UpDown+betaNegative)>=1)
    return(setPlusPlusR)
  }else{
    stop("Call of getSetR not correctly specified.")
  }
}

getSuperSetP<-function(paramsPerLead, alphaN, rIndex=which(substr(rownames(paramsPerLead[[1]]),1,1)=="R"),
                       alpha=paramsPerLead[[1]]$Alpha, unusedComp=NA, freeWaves=NA){

  # Parameters extraction
  omega<-paramsPerLead[[1]]$Omega; alphaR<-alpha[rIndex]
  availableRelevantLeads<-relevantLeads[relevantLeads %in% names(paramsPerLead)]
  relevantVar<-apply(sapply(availableRelevantLeads,function(x) paramsPerLead[[x]]$Var), 1, mean)
  distanceToAlphaR<-1-cos(alphaR-alpha); distanceToAlphaR[rIndex]<-NA
  distancesToAlphaN<-1-cos(alpha-alphaN)

  # SuperSetP's conditions
  leftOfAlphaR <- isAlpha1LeftOfAlpha2(alpha1 = alpha, alpha2 = alphaR); leftOfAlphaR[rIndex]<-NA

  ## Omega and Var condition
  omegaVarCondition<-omega < pMaxOmega[1] | (pMaxOmega[2] > omega & omega > pMaxOmega[1] & pMinVarMaxOmega < relevantVar)

  ## Wave position and distance conditions
  distCondition<-maxDistPR[2] > distanceToAlphaR & distanceToAlphaR>minDistPR[1] &
    distancesToAlphaN>minDistPN

  return(!getNoisySet(paramsPerLead) & getFreeSet(paramsPerLead, unusedComp, freeWaves) &
           leftOfAlphaR & distCondition & omegaVarCondition)
}

getSetPlusP<-function(paramsPerLead, alphaN, setPlusPlusP=FALSE, rIndex=which(substr(rownames(paramsPerLead[[1]]),1,1)=="R"),
                      alpha=paramsPerLead[[1]]$Alpha, unusedComp=NA, freeWaves=NA){

  # Parameters extraction
  alphaR<-alpha[rIndex]; omega<-paramsPerLead[[1]]$Omega
  availableRelevantLeads<-relevantLeads[relevantLeads %in% names(paramsPerLead)]
  relevantVar<-apply(sapply(availableRelevantLeads,function(x) paramsPerLead[[x]]$Var), 1, mean)
  relevantBetaMatrix<-sapply(availableRelevantLeads,function(x) paramsPerLead[[x]]$Beta)
  betaPositive_UpDown<-betaCharacteristics(relevantBetaMatrix)$betaPositive_UpDown
  distanceToAlphaR<-1-cos(alphaR-alpha); distanceToAlphaR[rIndex]<-NA

  # Minimum conditions
  superSetP<-getSuperSetP(paramsPerLead = paramsPerLead, alphaN = alphaN, alpha = alpha,
                          rIndex=rIndex, unusedComp=unusedComp, freeWaves = freeWaves)
  omegaVarCondition<-(0.01<omega & omega<0.2 & relevantVar>0.001) | (omega>0.2 & relevantVar>0.01)
  basicConditions<-superSetP & omegaVarCondition

  # Depending on the set, the distances restrictions between P and R are different
  doubtDistCondition<-distanceToAlphaR < minDistPR[2] & betaPositive_UpDown>=2

  if(!setPlusPlusP){
    safeDistCondition<-maxDistPR[1] > distanceToAlphaR & distanceToAlphaR > minDistPR[2]
  }else{
    safeDistCondition<-distanceToAlphaR > minDistPR[2]
  }

  return(basicConditions & (doubtDistCondition|safeDistCondition))
}

getSuperSetQ<-function(paramsPerLead, rIndex=which(substr(rownames(paramsPerLead[[1]]),1,1)=="R"),
                       alpha=paramsPerLead[[1]]$Alpha, unusedComp=NA, freeWaves=NA){

  # Parameters extraction
  omega<-paramsPerLead[[1]]$Omega; alphaR<-alpha[rIndex]
  availableRelevantLeads<-relevantLeads[relevantLeads %in% names(paramsPerLead)]
  relevantVar<-apply(sapply(availableRelevantLeads,function(x) paramsPerLead[[x]]$Var), 1, mean)
  distanceToAlphaR<-1-cos(alphaR-alpha); distanceToAlphaR[rIndex]<-NA

  # SuperSetQ's conditions
  leftOfAlphaR <- isAlpha1LeftOfAlpha2(alpha1 = alpha, alpha2 = alphaR); leftOfAlphaR[rIndex]<-NA

  ## Omega and Var condition
  omegaVarCondition<-omega < qMaxOmega[1] | (qMaxOmega[2] > omega & omega > qMaxOmega[1] & qMinVarMaxOmega < relevantVar)

  ## Distance conditions
  distCondition<-distanceToAlphaR<maxDistQR[2]

  return(!getNoisySet(paramsPerLead) & getFreeSet(paramsPerLead, unusedComp, freeWaves) &
           leftOfAlphaR & distCondition & omegaVarCondition)
}

getSetQ<-function(paramsPerLead, rIndex=which(substr(rownames(paramsPerLead[[1]]),1,1)=="R"),
                  alpha=paramsPerLead[[1]]$Alpha, unusedComp=NA, freeWaves=NA){

  # Parameters extraction
  alphaR<-alpha[rIndex]; distanceToAlphaR<-1-cos(alphaR-alpha); distanceToAlphaR[rIndex]<-NA

  return(getSuperSetQ(paramsPerLead=paramsPerLead, rIndex=rIndex, alpha=alpha,
                      unusedComp=unusedComp, freeWaves=freeWaves) &
           distanceToAlphaR < maxDistQR[1])
}

getSuperSetS<-function(paramsPerLead, alphaN, rIndex=which(substr(rownames(paramsPerLead[[1]]),1,1)=="R"),
                       alpha=paramsPerLead[[1]]$Alpha, unusedComp=NA, freeWaves=NA){

  # Parameters extraction
  omega<-paramsPerLead[[1]]$Omega; alphaR<-alpha[rIndex]
  availableRelevantLeads<-relevantLeads[relevantLeads %in% names(paramsPerLead)]
  relevantVar<-apply(sapply(availableRelevantLeads,function(x) paramsPerLead[[x]]$Var), 1, mean)
  distanceToAlphaR<-1-cos(alphaR-alpha); distanceToAlphaR[rIndex]<-NA
  distancesToAlphaN<-1-cos(alpha-alphaN); distancesToAlphaN[rIndex]<-NA

  # SuperSetS's conditions
  rightOfAlphaR <- isAlpha1LeftOfAlpha2(alpha1 = alphaR, alpha2 = alpha); rightOfAlphaR[rIndex]<-NA

  ## Omega and Var condition
  omegaVarCondition<-omega < sMaxOmega[1] | (sMaxOmega[2] > omega & omega > sMaxOmega[1] & sMinVarMaxOmega < relevantVar)

  ## Distance conditions
  distCondition<-distanceToAlphaR<maxDistRS | distancesToAlphaN<maxDistSN

  return(!getNoisySet(paramsPerLead) & getFreeSet(paramsPerLead, unusedComp, freeWaves) &
           rightOfAlphaR & distCondition & omegaVarCondition)
}

getSuperSetT<-function(paramsPerLead, alphaN, rIndex=which(substr(rownames(paramsPerLead[[1]]),1,1)=="R"),
                       alpha=paramsPerLead[[1]]$Alpha, unusedComp=NA, freeWaves=NA){

  # Parameters extraction
  alphaR<-alpha[rIndex]; distanceToAlphaR<-1-cos(alphaR-alpha); distanceToAlphaR[rIndex]<-NA
  distancesToAlphaN<-1-cos(alpha-alphaN)

  # SuperSetT's conditions
  rightOfAlphaR <- isAlpha1LeftOfAlpha2(alpha1 = alphaR, alpha2 = alpha); rightOfAlphaR[rIndex]<-NA

  ## Distance conditions
  distCondition<-distanceToAlphaR>minDistRT[1] & distancesToAlphaN>minDistNT

  return(!getNoisySet(paramsPerLead) & getFreeSet(paramsPerLead, unusedComp, freeWaves) &
           rightOfAlphaR & distCondition)
}

getSetT<-function(paramsPerLead, alphaN, rIndex=which(substr(rownames(paramsPerLead[[1]]),1,1)=="R"),
                  alpha=paramsPerLead[[1]]$Alpha, unusedComp=NA, freeWaves=NA){

  # Parameters extraction
  alphaR<-alpha[rIndex]; distanceToAlphaR<-1-cos(alphaR-alpha); distanceToAlphaR[rIndex]<-NA

  return(getSuperSetT(paramsPerLead=paramsPerLead, alphaN=alphaN, rIndex=rIndex, alpha=alpha,
                      unusedComp=unusedComp, freeWaves=freeWaves) & distanceToAlphaR > minDistRT[2])
}

getIndexU<-function(paramsPerLead, unusedComp=NA, freeWaves=NA){

  uCondition<-!getNoisySet(paramsPerLead) & getFreeSet(paramsPerLead, unusedComp, freeWaves)

  return(ifelse(any(uCondition, na.rm = TRUE), which(uCondition)[1], NA))
}

#### Wave assignment functions ####
getPairMaxVar<-function(set1, set2, relevantVar, alpha){

  # Get all possible combinations
  setsCombinations<-expand.grid(set1Index=c(which(set1),NA),set2Index=c(which(set2),NA))
  validIndexes<-ifelse(is.na(setsCombinations$set1Index) & is.na(setsCombinations$set2Index), FALSE,
                       ifelse(is.na(setsCombinations$set1Index) | is.na(setsCombinations$set2Index), TRUE,
                              setsCombinations$set1Index!=setsCombinations$set2Index))

  # alphaSet1 < alphaSet2 to be a valid combination
  validIndexes<-validIndexes & sapply(1:nrow(setsCombinations),
                                      function(x) isAlpha1LeftOfAlpha2(alpha1=alpha[setsCombinations[x,"set1Index"]],
                                                                       alpha2=alpha[setsCombinations[x,"set2Index"]]))

  # Evaluate the explained variance by the combinations
  setsCombinations<-setsCombinations[validIndexes,]
  explainedVarSets<-sapply(1:nrow(setsCombinations), function(x)
    ifelse(is.na(setsCombinations[x,"set1Index"]), relevantVar[setsCombinations[x,"set2Index"]],
           ifelse(is.na(setsCombinations[x,"set2Index"]), relevantVar[setsCombinations[x,"set1Index"]],
                  relevantVar[setsCombinations[x,"set1Index"]]+relevantVar[setsCombinations[x,"set2Index"]])))
  indexInSet1<-setsCombinations[which.max(explainedVarSets),"set1Index"]
  indexInSet2<-setsCombinations[which.max(explainedVarSets),"set2Index"]
  return(c(indexInSet1, indexInSet2))
}

rWaveAssignation<-function(paramsPerLead, alphaN, alpha=paramsPerLead[[1]]$Alpha, unusedComp=NA){

  # Parameters extraction
  availableRelevantLeads<-relevantLeads[relevantLeads %in% names(paramsPerLead)]
  relevantVar<-apply(sapply(availableRelevantLeads,function(x) paramsPerLead[[x]]$Var), 1, mean)

  # Definition of the R sets
  superSetR<-getSuperSetR(paramsPerLead=paramsPerLead, alphaN=alphaN, alpha=alpha, unusedComp=unusedComp)
  setR<-getSetR(paramsPerLead = paramsPerLead, alphaN = alphaN, alpha=alpha, unusedComp = unusedComp)
  setPlusR<-getSetR(paramsPerLead = paramsPerLead, setPlusR = 1, alphaN = alphaN, alpha=alpha, unusedComp = unusedComp)
  setPlusPlusR<-getSetR(paramsPerLead = paramsPerLead, setPlusR = 2, alphaN = alphaN, alpha=alpha, unusedComp = unusedComp)

  # Evaluate each set until any is not empty
  if(any(setR, na.rm = TRUE)){
    rIndex <- which(setR)[1]
  }else if(any(setPlusR, na.rm = TRUE)){
    rIndex <- which(setPlusR)[1]
  }else if(any(setPlusPlusR, na.rm = TRUE)){
    rIndex <- which(setPlusPlusR)[1]
  }else if(any(superSetR, na.rm = TRUE)){
    rIndex <- which(superSetR)[1]
  }else{ # Else, R is the most relevant wave
    rIndex<-which(relevantVar %in% sort(relevantVar, decreasing = TRUE)[1])[1]
  }

  return(rIndex)
}

uWaveCheck<-function(paramsPerLead, wavesNames, alphaN, uIndex, qIndex, rIndex, sIndex,
                     alpha=paramsPerLead[[1]]$Alpha){

  # Parameters extraction
  nSignals<-length(paramsPerLead); alphaR<-alpha[rIndex]
  availableRelevantLeads<-relevantLeads[relevantLeads %in% names(paramsPerLead)]
  relevantVarMatrix<-sapply(availableRelevantLeads,function(x) paramsPerLead[[x]]$Var)
  relevantVar<-apply(relevantVarMatrix, 1, mean); uChanges<-""

  ## Too close alphas to alphaR
  if(1-cos(alphaR-alpha[uIndex])<1-cos(0.01)){
    # This wave may be Q or S
    if(is.na(qIndex)){
      qIndex<-uIndex; wavesNames[qIndex]<-"Q"
      for(i in 1:nSignals) paramsPerLead[[i]]$Alpha[qIndex]<-alpha[rIndex]
      uChanges<-"U->Q"
    }else if(is.na(sIndex)){
      sIndex<-uIndex; wavesNames[sIndex]<-"S"
      for(i in 1:nSignals) paramsPerLead[[i]]$Alpha[sIndex]<-alpha[rIndex]
      uChanges<-"U->S"
    }else{
      # S must be deleted if: A) VarU>VarS and VarU<VarQ B) VarU>VarQ, VarU>VarS and VarS<VarQ
      # Q must be deleted if: A) VarU>VarQ and VarU<VarS B) VarU>VarQ, VarU>VarS and VarS>VarQ
      deleteConditionS<-(relevantVar[uIndex]<relevantVar[qIndex] & relevantVar[uIndex]>relevantVar[sIndex]) |
        (relevantVar[uIndex]>relevantVar[qIndex] & relevantVar[uIndex]>relevantVar[sIndex] & relevantVar[sIndex]<relevantVar[qIndex])
      deleteConditionQ<-(relevantVar[uIndex]>relevantVar[qIndex] & relevantVar[uIndex]<relevantVar[sIndex]) |
        (relevantVar[uIndex]>relevantVar[qIndex] & relevantVar[uIndex]>relevantVar[sIndex] & relevantVar[sIndex]>relevantVar[qIndex])
      waveToBeDeleted<-ifelse(deleteConditionS, "S", ifelse(deleteConditionQ, "Q", ""))

      if(waveToBeDeleted=="Q"){
        wavesNames[qIndex]<-"XQ1"; qIndex<-uIndex; wavesNames[qIndex]<-"Q"
        for(i in 1:nSignals) paramsPerLead[[i]]$Alpha[qIndex]<-alpha[rIndex]
        uChanges<-"U->Q"
      }else if(waveToBeDeleted=="S"){
        wavesNames[sIndex]<-"XS"; sIndex<-uIndex; wavesNames[sIndex]<-"S"
        for(i in 1:nSignals) paramsPerLead[[i]]$Alpha[sIndex]<-alpha[rIndex]
        uChanges<-"U->S"
      } # If neither, the U wave is left as it is
    }

  ## Can U be R?
  }else{
    uInSuperSetR<-getSuperSetR(paramsPerLead = paramsPerLead, alpha = alpha, alphaN=alphaN, minusSuperSetR=TRUE)[uIndex]
    rInSuperSetQ<-getSuperSetQ(paramsPerLead = paramsPerLead, alpha = alpha, rIndex = uIndex)[rIndex]
    rInSuperSetS<-getSuperSetS(paramsPerLead = paramsPerLead, alpha = alpha, alphaN = alphaN, rIndex = uIndex)[rIndex]
    minVarCondition<-relevantVar[uIndex]>minVarU_to_R

    # U is R, R is Q instead
    if(uInSuperSetR & rInSuperSetQ & minVarCondition){
      changeCondition<-relevantVar[uIndex]>relevantVar[rIndex] & ifelse(!is.na(qIndex), relevantVar[uIndex]>relevantVar[qIndex], TRUE) &
        ifelse(!is.na(sIndex), getSuperSetS(paramsPerLead = paramsPerLead, alpha = alpha, alphaN=alphaN, rIndex = uIndex)[sIndex], TRUE)
      if(changeCondition){
        if(!is.na(qIndex)) wavesNames[qIndex]<-"XQ1"
        wavesNames[uIndex]<-"R"; wavesNames[rIndex]<-"Q"
        uChanges<-"U->R"
      }

    # U is R, R is S instead
    }else if(uInSuperSetR & rInSuperSetS & minVarCondition){
      changeCondition<-relevantVar[uIndex]>relevantVar[rIndex] & ifelse(!is.na(sIndex), relevantVar[uIndex]>relevantVar[sIndex], TRUE) &
        ifelse(!is.na(qIndex), getSuperSetQ(paramsPerLead = paramsPerLead, alpha = alpha, rIndex = uIndex)[qIndex], TRUE)
      if(changeCondition){
        if(!is.na(sIndex)) wavesNames[sIndex]<-"XS"
        wavesNames[uIndex]<-"R"; wavesNames[rIndex]<-"S"
        uChanges<-"U->R"
      }
    }
  }
  return(list("wavesNames"=wavesNames, "uChanges"=uChanges, "paramsPerLead"=paramsPerLead))
}

waveAssignation<-function(paramsPerLead, alphaN, nObs, unusedComp=NA, inSearchAdditionalComponents=FALSE){
  nWaves<-nrow(paramsPerLead[[1]]); nSignals<-length(paramsPerLead)

  # Parameters extraction
  omega<-paramsPerLead[[1]]$Omega; alpha<-paramsPerLead[[1]]$Alpha
  availableRelevantLeads<-relevantLeads[relevantLeads %in% names(paramsPerLead)]
  relevantBetaMatrix<-sapply(availableRelevantLeads,function(x) paramsPerLead[[x]]$Beta)
  relevantVarMatrix<-sapply(availableRelevantLeads,function(x) paramsPerLead[[x]]$Var)
  relevantVar<-apply(relevantVarMatrix, 1, mean)

  ## In case of big omega values, tU is used instead
  isOmegaBig<-omega>0.2
  if(any(isOmegaBig)){
    bigOmegas<-which(isOmegaBig)
    for(i in bigOmegas){
      betaValues<-sapply(1:nSignals, function(x) paramsPerLead[[x]][i,"Beta"])
      if(any((2*pi/3)<=betaValues & betaValues<=(4*pi/3), na.rm=TRUE)){
        mostMaxBeta<-which.min(abs(betaValues-pi))
        fmmPeaks<-getFMMPeaksParams(currentBackResults = paramsPerLead[[mostMaxBeta]], nObs=nObs)
        alpha[i]<-((fmmPeaks$PeakU[i])+pi)%%(2*pi)
      }else if(any(((7*pi/4)<=betaValues & betaValues<=(2*pi))|(betaValues<=(pi/4)), na.rm=TRUE)){
        # If the wave has a clear minimum, use it
        mostMinBeta<-which.min(apply(cbind(betaValues,2*pi-betaValues),1,min))
        fmmPeaks<-getFMMPeaksParams(currentBackResults = paramsPerLead[[mostMinBeta]], nObs=nObs)
        alpha[i]<-((fmmPeaks$PeakL[i])+pi)%%(2*pi)
      }
    }
  }

  ## For components 6-7, at least RelExtraVar=15%
  if(inSearchAdditionalComponents & nWaves>5){
    alreadyExplainedVar<-sum(relevantVar[1:5])
    for(j in 6:(nWaves)){
      if((relevantVar[j]/(1-alreadyExplainedVar))<minExtraVar){
        unusedComp<-j
      }
    }
  }

  #### 1. R Assignment ####
  rIndex<-rWaveAssignation(paramsPerLead=paramsPerLead, alpha=alpha, alphaN=alphaN,
                           unusedComp=unusedComp)
  rIndexIsValid<-FALSE; rChangesCount<-0

  # While related to the R checks
  while(!rIndexIsValid & rChangesCount<2){
    wavesNames<-paste("X",1:nWaves,sep=""); freeWaves<-1:nWaves
    pIndex<-NA; qIndex<-NA; sIndex<-NA; tIndex<-NA
    wavesNames[rIndex]<-"R"
    freeWaves<-freeWaves[-match(rIndex, freeWaves)]

    #### 2. Assign the rest of the waves ####
    alphaR<-alpha[rIndex]
    leftOfAlphaR <- isAlpha1LeftOfAlpha2(alpha,alphaR)
    rightOfAlphaR <- isAlpha1LeftOfAlpha2(alphaR, alpha)

    #### Left Waves assignation ####
    if(any(leftOfAlphaR, na.rm = TRUE)){

      ## P and Q: assigned as tuple
      superSetQ<-getSuperSetQ(paramsPerLead=paramsPerLead, rIndex=rIndex, alpha=alpha,
                              unusedComp=unusedComp, freeWaves=freeWaves)
      setQ<-getSetQ(paramsPerLead=paramsPerLead, rIndex=rIndex, alpha=alpha,
                    unusedComp=unusedComp, freeWaves=freeWaves)
      superSetP<-getSuperSetP(paramsPerLead=paramsPerLead, rIndex=rIndex, alpha=alpha,
                              alphaN=alphaN, unusedComp=unusedComp, freeWaves=freeWaves)
      setPlusP<-getSetPlusP(paramsPerLead=paramsPerLead, alphaN=alphaN, setPlusPlusP = FALSE,
                            rIndex=rIndex, alpha=alpha, unusedComp=unusedComp, freeWaves=freeWaves)
      setPlusPlusP<-getSetPlusP(paramsPerLead=paramsPerLead, alphaN=alphaN, setPlusPlusP = TRUE,
                                rIndex=rIndex, alpha=alpha, unusedComp=unusedComp, freeWaves=freeWaves)

      if(any(superSetP, na.rm = TRUE) | any(setQ, na.rm = TRUE)){

        # If there is a single wave in both setP and setQ, it is decided what is
        if((sum(superSetP, na.rm = TRUE)==1 & sum(setQ, na.rm = TRUE)==1 &
            sum(superSetP & setQ, na.rm = TRUE)==1) |
           (sum(superSetP, na.rm = TRUE)==1 & sum(superSetQ, na.rm = TRUE)==1 &
            sum(superSetP & superSetQ, na.rm = TRUE)==1)){

          if(any(setPlusP, na.rm = TRUE)){
            pIndex<-which(setPlusP)
          }else{
            qIndex<-ifelse(any(setQ), which(setQ), which(superSetQ))
          }

        # Otherwise, select the pair P,Q that explain the most variability
        }else{
          ## The P set used is the most constrained one with at least a candidate
          if(any(setPlusPlusP, na.rm = TRUE) & !all(setPlusPlusP==setQ)){
            chosenSetP<-setPlusPlusP
          }else if(any(setPlusP, na.rm = TRUE) & !all(setPlusP==setQ)){
            chosenSetP<-setPlusP
          }else{
            chosenSetP<-superSetP
          }

          ## If no Q candidates are found, use set Q minus
          if(!any(setQ, na.rm = TRUE)){
            selectedPair<-getPairMaxVar(set1=chosenSetP, set2=superSetQ, relevantVar=relevantVar, alpha = alpha)
          }else{
            selectedPair<-getPairMaxVar(set1=chosenSetP, set2=setQ, relevantVar=relevantVar, alpha = alpha)
          }
          pIndex<-selectedPair[1]; qIndex<-selectedPair[2]
        }
        if(!is.na(pIndex)){
          wavesNames[pIndex]<-"P"; freeWaves<-freeWaves[-match(pIndex, freeWaves)]
        }
        if(!is.na(qIndex)){
          wavesNames[qIndex]<-"Q"; freeWaves<-freeWaves[-match(qIndex, freeWaves)]
        }
      }
    }

    #### Right Waves assignation ####
    if(any(rightOfAlphaR, na.rm = TRUE)){

      ## S and T: assigned as tuple
      superSetS<-getSuperSetS(paramsPerLead = paramsPerLead, alphaN = alphaN, rIndex = rIndex,
                              alpha=alpha, unusedComp = unusedComp, freeWaves = freeWaves)
      superSetT<-getSuperSetT(paramsPerLead = paramsPerLead, alphaN = alphaN, rIndex = rIndex,
                              alpha=alpha, unusedComp = unusedComp, freeWaves = freeWaves)
      setT<-getSetT(paramsPerLead = paramsPerLead, alphaN = alphaN, rIndex = rIndex,
                    alpha=alpha, unusedComp = unusedComp, freeWaves = freeWaves)

      if(any(superSetS, na.rm = TRUE) | any(setT, na.rm = TRUE)){

        # If there is a single wave in both setS and setT, it is decided what is
        if((sum(superSetS, na.rm = TRUE)==1 & sum(setT, na.rm = TRUE)==1 &
            sum(superSetS & setT, na.rm = TRUE)==1) |
           (sum(superSetS, na.rm = TRUE)==1 & sum(superSetT, na.rm = TRUE)==1 &
            sum(superSetS & superSetT, na.rm = TRUE)==1)){

          singleCandidate<-which(superSetS)
          if(omega[singleCandidate]>omegaS_or_T){
            tIndex<-singleCandidate
          }else{
            sIndex<-singleCandidate
          }

        # Otherwise, select the pair S, T that explain the most variability
        }else{
          ## If no T candidates are found, use set T minus
          if(!any(setT, na.rm = TRUE)){
            selectedPair<-getPairMaxVar(set1=superSetS, set2=superSetT, relevantVar=relevantVar, alpha = alpha)
          }else{
            selectedPair<-getPairMaxVar(set1=superSetS, set2=setT, relevantVar=relevantVar, alpha = alpha)
          }
          sIndex<-selectedPair[1]; tIndex<-selectedPair[2]
        }
        if(!is.na(sIndex)){
          wavesNames[sIndex]<-"S"; freeWaves<-freeWaves[-match(sIndex, freeWaves)]
        }
        if(!is.na(tIndex)){
          wavesNames[tIndex]<-"T"; freeWaves<-freeWaves[-match(tIndex, freeWaves)]
        }
      }
    }

    #### R Wave Check ####
    uIndex<-getIndexU(paramsPerLead = paramsPerLead, freeWaves = freeWaves, unusedComp = unusedComp)
    varR<-relevantVar[rIndex]
    superSetR<-getSuperSetR(paramsPerLead=paramsPerLead, alpha=alpha, alphaN=alphaN, unusedComp=unusedComp)

    # Can P really be R?
    canPbeR<-FALSE
    if(!is.na(pIndex)){
      canPbeR<-lowVarR_P>varR & superSetR[pIndex] & varR<relevantVar[pIndex]
    }

    # Can Q really be R?
    canQbeR<-FALSE
    if(!is.na(qIndex)){
      if((lowVarR_Q[1]>varR) | (lowVarR_Q[2]>varR & qIndex==which.max(relevantVar))){
        # Condition 1
        canQbeR<-superSetR[qIndex] & varR<relevantVar[qIndex]

        # Condition 2
        isQinSetPlusPlusR<-getSetR(paramsPerLead = paramsPerLead, alpha=alpha, alphaN = alphaN, setPlusR = 2, unusedComp = unusedComp)[qIndex]
        canQbeR<-canQbeR | (superSetR[qIndex] & is.na(sIndex) & relevantVar[qIndex]>minVarQ_to_R & isQinSetPlusPlusR)
      }
    }

    # Can S really be R?
    canSbeR<-FALSE
    if(!is.na(sIndex) & lowVarR_S>varR){
      betaCharacs<-betaCharacteristics(relevantBetaMatrix)
      betaPositive_UpDown<-betaCharacs$betaPositive_UpDown; betaNegative<-betaCharacs$betaNegative

      ## Point 3.1.1. and 3.1.2.
      commonConditions<-superSetR[sIndex] & varR<relevantVar[sIndex]

      # Condition 1
      canSbeR<-betaPositive_UpDown[sIndex]>=2

      # Condition 2
      existsAnotherS<-FALSE
      updatedSuperSetS<-getSuperSetS(paramsPerLead = paramsPerLead, alpha = alpha, alphaN=alphaN, rIndex = sIndex)
      sOtherIndex<-ifelse(any(updatedSuperSetS, na.rm=TRUE) & which(updatedSuperSetS)[1]<=5,
                          which(updatedSuperSetS)[1], NA)

      if(!is.na(sOtherIndex)){
        if(!is.na(qIndex)){
          existsAnotherS<-existsAnotherS & relevantVar[sOtherIndex]>relevantVar[qIndex]
        }
      }
      canSbeR<-commonConditions & (canSbeR | existsAnotherS)

      ## Point 3.1.3.
      nPositiveRelevantLeads<-length(availableRelevantLeads[availableRelevantLeads %in% c("I","II","V5")])
      # omegaVarConditionR<-omega[sIndex]<rMaxOmega[1] |
      #   (rMaxOmega[1] < omega[sIndex] & omega[sIndex] < rMaxOmega[2] & rMinVarMaxOmega < relevantVar[sIndex])
      # This change is evaluated only if 3 or 2 positive relevant leads are available
      positiveBetasCondition<-betaPositive_UpDown[sIndex]>=ifelse(nPositiveRelevantLeads==2, 1, 2)

      canSbeR<-canSbeR | (varR<relevantVar[sIndex] & positiveBetasCondition & betaNegative[sIndex]==1)
    }

    canUbeQ<-FALSE; canUbeS<-FALSE
    if(!is.na(uIndex) & (!is.na(qIndex) | !is.na(sIndex))){
      if(!is.na(qIndex)){
        # Can Q be R, U be Q?
        rInSuperSetS<-getSuperSetS(paramsPerLead = paramsPerLead, alpha = alpha, alphaN=alphaN, rIndex = qIndex)[rIndex]
        uInSuperSetQ<-getSuperSetQ(paramsPerLead=paramsPerLead, rIndex=qIndex, alpha=alpha)[uIndex]
        conditionVarS<-ifelse(!is.na(sIndex), relevantVar[sIndex]<relevantVar[uIndex], TRUE)
        canUbeQ<-superSetR[qIndex] & rInSuperSetS & uInSuperSetQ &
          relevantVar[rIndex]<relevantVar[qIndex] & conditionVarS & highVarQ<relevantVar[qIndex]
      }
      if(!is.na(sIndex)){
        # Can S be R, U be S?
        rInSuperSetQ<-getSuperSetQ(paramsPerLead = paramsPerLead, alpha = alpha, rIndex = sIndex)[rIndex]
        uInSuperSetS<-getSuperSetS(paramsPerLead = paramsPerLead, alpha = alpha, alphaN=alphaN, rIndex = sIndex)[uIndex]
        conditionVarQ<-ifelse(!is.na(qIndex), relevantVar[qIndex]<relevantVar[uIndex], TRUE)
        canUbeS<-superSetR[sIndex] & rInSuperSetQ & uInSuperSetS &
          relevantVar[rIndex]<relevantVar[sIndex] & conditionVarQ & highVarS<relevantVar[sIndex]
      }
    }

    # Register changes related to R wave
    if(canSbeR){
      rIndex<-sIndex; rChangesCount<-rChangesCount+1; rChanges<-c(rChanges,"S->R")
    }else if(canUbeS){
      rIndex<-sIndex; rChangesCount<-rChangesCount+1; rChanges<-c(rChanges,"U->S")
    }else if(canQbeR){
      rIndex<-qIndex; rChangesCount<-rChangesCount+1; rChanges<-c(rChanges,"Q->R")
    }else if(canUbeQ){
      rIndex<-qIndex; rChangesCount<-rChangesCount+1; rChanges<-c(rChanges,"U->Q")
    }else if(canPbeR){
      rIndex<-pIndex; rChangesCount<-rChangesCount+1; rChanges<-c(rChanges,"P->R")
    }else{
      rIndexIsValid<-TRUE
    }
  }

  #### U Wave Check ####
  if(!is.na(uIndex)){
    u1WaveCheck<-uWaveCheck(paramsPerLead=paramsPerLead, wavesNames=wavesNames, alpha=alpha,
                            alphaN=alphaN, uIndex=uIndex, qIndex=qIndex, rIndex=rIndex, sIndex=sIndex)
    if(u1WaveCheck$uChanges!=""){
      wavesNames<-u1WaveCheck$wavesNames; uChanges<-c(uChanges,u1WaveCheck$uChanges)
      paramsPerLead<-u1WaveCheck$paramsPerLead
    }else{
      # Sometimes, it may be needed to search for more Us
      continueU<-FALSE
      unusedCompU<-unusedComp
      while(!continueU){
        if(any(is.na(unusedCompU))){
          unusedCompU<-uIndex
        }else{
          unusedCompU<-c(unusedCompU, uIndex)
        }
        uIndex<-getIndexU(paramsPerLead = paramsPerLead, freeWaves = freeWaves,
                          unusedComp = unusedCompU)
        if(!is.na(uIndex)){
          uWaveCheck<-uWaveCheck(paramsPerLead=paramsPerLead, wavesNames=wavesNames, alpha=alpha,
                                 alphaN=alphaN, uIndex=uIndex, qIndex=qIndex, rIndex=rIndex, sIndex=sIndex)
          if(uWaveCheck$uChanges!=""){
            wavesNames<-uWaveCheck$wavesNames; uChanges<-c(uChanges,uWaveCheck$uChanges)
            paramsPerLead<-uWaveCheck$paramsPerLead
            continueU<-TRUE
          }
        }else{
          continueU<-TRUE
        }
      }
    }
  }

  # Waves are named accordingly and rChanges and uChanges are returned to global environment
  for(i in 1:nSignals) rownames(paramsPerLead[[i]])<-wavesNames
  rChanges<<-rChanges; uChanges<<-uChanges

  return(paramsPerLead)
}

#### Functions to evaluate and find better P waves candidates ####
doubtfulP<-function(paramsPerLead, alphaN, strict=TRUE, inSetPlusPlusP=FALSE){
  nWaves<-nrow(paramsPerLead[[1]])

  if("P" %in% rownames(paramsPerLead[[1]])){
    if(which(rownames(paramsPerLead[[1]])=="P")[1]!=7){
      rWave<-rownames(paramsPerLead[[1]])=="R"; pWave<-rownames(paramsPerLead[[1]])=="P"
      distancePR<-1-cos(paramsPerLead[[1]]$Alpha[rWave]-paramsPerLead[[1]]$Alpha[pWave])

      # Check omega conditions
      doubtOmegaCondition<-paramsPerLead[[1]]$Omega[pWave]<typMinOmegaP |
        paramsPerLead[[1]]$Omega[pWave]>typMaxOmegaP

      if(inSetPlusPlusP){
        # If needed, check if P is in setPlusPlusP
        isInSetPlusPlusP<-ifelse(!inSetPlusPlusP, TRUE,
                                 getSetPlusP(paramsPerLead = paramsPerLead, alphaN = alphaN,
                                             setPlusPlusP = TRUE, freeWaves = 1:nWaves)[pWave])

        return(doubtOmegaCondition | !isInSetPlusPlusP)
      }else{
        # If needed, check distance conditions
        if(strict){
          doubtDistCondition<-distancePR<typMinDistPR[1] | distancePR>maxDistPR[1]
        }else{
          doubtDistCondition<-distancePR<typMinDistPR[2] | distancePR>maxDistPR[1]
        }
        return(doubtDistCondition | doubtOmegaCondition)
      }
    }else{
      # If P is the seventh component, no need to search more
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}

# Function to determine if the new found wave J in SetPlusPlusP should substitute P
shouldDiscardPJ<-function(paramsPerLead, jIndexInSetPlusPlusP, alreadyExplainedVarWithoutP){
  pIndex<-which(rownames(paramsPerLead[[1]])=="P")
  availableRelevantLeads<-relevantLeads[relevantLeads %in% names(paramsPerLead)]

  # Extract measures
  relevantVarMatrixJ<-sapply(availableRelevantLeads,function(x) paramsPerLead[[x]][jIndexInSetPlusPlusP,"Var"])
  relevantVarMatrixP<-sapply(availableRelevantLeads,function(x) paramsPerLead[[x]][pIndex,"Var"])
  varCocientJP<-relevantVarMatrixJ/relevantVarMatrixP
  extraVarJ<-(relevantVarMatrixJ/(1-alreadyExplainedVarWithoutP))

  return(any(varCocientJP>maxRelativeVarDecrease) & any(extraVarJ>minExtraVarP) &
           (mean(relevantVarMatrixP)-mean(relevantVarMatrixJ))<maxVarDecrease)
}

# The new wave J' is worth it: should Q left as it is or should old P become Q?
shouldDiscardQP<-function(paramsPerLead, pNewIndex){
  pIndex<-which(rownames(paramsPerLead[[1]])=="P")
  qIndex<-ifelse("Q" %in% rownames(paramsPerLead[[1]]), which(rownames(paramsPerLead[[1]])=="Q"), NA)
  availableRelevantLeads<-relevantLeads[relevantLeads %in% names(paramsPerLead)]
  alphaCondition<-isAlpha1LeftOfAlpha2(paramsPerLead[[1]][pNewIndex,"Alpha"],
                                       paramsPerLead[[1]][pIndex,"Alpha"])

  if(is.na(qIndex)){
    return(alphaCondition)
  }else{
    # Extract measures
    relevantVarQ<-mean(sapply(availableRelevantLeads,function(x) paramsPerLead[[x]][qIndex,"Var"]))
    relevantVarP<-mean(sapply(availableRelevantLeads,function(x) paramsPerLead[[x]][pIndex,"Var"]))
    isOldPinSetQ<-pIndex %in% which(getSetQ(paramsPerLead=paramsPerLead))

    return(isOldPinSetQ & alphaCondition & relevantVarP> relevantVarQ)
  }
}

# Continue backfitting, even though R2 worsens, if the P found is more typical than previous
isImprovedP<-function(paramsPerLead, lastParamsPerLead, alphaN){
  currentIndexP<-which(rownames(paramsPerLead[[1]])=="P")
  lastIndexP<-which(rownames(lastParamsPerLead[[1]])=="P")
  pImproves<-FALSE

  # This should only be evaluated beginning on the second backfitting,
  # also, P must be assigned on both
  if(length(currentIndexP)+length(lastIndexP)>1){
    lastPinSetPlusPlusP<-getSetPlusP(paramsPerLead = lastParamsPerLead, alphaN = alphaN,
                                     setPlusPlusP = TRUE)[lastIndexP]
    currentPinSetPlusPlusP<-getSetPlusP(paramsPerLead = paramsPerLead, alphaN = alphaN,
                                        setPlusPlusP = TRUE)[currentIndexP]
    pImproves<-!lastPinSetPlusPlusP & currentPinSetPlusPlusP

    # If neither are in setPlusPlusP, check in setPlusP
    if(!lastPinSetPlusPlusP & !currentPinSetPlusPlusP){
      lastPinSetPlusP<-getSetPlusP(paramsPerLead = lastParamsPerLead, alphaN = alphaN,
                                   setPlusPlusP = FALSE)[lastIndexP]
      currentPinSetPlusP<-getSetPlusP(paramsPerLead = paramsPerLead, alphaN = alphaN,
                                      setPlusPlusP = FALSE)[currentIndexP]
      pImproves<-!lastPinSetPlusP & currentPinSetPlusP
    }
  }

  return(pImproves)
}

# In cases where the P maybe more distant than typical P waves
lastChanceP<-function(paramsPerLead, alreadyExplainedVar, alphaN){
  nWaves<-nrow(paramsPerLead[[1]]); wavesNames<-rownames(paramsPerLead[[1]])
  availableRelevantLeads<-relevantLeads[relevantLeads %in% names(paramsPerLead)]

  # Parameters extraction
  alpha<-paramsPerLead[[1]][,"Alpha"]
  relevantVarMatrix<-sapply(availableRelevantLeads,function(x) paramsPerLead[[x]]$Var)

  # P must be in SetPlusPlusP
  pCandidates<-getSetPlusP(paramsPerLead=paramsPerLead, setPlusPlusP=FALSE, alphaN = alphaN,
                           freeWaves=which(substr(rownames(paramsPerLead[[1]]),1,1)=="X"))

  # P must be left of Q, if assigned
  qIndex<-ifelse("Q" %in% rownames(paramsPerLead[[1]]), which(rownames(paramsPerLead[[1]])=="Q"), NA)
  if(!is.na(qIndex)){
    pCandidates <- pCandidates & isAlpha1LeftOfAlpha2(alpha1 = alpha, alpha2 = alpha[qIndex])
  }

  # P must explain extra variance
  explainsExtraVar<-sapply(1:nWaves, function(x){
    any((relevantVarMatrix[x,]/(1-alreadyExplainedVar))>minExtraVarP)})

  pCandidates<-pCandidates & explainsExtraVar

  # Assign if any fulfills the conditions
  if(any(pCandidates, na.rm = TRUE)){
    pIndex<-which(pCandidates)[1]; wavesNames[pIndex]<-"P"
    for(i in 1:length(paramsPerLead)) rownames(paramsPerLead[[i]])<-wavesNames
  }

  return(paramsPerLead)
}

#### Data loading functions ####
loadData<-function(dataPath="./01 Data/Preprocessed PTBXL Database/", leads=8){

  if(leads==8){
    cat("Loading leads... ")

    # Load 8-leads datafiles
    load(file=paste(dataPath,"detrend_I.RData", sep="")); cat("I ")
    load(file=paste(dataPath,"detrend_II.RData", sep="")); cat("II ")
    load(file=paste(dataPath,"detrend_V1.RData", sep="")); cat("V1 ")
    load(file=paste(dataPath,"detrend_V2.RData", sep="")); cat("V2 ")
    load(file=paste(dataPath,"detrend_V3.RData", sep="")); cat("V3 ")
    load(file=paste(dataPath,"detrend_V4.RData", sep="")); cat("V4 ")
    load(file=paste(dataPath,"detrend_V5.RData", sep="")); cat("V5 ")
    load(file=paste(dataPath,"detrend_V6.RData", sep="")); cat("V6 \n")

    # Annotations Matrix loading
    annotationsMatrix<-read.table(file = paste(dataPath,"AnnotationsMatrixV6.txt", sep=""), header = TRUE)
    annotationsMatrix[,6:13]<-sapply(c(6:13), function(x) as.logical(annotationsMatrix[,x]))

    # Pathology data
    classesLabels<-read.csv(file = paste(dataPath,"Complete_PTBXL_Classes.csv", sep=""), header = TRUE)

    return(list("dataPTBXL_I"=detrend_I, "dataPTBXL_II"=detrend_II,
                "dataPTBXL_V1"=detrend_V1, "dataPTBXL_V2"=detrend_V2, "dataPTBXL_V3"=detrend_V3,
                "dataPTBXL_V4"=detrend_V4, "dataPTBXL_V5"=detrend_V5, "dataPTBXL_V6"=detrend_V6,
                "annotationsMatrix"=annotationsMatrix, "classesLabels"=classesLabels))

  }else if(leads==12){

    cat("Loading leads... ")

    # Load all 12-leads datafiles
    load(file=paste(dataPath,"detrend_I.RData", sep="")); cat("I ")
    load(file=paste(dataPath,"detrend_II.RData", sep="")); cat("II ")
    load(file=paste(dataPath,"detrend_III.RData", sep="")); cat("III ")
    load(file=paste(dataPath,"detrend_AVR.RData", sep="")); cat("aVR ")
    load(file=paste(dataPath,"detrend_AVL.RData", sep="")); cat("aVL ")
    load(file=paste(dataPath,"detrend_AVF.RData", sep="")); cat("aVF ")
    load(file=paste(dataPath,"detrend_V1.RData", sep="")); cat("V1 ")
    load(file=paste(dataPath,"detrend_V2.RData", sep="")); cat("V2 ")
    load(file=paste(dataPath,"detrend_V3.RData", sep="")); cat("V3 ")
    load(file=paste(dataPath,"detrend_V4.RData", sep="")); cat("V4 ")
    load(file=paste(dataPath,"detrend_V5.RData", sep="")); cat("V5 ")
    load(file=paste(dataPath,"detrend_V6.RData", sep="")); cat("V6 \n")

    # Annotations Matrix loading
    annotationsMatrix<-read.table(file = paste(dataPath,"AnnotationsMatrixV6.txt", sep=""), header = TRUE)
    annotationsMatrix[,6:13]<-sapply(c(6:13), function(x) as.logical(annotationsMatrix[,x]))

    # Pathology data
    classesLabels<-read.csv(file = paste(dataPath,"Complete_PTBXL_Classes.csv", sep=""), header = TRUE)

    return(list("dataPTBXL_I"=detrend_I, "dataPTBXL_II"=detrend_II, "dataPTBXL_III"=detrend_III,
                "dataPTBXL_aVR"=detrend_AVR, "dataPTBXL_aVL"=detrend_AVL, "dataPTBXL_aVF"=detrend_AVF,
                "dataPTBXL_V1"=detrend_V1, "dataPTBXL_V2"=detrend_V2, "dataPTBXL_V3"=detrend_V3,
                "dataPTBXL_V4"=detrend_V4, "dataPTBXL_V5"=detrend_V5, "dataPTBXL_V6"=detrend_V6,
                "annotationsMatrix"=annotationsMatrix, "classesLabels"=classesLabels))

  }else{
    stop("Just 8 or 12 leads can be selected")
  }
}

readData<-function(loadedData, ecgId, beatId, usedLeads, otherLeadsCalculus=FALSE){

  annotationsMatrix<-loadedData$annotationsMatrix
  currentBeat<-annotationsMatrix[annotationsMatrix$ecg_id==ecgId & annotationsMatrix$beat_id==beatId,]
  eightLeads<-c("I", "II", "V1", "V2", "V3", "V4", "V5", "V6"); fourLeads<-c("I", "II", "V2", "V5")

  if(nrow(currentBeat)!=1){
    stop("Unanalyzed beat")
  }else if(apply(currentBeat[,eightLeads], 1, sum)[1]<2 & !otherLeadsCalculus){
    stop("Beat with single estimation lead")
  }else if(apply(currentBeat[,fourLeads], 1, sum)[1]==0 & !otherLeadsCalculus){
    stop("Beat without assignation leads")
  }else{
    segmentBegin<-currentBeat$iniRef; segmentEnd<-currentBeat$finRef
    segmentLength<-segmentEnd-segmentBegin
    annoActB<-currentBeat$annoRef-segmentBegin+1

    # Get signals from each lead
    iSignal<-NA; iiSignal<-NA; iiiSignal<-NA; avrSignal<-NA; avlSignal<-NA; avfSignal<-NA
    v1Signal<-NA; v2Signal<-NA; v3Signal<-NA; v4Signal<-NA; v5Signal<-NA; v6Signal<-NA
    if(!otherLeadsCalculus){
      if("I" %in% usedLeads & currentBeat$I) iSignal<-loadedData$dataPTBXL_I[ecgId,segmentBegin:segmentEnd]
      if("II" %in% usedLeads & currentBeat$II) iiSignal<-loadedData$dataPTBXL_II[ecgId,segmentBegin:segmentEnd]
      if("III" %in% usedLeads) iiiSignal<-loadedData$dataPTBXL_III[ecgId,segmentBegin:segmentEnd]
      if("aVR" %in% usedLeads) avrSignal<-loadedData$dataPTBXL_aVR[ecgId,segmentBegin:segmentEnd]
      if("aVL" %in% usedLeads) avlSignal<-loadedData$dataPTBXL_aVL[ecgId,segmentBegin:segmentEnd]
      if("aVF" %in% usedLeads) avfSignal<-loadedData$dataPTBXL_aVF[ecgId,segmentBegin:segmentEnd]
      if("V1" %in% usedLeads & currentBeat$V1) v1Signal<-loadedData$dataPTBXL_V1[ecgId,segmentBegin:segmentEnd]
      if("V2" %in% usedLeads & currentBeat$V2) v2Signal<-loadedData$dataPTBXL_V2[ecgId,segmentBegin:segmentEnd]
      if("V3" %in% usedLeads & currentBeat$V3) v3Signal<-loadedData$dataPTBXL_V3[ecgId,segmentBegin:segmentEnd]
      if("V4" %in% usedLeads & currentBeat$V4) v4Signal<-loadedData$dataPTBXL_V4[ecgId,segmentBegin:segmentEnd]
      if("V5" %in% usedLeads & currentBeat$V5) v5Signal<-loadedData$dataPTBXL_V5[ecgId,segmentBegin:segmentEnd]
      if("V6" %in% usedLeads & currentBeat$V6) v6Signal<-loadedData$dataPTBXL_V6[ecgId,segmentBegin:segmentEnd]
    }else{
      iSignal<-loadedData$dataPTBXL_I[ecgId,segmentBegin:segmentEnd]
      iiSignal<-loadedData$dataPTBXL_II[ecgId,segmentBegin:segmentEnd]
      iiiSignal<-loadedData$dataPTBXL_III[ecgId,segmentBegin:segmentEnd]
      avrSignal<-loadedData$dataPTBXL_aVR[ecgId,segmentBegin:segmentEnd]
      avlSignal<-loadedData$dataPTBXL_aVL[ecgId,segmentBegin:segmentEnd]
      avfSignal<-loadedData$dataPTBXL_aVF[ecgId,segmentBegin:segmentEnd]
      iiiSignal<-loadedData$dataPTBXL_III[ecgId,segmentBegin:segmentEnd]
      avrSignal<-loadedData$dataPTBXL_aVR[ecgId,segmentBegin:segmentEnd]
      avlSignal<-loadedData$dataPTBXL_aVL[ecgId,segmentBegin:segmentEnd]
      avfSignal<-loadedData$dataPTBXL_aVF[ecgId,segmentBegin:segmentEnd]
      v1Signal<-loadedData$dataPTBXL_V1[ecgId,segmentBegin:segmentEnd]
      v2Signal<-loadedData$dataPTBXL_V2[ecgId,segmentBegin:segmentEnd]
      v3Signal<-loadedData$dataPTBXL_V3[ecgId,segmentBegin:segmentEnd]
      v4Signal<-loadedData$dataPTBXL_V4[ecgId,segmentBegin:segmentEnd]
      v5Signal<-loadedData$dataPTBXL_V5[ecgId,segmentBegin:segmentEnd]
      v6Signal<-loadedData$dataPTBXL_V6[ecgId,segmentBegin:segmentEnd]
    }

    if(annoActB>segmentLength|annoActB<1){
      annoActB<-round(segmentLength*0.4)
      warning("Annotation is not correct, it has been put in the 40% of the signal")
    }

    # Get annotation in [0,2pi]
    alphaN<-getAlphaN(annoActB=annoActB, nObs=segmentLength)

    # Get class and superclass
    currentPatient<-loadedData$classesLabels[loadedData$classesLabels$EcgId==ecgId,]
    superclass<-currentPatient$Superclass; class<-currentPatient$Class

    # Save extracted Data
    return(list("vDataMatrix"=cbind("I"=iSignal, "II"=iiSignal, "III"=iiiSignal,
                                    "aVR"=avrSignal, "aVL"=avlSignal, "aVF"=avfSignal,
                                    "V1"=v1Signal, "V2"=v2Signal, "V3"=v3Signal,
                                    "V4"=v4Signal, "V5"=v5Signal, "V6"=v6Signal),
                "alphaN"=alphaN, "annoActB"=annoActB, "superclass"=superclass, "class"=class))
  }
}

checkData<-function(vDataMatrix, otherLeadsCalculus=FALSE){
  eightLeadsOrder<-c("I", "II", "V1", "V2", "V3", "V4", "V5", "V6")
  twelveLeadsOrder<-c("I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6")

  if(!otherLeadsCalculus){
    # Is there any missing column?
    missingColumns<-eightLeadsOrder[!eightLeadsOrder %in% colnames(vDataMatrix)]
    nMissing<-length(missingColumns)

    if(nMissing>0){
      allColumnsNames<-c(colnames(vDataMatrix),missingColumns)
      for(i in 1:nMissing) vDataMatrix<-cbind(vDataMatrix,NA)
      colnames(vDataMatrix)<-allColumnsNames
    }

    # Order columns
    vDataMatrix<-vDataMatrix[,eightLeadsOrder]
  }else{
    # Is there any missing column?
    missingColumns<-twelveLeadsOrder[!twelveLeadsOrder %in% colnames(vDataMatrix)]
    nMissing<-length(missingColumns)

    if(nMissing>0){
      allColumnsNames<-c(colnames(vDataMatrix),missingColumns)
      for(i in 1:nMissing) vDataMatrix<-cbind(vDataMatrix,NA)
      colnames(vDataMatrix)<-allColumnsNames
    }

    # Order columns
    vDataMatrix<-vDataMatrix[,twelveLeadsOrder]
  }


  return(list("vDataMatrix"=vDataMatrix,
              "denseDataMatrix"=vDataMatrix[,!sapply(1:ncol(vDataMatrix),function(x)
                any(is.na(vDataMatrix[,x])))]))
}

getBeatsPerPatient<-function(loadedData){
  return(loadedData$annotationsMatrix %>% group_by(ecg_id) %>%
           summarise(nbeats=max(beat_id), .groups = 'drop'))
}

getAlphaN<-function(annoActB, nObs){

  segmentBegin<-1; segmentEnd<-nObs
  # Annotation in [0,2pi], used for the assignment
  alphaN<-match(segmentBegin+annoActB,segmentBegin:segmentEnd)/(segmentEnd-segmentBegin+1)*2*pi
  if(is.na(alphaN)) alphaN<-which.min(abs(segmentBegin:segmentEnd-(segmentBegin+annoActB)))/(segmentEnd-segmentBegin+1)*2*pi
  #Sometimes, the comparison fails due to Floating-point Arithmetic

  return((alphaN+pi)%%(2*pi))
}

#### Other functions related to output ####
plotUnavailableLead<-function(leadName=NA){
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  if(!is.na(leadName)) text(x = 0.5, y = 0.5, paste0("Lead ",leadName," unavailable"), cex = 1.6)
}

fmmConstructor<-function(arguments){
  nObs<-length(arguments$data)
  fittedValues<-generateFMM(M=arguments$M, A=arguments$A, alpha=arguments$alpha, beta=arguments$beta,
                            omega=arguments$omega, length.out = nObs, plot = FALSE)$y
   FMM(
    timePoints=FMM::seqTimes(nObs),
    data=arguments$data, summarizedData=arguments$data, nPeriods=1,
    fittedValues=fittedValues,
    M=arguments$M, A=arguments$A,
    alpha=arguments$alpha, beta=arguments$beta, omega=arguments$omega,
    SSE=1, R2=arguments$R2, nIter=1
  )
}

plotMultiFMM_ECG<-function(vDataMatrix, fittedWaves, currentBack, paramsPerLead,
                           filename=NA, leadNames=1:length(paramsPerLead), path="./03 Results",
                           plotToFile=TRUE, unassigned=FALSE, extra=FALSE){
  
  # print(paramsPerLead)
  # If some of the 8 leads are unavailable, it should be noted on the plot
  allEightLeads<-c("I", "II", "V1", "V2", "V3", "V4", "V5", "V6")
  nSignals<-length(allEightLeads); halfSignals<-round(nSignals/2); signalIndex<-0
  plotLayout <- matrix(c(1:(2*nSignals), rep((2*nSignals)+1, halfSignals)),
                       nrow = 5, ncol = halfSignals, byrow = TRUE)

  if(plotToFile){
    if(!is.na(filename)){
      # Check if subdirectory exists; if not, create it
      if(!dir.exists(file.path(".", "Plots"))){
        dir.create(file.path(".", "Plots"))
      }
      png(filename=paste(path,"/Plots/",filename,"_Back_",currentBack,".png",sep=""),
          width=150*nSignals, height=600, type="cairo")
    }else{
      stop("Filename must be provided if plotToFile=TRUE")
    }
  }
  layout(mat = plotLayout, heights = c(0.225,0.225,0.225,0.225,0.1))
  par(mar=c(0.1,0.1,0.1,0.1))

  fittedWaves<-sapply(1:length(paramsPerLead), simplify = FALSE, function(y){
    leadParams<-paramsPerLead[[y]]
    arguments<-list("data"=vDataMatrix[,y], "M"=as.numeric(leadParams$M[which(substr(rownames(paramsPerLead[[1]]),1,1)!="X")[1]]),
                    "A"=as.numeric(leadParams$A), "alpha"=as.numeric(leadParams$Alpha), "beta"=as.numeric(leadParams$Beta),
                    "omega"=as.numeric(leadParams$Omega), "R2"=as.numeric(leadParams$Var))

    matrix(unlist(extractWaves(fmmConstructor(arguments))), nrow = nrow(vDataMatrix))}
  )

  sapply(1:halfSignals,function(x){
    if(allEightLeads[x] %in% leadNames){
      signalIndex<<-signalIndex+1
      plotMultiFMM_Sum(vDatai=vDataMatrix[!is.na(vDataMatrix[signalIndex,]),signalIndex], fittedWaves = fittedWaves[[signalIndex]],
                       currentBackResults = paramsPerLead[[signalIndex]],
                       currentBack=currentBack, filename=filename, leadName=leadNames[signalIndex],
                       unassigned=unassigned, extra=extra, plotToFile = FALSE)
    }else{plotUnavailableLead(leadName = allEightLeads[x])}
  })
  signalIndex<-0
  sapply(1:halfSignals,function(y){
    if(allEightLeads[y] %in% leadNames){
      signalIndex<<-signalIndex+1
      plotMultiFMM_Comps(vDatai=vDataMatrix[!is.na(vDataMatrix[signalIndex,]),signalIndex], fittedWaves = fittedWaves[[signalIndex]],
                         currentBackResults = paramsPerLead[[signalIndex]],
                         currentBack=currentBack, filename=filename, leadName=leadNames[signalIndex],
                         plotLegend=FALSE, unassigned=unassigned, plotToFile = FALSE)
    }else{plotUnavailableLead()}
  })
  maxHalfSignals<-signalIndex
  sapply((halfSignals+1):nSignals,function(x){
    if(allEightLeads[x] %in% leadNames){
      signalIndex<<-signalIndex+1
      plotMultiFMM_Sum(vDatai=vDataMatrix[!is.na(vDataMatrix[signalIndex,]),signalIndex], fittedWaves = fittedWaves[[signalIndex]],
                       currentBackResults = paramsPerLead[[signalIndex]],
                       currentBack=currentBack, filename=filename, leadName=leadNames[signalIndex],
                       unassigned=unassigned, extra=extra, plotToFile = FALSE)
    }else{plotUnavailableLead(leadName = allEightLeads[x])}
  })
  signalIndex<-maxHalfSignals
  sapply((halfSignals+1):nSignals,function(y){
    if(allEightLeads[y] %in% leadNames){
      signalIndex<<-signalIndex+1
      plotMultiFMM_Comps(vDatai=vDataMatrix[!is.na(vDataMatrix[signalIndex,]),signalIndex], fittedWaves = fittedWaves[[signalIndex]],
                         currentBackResults = paramsPerLead[[signalIndex]],
                         currentBack=currentBack, filename=filename, leadName=leadNames[signalIndex],
                         plotLegend=FALSE, unassigned=unassigned, plotToFile = FALSE)
    }else{plotUnavailableLead()}
  })

  # Add legend to the plot
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  if(unassigned){
    usedColors<-c(brewer.pal(n = 9, name = "Set1"), "aquamarine", "saddlebrown", "darkolivegreen1")
    usedColors<-usedColors[1:nrow(paramsPerLead[[1]])]
    names(usedColors)<-1:nrow(paramsPerLead[[1]])
  }else{
    usedColors<-getUsedColors(paramsPerLead[[1]])
    # Correct order
    if("X10" %in% names(usedColors)){
      x10Index<-which(names(usedColors)=="X10")
      usedColors<-c(usedColors[-x10Index],usedColors[x10Index])
    }
  }

  legend(x = "top",inset = 0, legend = names(usedColors), col=usedColors,
         lwd=5, cex=1.2, horiz = TRUE)

  if(plotToFile) dev.off()


}

prepareResults<-function(paramsPerLead, leadNames,
                         ## Arguments to be suppressed in app version
                         ecgId=NA, beatId=NA){
  singleResultNames<-c("M","A","Alpha","Beta","Omega","Var")
  allWavesNames<-c("P","Q","R","S","T")
  nSignals<-length(paramsPerLead)

  allLeadsNames<-c("I", "II", "V1", "V2", "V3", "V4", "V5", "V6")
  names(paramsPerLead)<-leadNames
  availableRelevantLeads<-relevantLeads[relevantLeads %in% names(paramsPerLead)]
  wavesNames<-rownames(paramsPerLead[[1]])
  nSignals<-length(leadNames)

  # The most important unused wave must be saved
  uIndex<-getIndexU(paramsPerLead=paramsPerLead, freeWaves = which(substr(rownames(paramsPerLead[[1]]),1,1)=="X"))
  alphaU<-ifelse(!is.na(uIndex),paramsPerLead[[1]]$Alpha[uIndex],NA)
  omegaU<-ifelse(!is.na(uIndex),paramsPerLead[[1]]$Omega[uIndex],NA)
  relVarU<-ifelse(!is.na(uIndex),mean(sapply(availableRelevantLeads, function(x) paramsPerLead[[x]]$Var[uIndex]), na.rm=TRUE),NA)
  meanVarU<-ifelse(!is.na(uIndex),mean(sapply(1:nSignals, function(x) paramsPerLead[[x]]$Var[uIndex]), na.rm=TRUE),NA)

  # There may be unassigned waves in the results
  if(any(substr(rownames(paramsPerLead[[1]]),1,1)=="X")){

    # Unassigned waves are deleted
    unassignedWaves<-substr(rownames(paramsPerLead[[1]]),1,1)=="X"
    paramsPerLead<-lapply(1:nSignals, function(x) paramsPerLead[[x]][!unassignedWaves,])

    # If some waves have not been observed
    if(!all(allWavesNames %in% rownames(paramsPerLead[[1]]))){
      # Unobserved waves parameters are added as empty rows
      unobservedWaves<-allWavesNames[!allWavesNames %in% rownames(paramsPerLead[[1]])]
      paramsPerLead<-lapply(1:nSignals, function(x){
        for(j in unobservedWaves){
          emptyRow<-data.frame(matrix(ncol = length(singleResultNames), dimnames = list(j,singleResultNames)))
          paramsPerLead[[x]]<-rbind(paramsPerLead[[x]], emptyRow)
        }
        return(paramsPerLead[[x]])
      })
    }
  }

  # Some leads may not have not been used
  if(nSignals<8){
    unusedLeads<-allLeadsNames[!allLeadsNames %in% leadNames]
    for(i in unusedLeads){
      nSignals<-nSignals+1
      paramsPerLead[[nSignals]]<-as.data.frame(matrix(NA, nrow = nrow(paramsPerLead[[1]]), ncol = ncol(paramsPerLead[[1]]),
                                               dimnames = dimnames(paramsPerLead[[1]])))
      leadNames<-c(leadNames, i)
    }
  }

  # Reorder rows of paramsPerLead
  paramsPerLead<-lapply(1:nSignals, function(x) paramsPerLead[[x]][order(rownames(paramsPerLead[[x]])),])
  wavesNames<-rownames(paramsPerLead[[1]])

  nBack<-nrow(paramsPerLead[[1]])
  matrixNames<-list(leadNames, paste(rep(singleResultNames,each=length(wavesNames)),
                                     rep(wavesNames,length(singleResultNames)), sep=""))
  result <- data.frame(matrix(unlist(paramsPerLead), nrow=length(paramsPerLead),
                              byrow=TRUE, dimnames = matrixNames))
  # M parameter is repeated
  result <- result[,-c(1:2,4:nBack)]; colnames(result)[1]<-"M"

  ## Calculate R2
  result$R2<-sapply(1:nSignals, function(x){
    if(all(is.na(paramsPerLead[[x]]$Var))){
      return(NA)
    }else{
      sum(paramsPerLead[[x]]$Var, na.rm = TRUE)
    }
  })

  ## Transformed from wide to long
  resultLong<-result %>%rownames_to_column() %>%
    pivot_longer(cols = -rowname) %>%
    unite(name, name, rowname, sep = "_")

  ## Drop repeated measurements
  dropAlphaIndexes<-which(substr(resultLong$name, 1,5)=="Alpha")[-c(1:nBack)]
  dropOmegaIndexes<-which(substr(resultLong$name, 1,5)=="Omega")[-c(1:nBack)]
  resultLong<-resultLong[-c(dropAlphaIndexes, dropOmegaIndexes),]

  ## Alpha and Omega should be before
  alphaIndexes<-which(substr(resultLong$name, 1,5)=="Alpha")
  omegaIndexes<-which(substr(resultLong$name, 1,5)=="Omega")
  resultLong<-resultLong[c(alphaIndexes,omegaIndexes,
                           (1:nrow(resultLong))[-c(alphaIndexes,omegaIndexes)]),]
  resultLong$name[1:(2*nBack)]<-substr(resultLong$name[1:(2*nBack)],1,6)

  ## Return to wide dataframe
  resultWide<-resultLong %>% pivot_wider()
  resultWide$EcgId<-ecgId; resultWide$BeatId<-beatId
  resultWide<-resultWide[,c(c(-1,0)+ncol(resultWide),1:(ncol(resultWide)-2))]

  # Add U wave parameters
  resultWide<-cbind(resultWide, "AlphaU"=alphaU, "OmegaU"=omegaU, "RelVarU"=relVarU, "MeanVarU"=meanVarU)

  ## Columns must be correctly ordered to avoid incongruities
  colOrder<-c("EcgId","BeatId", paste0(rep(c("Alpha","Omega"), each=5),allWavesNames),
              paste(rep(c("M",paste0(rep(c("A","Beta","Var"), each=5),allWavesNames), "R2"), nSignals),
                    rep(allLeadsNames,each=17),sep="_"),"AlphaU","OmegaU","RelVarU","MeanVarU")
  resultWide<-resultWide[,colOrder]

  ## Results can be saved
  return(resultWide)
}

renameUnusedBackplots<-function(filename, totalBack, unusedBackfittings,
                                resultsPath="./03 Results/Plots/"){
  for(i in unusedBackfittings){
    lastFilename<-paste0(resultsPath,filename,"_Back_",i)
    newName<-paste0(lastFilename,"_UNUSED.png"); lastFilename<-paste0(lastFilename,".png")
    file.rename(lastFilename, newName)
  }
}

moveBackplots<-function(filename, totalBack, unusedBackfittings=NA,
                        resultsPath="./03 Results/Plots/"){
  backfittingPath<-paste(resultsPath,"Supplementary/Backfitting/",sep="")

  # Check if subdirectories exist
  if(!dir.exists(file.path(resultsPath, "Supplementary"))){
    dir.create(file.path(resultsPath, "Supplementary"))
  }
  if(!dir.exists(file.path(paste(resultsPath,"Supplementary/",sep=""), "Backfitting"))){
    dir.create(file.path(paste(resultsPath,"Supplementary/",sep=""), "Backfitting"))
  }

  # Last used backfitting plot: in both the backfitting directory and the parent Plots directory
  fromFilename<-paste(resultsPath,filename,"_Back_",totalBack,".png",sep="")
  toFilename<-paste(backfittingPath,filename,"_Back_",totalBack,".png",sep="")
  file.copy(from = fromFilename, to = toFilename)

  # The rest of the backfittings are moved to the Backfitting subdirectory
  if(totalBack>1){
    for(i in 1:(totalBack-1)){
      fromFilename<-paste(resultsPath,filename,"_Back_",i,".png",sep="")
      toFilename<-paste(backfittingPath,filename,"_Back_",i,".png",sep="")
      file.rename(from = fromFilename, to = toFilename)
    }
  }

  # The unused backfitting is also stored in the Backfitting subdirectory
  if(any(!is.na(unusedBackfittings))){
    for(i in unusedBackfittings){
      fromFilename<-paste(resultsPath,filename,"_Back_",i,"_UNUSED.png",sep="")
      toFilename<-paste(backfittingPath,filename,"_Back_",i,"_UNUSED.png",sep="")
      file.rename(from = fromFilename, to = toFilename)
    }
  }
}


# Function to save the results in the appropriate CSV.
saveResultsCSV<-function(ecgId, beatId, leadNames, paramsPerLead, superclass=NA,
                         class=NA, totalBack=NA, stopCriteria=NA, rChanges="",uChanges="",
                         annotation=NA, alphaN=NA, nObs=NA, errorMessage=NA,
                         path="./03 Results", resultsFile="multiLeadResults.csv"){

  # Results are prepared only if they are not an error
  if(is.na(errorMessage)){
    preparedResults<-prepareResults(paramsPerLead=paramsPerLead, leadNames=leadNames,
                                    ecgId=ecgId, beatId=beatId)
  }else{
    preparedResults<-paramsPerLead
  }

  # Set classes in the correct position
  preparedResults$Superclass<-superclass; preparedResults$Class<-class
  nCurrentColumns<-ncol(preparedResults)
  preparedResults<-preparedResults[,c(1:2,(nCurrentColumns-1):nCurrentColumns,
                                      3:(nCurrentColumns-2))]

  # Set annotations and fitting process columns
  preparedResults$nBack<-totalBack; preparedResults$StopCriteria<-stopCriteria
  preparedResults$Annotation<-annotation; preparedResults$AlphaN<-alphaN
  preparedResults$nObs<-nObs; preparedResults$RChanges<-rChanges; preparedResults$UChanges<-uChanges
  preparedResults$Error<-errorMessage

  destFile<-paste(path,resultsFile,sep="/")
  if(file.exists(destFile)){
    write.table(preparedResults, file = destFile,sep = ",",
                append = TRUE, quote = FALSE,col.names = FALSE,
                row.names = FALSE)
  } else{
    write.table(preparedResults, file = destFile,sep = ",",
                append = TRUE, quote = FALSE,col.names = TRUE,
                row.names = FALSE)

  }
}

# Function to add other derived measures like explained vars, distances,...
addDerivedMeasures<-function(completeResults,
                             loadedData=loadData(dataPath = "../../01 Data/Preprocessed PTBXL Database/", leads = 8),
                             extraInfo=read.csv("../../01 Data/Preprocessed PTBXL Database/PTBXL_AdditionalInfoV2.csv")){

  # Add beatIdNew
  contBeat<-1
  completeResults$BeatIdNew<-sapply(1:nrow(completeResults), function(x){
    if(x>1){
      contBeat<<-ifelse(completeResults[x,"EcgId"]!=completeResults[x-1,"EcgId"] & !is.na(completeResults[x,"Superclass"]), 1,
                        # Patient change, first beat is valid
                        ifelse(completeResults[x,"EcgId"]!=completeResults[x-1,"EcgId"] & is.na(completeResults[x,"Superclass"]), 0,
                               # Patient change, first beat is not valid
                               ifelse(!is.na(completeResults[x,"Superclass"]), contBeat+1, contBeat)))
      return(ifelse(!is.na(completeResults[x,"Superclass"]), contBeat, NA))
    }else{
      return(contBeat)
    }
  })

  # Add MeanVar and RelVar variables
  twelveLeads<-c("I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6")
  eightLeads<-c("I", "II", "V1", "V2", "V3", "V4", "V5", "V6"); fourLeads<-c("I", "II", "V2", "V5")
  completeResults$MeanVarP<-apply(completeResults[paste0("VarP_",eightLeads)],1,mean,na.rm=TRUE); completeResults$RelVarP<-apply(completeResults[paste0("VarP_",fourLeads)],1,mean,na.rm=TRUE)
  completeResults$MeanVarQ<-apply(completeResults[paste0("VarQ_",eightLeads)],1,mean,na.rm=TRUE); completeResults$RelVarQ<-apply(completeResults[paste0("VarQ_",fourLeads)],1,mean,na.rm=TRUE)
  completeResults$MeanVarR<-apply(completeResults[paste0("VarR_",eightLeads)],1,mean,na.rm=TRUE); completeResults$RelVarR<-apply(completeResults[paste0("VarR_",fourLeads)],1,mean,na.rm=TRUE)
  completeResults$MeanVarS<-apply(completeResults[paste0("VarS_",eightLeads)],1,mean,na.rm=TRUE); completeResults$RelVarS<-apply(completeResults[paste0("VarS_",fourLeads)],1,mean,na.rm=TRUE)
  completeResults$MeanVarT<-apply(completeResults[paste0("VarT_",eightLeads)],1,mean,na.rm=TRUE); completeResults$RelVarT<-apply(completeResults[paste0("VarT_",fourLeads)],1,mean,na.rm=TRUE)

  # Add distance to R variables
  completeResults$dRP<-1-cos(completeResults$AlphaP-completeResults$AlphaR)
  completeResults$dRQ<-1-cos(completeResults$AlphaQ-completeResults$AlphaR)
  completeResults$dRS<-1-cos(completeResults$AlphaS-completeResults$AlphaR)
  completeResults$dRT<-1-cos(completeResults$AlphaT-completeResults$AlphaR)

  # Add RPrime measures
  cat("Calculating RPrime measure...\n")
  diffSquares<-t(sapply(1:nrow(completeResults), function(x){
    ecgId<-completeResults[x, "EcgId"]; beatId<-completeResults[x, "BeatId"]
    if(is.na(completeResults[x, "Error"])){
      vDataMatrix<-readData(loadedData = loadedData, ecgId = ecgId, beatId = beatId,
                            usedLeads = eightLeads)$vDataMatrix
      meanMatrix<-matrix(rep(apply(vDataMatrix, 2, mean), each=nrow(vDataMatrix)), nrow = nrow(vDataMatrix), ncol = ncol(vDataMatrix))
      return(t(apply(((vDataMatrix-meanMatrix)^2), 2, sum)/apply((vDataMatrix^2), 2, sum)))
    }else{
      return(rep(NA,12))
    }
  })); colnames(diffSquares)<-twelveLeads
  completeResults[paste0("RPrime_",eightLeads)]<-1-((1-completeResults[,paste0("R2_",eightLeads)])*diffSquares[,eightLeads])

  # Order columns
  resultsPerLead<-NULL; for(i in 0:7){
    resultsPerLead<-c(resultsPerLead, c((14:30)+((17*i)+1),ncol(completeResults)-(7-i)))
  }
  colOrder<-c(1:2,ncol(completeResults)-22,3:14,
              seq(ncol(completeResults)-20,ncol(completeResults)-12,2),
              seq(ncol(completeResults)-21,ncol(completeResults)-12,2),
              (ncol(completeResults)-3):ncol(completeResults)-8,
              resultsPerLead, (ncol(completeResults)-34):(ncol(completeResults)-23))
  completeResults<-completeResults[,colOrder]

  # Order rows by EcgId and BeatId
  completeResults<-completeResults[with(completeResults, order(EcgId, BeatId)),]

  # Add extra information
  completeResults$Class[completeResults$Class=="AF"]<-"OTRA"
  othersClass<-sapply(completeResults$EcgId, function(x) extraInfo[extraInfo$EcgId==x, "OthersLabel"])
  age<-sapply(completeResults$EcgId, function(x) extraInfo[extraInfo$EcgId==x, "Age"])
  sex<-sapply(completeResults$EcgId, function(x) extraInfo[extraInfo$EcgId==x, "Sex"])
  isPace<-sapply(completeResults$EcgId, function(x) extraInfo[extraInfo$EcgId==x, "PACE"])
  likelihood<-sapply(completeResults$EcgId, function(x) extraInfo[extraInfo$EcgId==x, "Likelihood"])
  isHighLikelihood<-sapply(completeResults$EcgId, function(x) extraInfo[extraInfo$EcgId==x, "IsLikelihoodHigh"])
  completeResults<-cbind(completeResults, "OthersClass"=othersClass, "Age"=age, "Sex"=sex, "PACE"=isPace,
                         "IsLikelihoodHigh"=isHighLikelihood,"Likelihood"=likelihood)

  # Reorder columns
  completeResults<-completeResults[,c(1:5, ncol(completeResults)-c(5,2,4,3,0,1),
                                      6:(ncol(completeResults)-6))]

  return(completeResults)
}

# Function to load current results
loadResults<-function(resultsPath="../03 Results/", resultsFile="completeMultiLeadResults.csv"){
  cat("Loading results...")
  results<-read.csv(paste0(resultsPath,resultsFile))
  cat("\n")
  return(results)
}

# Function to simulate multiFMM signals
generateMultiFMM<-function(alpha, omega, mVector, aMatrix, betaMatrix,
                           colMatrixNames=1:ncol(aMatrix),
                           from = 0, to = 2*pi, length.out = 200){
  timePoints = seq(from = from, to = to, length.out = length.out)
  nSignals<-ncol(aMatrix)
  vDataMatrix<-matrix(nrow = length.out, ncol = nSignals)

  for(i in 1:nSignals){
    vDataMatrix[,i] <- as.vector(mVector[i] + aMatrix[,i]%*%
                               t( calculateCosPhi(alpha = alpha, beta = betaMatrix[,i],
                                                 omega = omega, timePoints = timePoints)))
  }
  colnames(vDataMatrix)<-colMatrixNames

  par(mfrow=c(4,3))
  dummyVar<-sapply(1:nSignals, function(x){plot(vDataMatrix[,x], type="l", ylab = NA, xlab = NA, main = colMatrixNames[x])})

  return(vDataMatrix)
}



splitData<-function(data=csvReaded,mAnnot=annots,idBeat=beatId){
  nLeads<-twelveLeads<-c("I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6")
  a<-mAnnot[idBeat,2]
  b<-mAnnot[idBeat,3]
  d<-mAnnot[idBeat,1]
  if(length(nLeads)==12){
    vDataBeat<-data[a:b,nLeads]
    #plot(data[mAnnot[idBeat,2]:mAnnot[idBeat,3],nLeads][,2],type="l")
  }else{
    vDataBeatAux<-data[a:b,nLeads]
    putNA<-rep(NA,nrow(vDataBeatAux))
    vDataBeat<-cbind(vDataBeatAux[,1:2],putNA,putNA,putNA,putNA,vDataBeatAux[,3:8])
    colnames(vDataBeat)<-twelveLeads
  }
  alphaN<-(d-a)/(b-a+1)*(2*pi)
  return(list("beatMatrix"=vDataBeat ,"alphaN"=alphaN, "annoActB"=d-a+1))
}

giveSegmentsECG<-function(allECG,allRPeaks,fr){
  sucBeats<-giveSucesiveECG(allRPeaks,fr)[[1]]
  medRR<-giveSucesiveECG(allRPeaks,fr)[[2]]
  if(!is.na(medRR)){
    datpac<-allECG
    boundBeats<-matrix(0,length(allRPeaks),2)
    conti<-c()
    ya<-FALSE
    for(j in 1:length(allRPeaks)){
      
      if(sum(sucBeats[j,]==c(1,0,1))==3){
        conti[j]<-"R&L Cont"
      }else{
        if(sum(sucBeats[j,]==c(0,0,1))==3){
          conti[j]<-"R Cont"
        }else{
          if(sum(sucBeats[j,]==c(1,0,0))==3){
            conti[j]<-"L Cont"
          }else{
            if(sum(sucBeats[j,]==c(0,1,0))==3){
              conti[j]<-"Alone"
            }
          }
        }
      }
      
      if(sum(sucBeats[j,]==c(0,0,1))==3){
        RRbeat<-allRPeaks[j+1]-allRPeaks[j]
        if(sum(sucBeats[j+1,]==c(0,1,0))!=3){
          if(j==1 & floor(allRPeaks[j]-0.4*RRbeat)<(1:length(datpac))[1]){
            boundBeats[j,]<-c((1:length(datpac))[1],ceiling(allRPeaks[j]+0.6*RRbeat))
          }else{
            boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*RRbeat),ceiling(allRPeaks[j]+0.6*RRbeat))
          }
        }else{
          boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*medRR),ceiling(allRPeaks[j]+0.6*medRR))
        }
      }
      if(sum(sucBeats[j,]==c(1,0,0))==3 | (sum(sucBeats[j,]==c(0,1,0))==3 & j==length(allRPeaks))){
        RRbeat<-allRPeaks[j]-allRPeaks[j-1]
        if(sum(sucBeats[j-1,]==c(0,1,0))!=3 & (sum(sucBeats[j,]==c(0,1,0))!=3)){
          if(j==length(allRPeaks) & ceiling(allRPeaks[j]+0.6*RRbeat)>(1:length(datpac))[length((1:length(datpac)))]){
            boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*RRbeat),(1:length(datpac))[length((1:length(datpac)))])
            ya<-TRUE
          }else{
            boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*RRbeat),ceiling(allRPeaks[j]+0.6*RRbeat))
          }
        }else{
          boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*medRR),ceiling(allRPeaks[j]+0.6*medRR))
        }
      }
      if(sum(sucBeats[j,]==c(1,0,1))==3){
        RRbeatL<-allRPeaks[j]-allRPeaks[j-1]
        RRbeatR<-allRPeaks[j+1]-allRPeaks[j]
        boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*RRbeatL),
                          ceiling(allRPeaks[j]+0.6*RRbeatR))
      }
      if(sum(sucBeats[j,]==c(0,1,0))==3 & ya==FALSE){
        RRbeat<-medRR
        boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*medRR),ceiling(allRPeaks[j]+0.6*medRR))
      }
    }
    #corregimos para no ajustar dos veces los mismos puntos
    for(j in 2:nrow(boundBeats)){
      if(boundBeats[j,1]<=boundBeats[j-1,2])boundBeats[j,1]<-boundBeats[j-1,2]+1
    }
  }else{
    boundBeats<-c()
  }
  
  return(boundBeats)
}

giveSucesiveECG<-function(vAnnotN,frec){
  if(length(vAnnotN)>1){
    lRR<-c()
    for(i in 2:length(vAnnotN)){
      lRR[i-1]<-vAnnotN[i]-vAnnotN[i-1]
    }
    #si en al menos el 75% de ellos la longitud es mayor que la frequencia y tenemos mas de 3 latinos	
    if(sum(lRR>frec)>0.75*length(lRR) & length(lRR)>3){
      lRR<-sort(lRR)[-c(length(lRR):(length(lRR)-2))]
    }else{
      #si tenemos mas de tres latidos y en al menos uno la long es mayor que la frequencia 
      if(sum(lRR>frec)>0 & length(lRR)>3){
        lRR<-(lRR)[-which(lRR>frec)]##
      }else{
        #si tenemos tres latidos y  en alguno la longitud mayor a la frequencia
        if(sum(lRR>frec)>0 & length(lRR)==3){
          lRR<-sort(lRR)[-length(lRR)]
        }
      }
    }
    
    M<-median(lRR)
    #print(paste("es m",M))
    fin<-c()
    for(i in 1:length(vAnnotN)){
      if(i==1 &
         (vAnnotN[2]-vAnnotN[1])<1.5*M)idSuc<-c(0,0,1)#podr?a se problema si RR_12 != RR_23, asumimos que no ocurre
      if(i==1 &
         (vAnnotN[2]-vAnnotN[1])>=1.5*M)idSuc<-c(0,1,0)
      if(i==length(vAnnotN) &
         (vAnnotN[length(vAnnotN)]-vAnnotN[length(vAnnotN)-1])<1.5*M)idSuc<-c(1,0,0)
      if(i==length(vAnnotN) &
         (vAnnotN[length(vAnnotN)]-vAnnotN[length(vAnnotN)-1])>=1.5*M)idSuc<-c(0,1,0)
      if( i >=2 & i<=(length(vAnnotN)-1)){
        RRl<-vAnnotN[i]-vAnnotN[i-1]
        RRd<-vAnnotN[i+1]-vAnnotN[i]
        if(RRl<1.5*M & RRd<1.5*M)idSuc<-c(1,0,1)
        if(RRl<1.5*M & RRd>1.5*M)idSuc<-c(1,0,0)
        if(RRl>1.5*M & RRd<1.5*M)idSuc<-c(0,0,1)
        if(RRl>1.5*M & RRd>1.5*M)idSuc<-c(0,1,0)
      }
      fin<-rbind(fin,idSuc)
    }
  }else{
    M<-NA
    fin<-c()
  }
  
  return(list(fin,M))
}
