################################################################################
# Internal function: fit monocomponent FMM model
# Arguments:
#   vData: data to be fitted an FMM model.
#   timePoints: one single period time points.
#   lengthAlphaGrid, lengthOmegaGrid: precision of the grid of
#                                     alpha and omega parameters.
#   alphaGrid, omegaGrid: grids of alpha and omega parameters.
#   omegaMax: max value for omega.
#   numReps: number of times the alpha-omega grid search is repeated.
# Returns an object of class FMM.
################################################################################
fitFMM_unit <- function(vData, timePoints = seqTimes(length(vData)),
                      lengthAlphaGrid = 48, lengthOmegaGrid = 18,
                      alphaGrid = seq(0, 2*pi, length.out = lengthAlphaGrid),
                      omegaMin = 0.0001, omegaMax = 1,
                      omegaGrid = exp(seq(log(0.01), log(omegaMax),
                                          length.out = lengthOmegaGrid)),
                      numReps = 3, usedApply = getApply(FALSE)[[1]]){
  nObs <- length(vData)
  grid <- expand.grid(alphaGrid, omegaGrid)
  step1OutputNames <- c("M","A","alpha","beta","omega","RSS")

  ## Step 1: initial values of M, A, alpha, beta and omega. Parameters alpha and
  # omega are initially fixed and cosinor model is used to calculate the rest of the parameters.
  # step1FMM function is used to make this estimate.
  step1 <- usedApply(FUN = step1FMM, X = grid, vData = vData,
                     timePoints = timePoints)
  colnames(step1) <- step1OutputNames

  # We use bestStep1 internal function
  bestPar <- bestStep1(vData, step1)
  postOpt <- optim(par = bestPar[c(3, 5)], fn = step2FMMExp, method = "L-BFGS-B",
                   lower = c(0, 0.001), upper = c(2*pi, 0.99), vData = vData)

  parametersRSSandFittedValues <- determineNoiseParameters(postOpt$par, vData)
  parFinal <- parametersRSSandFittedValues[["par"]]

  names(parFinal) <- step1OutputNames[-6]
  # Returns an object of class FMM.
  fittedFMMvalues <- parFinal["M"] + parFinal["A"]*cos(parFinal["beta"] +
        2*atan(parFinal["omega"]*tan((timePoints-parFinal["alpha"])/2)))
  SSE <- sum((fittedFMMvalues-vData)^2)

  return(FMM(
    M = parFinal[[1]],
    A = parFinal[[2]],
    alpha = parFinal[[3]],
    beta = parFinal[[4]],
    omega = parFinal[[5]],
    timePoints = timePoints,
    summarizedData = vData,
    fittedValues = fittedFMMvalues,
    SSE = SSE,
    R2 = PV(vData, fittedFMMvalues),
    nIter = 0))
}

