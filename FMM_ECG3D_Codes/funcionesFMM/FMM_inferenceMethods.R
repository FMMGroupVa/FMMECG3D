


################################################################################
# Log-likelihood function for a single wave
# Arguments:
#   parameters:  alpha, omega
#   vData: data to be fitted an FMM model.
#   varErr: estimated variance of the error (variance of residuals when omega
#           and alpha are ML).
# Returns: logLikelihood of FMM model (given an estimation of the variance)
################################################################################
logLikFMMWave <- function(alphaOmegaParameters, vData, varErr){
  alpha <- alphaOmegaParameters[1]%%(2*pi)
  omega <- max(min(alphaOmegaParameters[2], 1), 0.0001)
  fullParams <- determineNoiseParameters(alphaOmegaParameters = c(alpha, omega),
                                         vData = vData)
  return(-sum(0.5*(1/varErr)*(vData-fullParams[["yFit"]])^2))
}



################################################################################
# Calculate FMM1 parameters from (alpha, omega).
# Arguments:
#   alphaOmegaParameters:  vector (alpha, omega).
#   vData: data to be fitted an FMM model.
#   timePoints:
# Returns: List[ par = "FMM1 parameters",
#                RSS = Residual Sum of Squares,
#                yFit = fitted values ]
################################################################################
determineNoiseParameters <- function(alphaOmegaParameters, vData,
                                     timePoints = seqTimes(length(vData))){
  alpha <- alphaOmegaParameters[1]
  omega <- alphaOmegaParameters[2]
  tStar <- alpha + 2*atan2(omega*sin((timePoints-alpha)/2), cos((timePoints-alpha)/2))
  
  dM <- cbind(rep(1, length(vData)), cos(tStar), sin(tStar))
  mDeltaGamma <- solve(t(dM)%*%dM)%*%t(dM)%*%vData
  
  # Calculus of explicit parameters, RSS and fitted values
  optM <- mean(vData) - mDeltaGamma[2]*mean(cos(tStar)) - mDeltaGamma[3]*
    mean(sin(tStar))
  yFit <- optM + mDeltaGamma[2]*cos(tStar) + mDeltaGamma[3]*sin(tStar)
  optBeta <- atan2(-mDeltaGamma[3], mDeltaGamma[2]) + alpha
  optA <- sqrt(sum(mDeltaGamma[2:3]^2))
  optRSS <- sum((vData - yFit)^2)/length(vData)
  return(list(par = c(optM, optA, alpha%%(2*pi), optBeta%%(2*pi), omega),
              RSS = optRSS, yFit = yFit))
}


