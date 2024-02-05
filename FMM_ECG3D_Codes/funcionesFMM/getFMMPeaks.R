#' Peak and trough times and signal values
#'
#'
#' \code{getFMMPeaks()} is used to estimate peak and trough times and signal values at those times for each component of the model. These parameters result to be useful in multiple applications.
#'
#'
#' @param objFMM Object of class \code{'FMM'}
#' @param timePointsIn2pi \code{TRUE} to return peak and trough times in the \eqn{[0, 2\pi]} interval. When \code{timePointsIn2pi = FALSE} the positions of peak and trough times are returned. Its default value is \code{TRUE}.
#'
#'
#' @return
#' A list with the following components is returned:
#' \item{tpeakU}{a numeric vector with the time points at which the peak of each wave is estimated.}
#' \item{tpeakL}{a numeric vector with the time points at which the trough of each wave is estimated.}
#' \item{ZU}{a numeric vector with the estimated signal peak values of each wave.}
#' \item{ZL}{a numeric vector with the estimated signal trough values of each wave.}
#'
#'
#' @references
#' Rueda C, Larriba Y, Peddada SD (2019).
#' Frequency Modulated Moebius Model Accurately Predicts Rhythmic Signals in Biological and Physical Sciences.
#' \emph{Scientific reports}, \bold{9} (1), 18701. \url{https://www.nature.com/articles/s41598-019-54569-1}
#'
#'
#' @examples
#' ## Generate example data:
#' fmm2.data <- generateFMM(0, rep(2, 2), c(1.5, 3.4), c(0.2, 2.3), c(0.1, 0.2),
#'                          plot = FALSE, outvalues = TRUE,
#'                          sigmaNoise = 0.5) # add a gaussian noise with sigma = 0.5
#'
#' ## Fit the FMM model with nback = 2 components
#' ## fit is an object of S4 class 'FMM'
#' fit <- fitFMM(fmm2.data$y,timePoints = fmm2.data$t,nback = 2,
#'               lengthAlphaGrid = 24,lengthOmegaGrid = 10)
## estimate peakes and trough times and signals
#' getFMMPeaks(fit, timePointsIn2pi = TRUE) # times in the [0,2*pi] interval
#'
getFMMPeaks <- function(objFMM, timePointsIn2pi = TRUE) {

  M <- getM(objFMM)
  A <- getA(objFMM)
  alpha <- getAlpha(objFMM)
  beta <- getBeta(objFMM)
  omega <- getOmega(objFMM)

  data<- getData(objFMM)
  if(getNPeriods(objFMM) > 1) data <- getSummarizedData(objFMM)

  nData <- length(data)
  nComp <- length(alpha)

  # timePoints estimation
  peakU <- (alpha + 2*atan2(1/omega*sin(-beta/2), cos(-beta/2))) %% (2*pi)
  peakL <- (alpha + 2*atan2(1/omega*sin((pi-beta)/2), cos((pi-beta)/2))) %% (2*pi)

  if(timePointsIn2pi){
    tpeakU <- peakU
    tpeakL <- peakL
  }else{
    tpeakU <- peakU*nData/(2*pi) + 1
    tpeakL <- peakL*nData/(2*pi) + 1
  }

  # signal estimation
  phU <- lapply(1:nComp, function(k) (peakU[k]-alpha)/2)
  phL <- lapply(1:nComp, function(k) (peakL[k]-alpha)/2)
  ZU <- sapply(1:nComp, function(k) M+sum(A*cos(beta + 2*atan(omega*tan(phU[[k]])))))
  ZL <- sapply(1:nComp, function(k) M+sum(A*cos(beta + 2*atan(omega*tan(phL[[k]])))))
  names(ZU) <- NULL; names(ZL) <- NULL

  return(list(tpeakU = tpeakU, tpeakL = tpeakL, ZU = ZU, ZL = ZL))
}
