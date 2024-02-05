#'
#' Fitting FMM models
#'
#'
#' \code{fitFMM()} is used to fit FMM models. The only required argument to fit FMM models is the input data.
#' By default it is assumed that time points, corresponding to a single time period, are equally spaced from 0 to \eqn{2\pi}.
#'
#'
#' @param vData A numeric vector containing the data to be fitted a FMM model.
#'
#' @param nPeriods A numeric value specifying the number of periods at which \code{vData} is observed.
#'
#' @param timePoints A numeric vector containing the time points at which each data of one single period is observed. The default value is \code{NULL},
#'   in which case they are equally spaced in range \eqn{[0, 2\pi]}. It must be between 0 and \eqn{2\pi}.
#'
#' @param nback Number of FMM components to be fitted. Its default value is 1.
#'
#' @param betaOmegaRestrictions An integer vector of length \code{nback} indicating which FMM waves are constrained to have equal \code{beta} and \code{omega} parameters. For example, \code{c(1,1,1,2,2)} indicates that \code{beta1=beta2=beta3} and
#'   \code{beta4=beta5} as well as \code{omega1=omega2=omega3} and \code{omega4=omega5}. In brief, some waves are restricted to have the same shape. Its default value is the sequence \code{1:nback} to fit the FMM model without restrictions on shape parameters (\code{beta} and \code{omega}).
#'
#' @param maxiter Maximum number of iterations for the backfitting algorithm. By default, it is set at \code{nback}.
#'
#' @param stopFunction Function to check the convergence criterion for the backfitting algorithm (see Details).
#'
#' @param lengthAlphaGrid Precision of the grid of alpha in the search of the best model. If it is increased, more possible values of alpha will be considered, resulting in an increasing in the computation time too.
#'   By default, it is set to 48 possible values of alpha, equally spaced between 0 and \eqn{2\pi}.
#'
#' @param lengthOmegaGrid Precision of the grid of omega in the search of the best model. If it is increased, more possible values of omega will be considered, resulting in an increasing in the computation time too.
#'   By default it is set to 24 possible values of omega, equally spaced between 0 and 1 in a logarithmic way.
#'
#' @param numReps Number of times (alpha, omega) parameters are refined.
#'
#' @param showProgress \code{TRUE} to display a progress indicator on the console.
#'
#' @param showTime \code{TRUE} to display execution time on the console.
#'
#' @param parallelize \code{TRUE} to use parallelized procedure to fit FMM model. Its default value is \code{FALSE}. When it is \code{TRUE}, the number of cores to be used is equal to 12, or if the machine has less, the number of cores - 1.
#'
#' @param restrExactSolution \code{FALSE} to use an aproximated algorithm to fit the FMM model with restrictions (default). If \code{TRUE} is specified, the exact solution is computed.
#'
#' @details
#' Data will be collected over \code{nPeriods} periods. When \code{nPeriods > 1} the fitting is carried out by averaging the data collected
#' at each time point across all considered periods. The model is fitting to summarized data.
#' \code{timePoints} is a \code{n}-length numeric vector where \code{n} is the number of different time points per period.
#'
#' The \code{stopFunction} argument can either be the functions \code{alwaysFalse} or \code{R2} included in the package or user-defined functions that have the same arguments. The included functions serve for the following:
#' \itemize{
#'   \item{\code{alwaysFalse()}, its default value, which returns \code{FALSE} to force \code{maxiter} iterations; and}
#'   \item{\code{R2(vData,pred,prevPred,difMax = 0.001)}, a function that computes the difference between the explained variability in two consecutive iterations returning \code{TRUE} when the convergence criterion is reached.
#'   To calculate the explained variability difference, the data and the fitted values from the current and previous iteration are passed as arguments \code{vData}, \code{pred} and \code{prevPred}, respectively. The convergence
#'   criterion is fulfilled when the explained variability difference is less than the argument \code{difMax} (by default 0.001).}
#' }
#'
#'
#' @return
#'  An S4 object of class \code{'FMM'} with information about the fitted model. The object contains the following slots:
#'  \describe{
#'    \item{@timePoints}{The time points as specified by the input argument. It is a numeric vector containing the time points at which each data of one single period is observed.}
#'    \item{@data}{The data as specified by the input argument. It is a numeric vector containing the data to be fitted a FMM model. Data could be collected over multiple periods.}
#'    \item{@summarizedData}{When the data has more than one period, a numeric vector containing \code{data} averaging the data at each time point across all considered periods.}
#'    \item{@nPeriods}{A numeric value containing the number of periods in data as specified by the input argument.}
#'    \item{@fittedValues}{A numeric vector of the fitted values by the FMM model.}
#'    \item{@M}{A numeric value of the estimated intercept parameter \eqn{M}.}
#'    \item{@A}{A numeric value or vector of the estimated FMM wave amplitude parameter(s) \eqn{A}.}
#'    \item{@alpha}{A numeric value or vector of the estimated FMM wave phase translation parameter(s) \eqn{\alpha}.}
#'    \item{@beta}{A numeric value or vector of the estimated FMM wave skewness parameter(s) \eqn{\beta}.}
#'    \item{@omega}{A numeric value or vector of the estimated FMM wave kurtosis parameter(s) \eqn{\omega}.}
#'    \item{@SSE}{A numeric value of the sum of the residual squares values.}
#'    \item{@R2}{A numeric vector specifying the explained variance by each of the fitted FMM components.}
#'    \item{@nIter}{A numeric value specifying the number of iterations of the fitting algorithm.}
#'  }
#'
#'
#' @references
#' Rueda C, Larriba Y, Peddada SD (2019). Frequency Modulated Moebius Model Accurately Predicts Rhythmic Signals in Biological and Physical Sciences.
#' \emph{Scientific reports}, \bold{9} (1), 18701. \url{https://www.nature.com/articles/s41598-019-54569-1}
#'
#'
#' @examples
#' # A monocomponent FMM model is fitted.
#' FMM_data <- generateFMM(2, 3, 1.5, 2.3, 0.1,
#'                         from = 0, to = 2*pi, length.out = 100,
#'                         outvalues = TRUE, sigmaNoise = 0.3, plot = FALSE)
#' fit <- fitFMM(FMM_data$y, lengthAlphaGrid = 10, lengthOmegaGrid = 10)
#' summary(fit)
#'
#' # To see the differences between number of repetitions.
#' FMM_data <- generateFMM(2, 1, 1.5, 1.1 ,0.14 , outvalues = TRUE, sigmaNoise = 0.15, plot=TRUE)
#' fit1 <- fitFMM(FMM_data$y, lengthAlphaGrid = 6, lengthOmegaGrid = 3, numReps = 1)
#' fit2 <- fitFMM(FMM_data$y, lengthAlphaGrid = 6, lengthOmegaGrid = 3, numReps = 6,
#'                showProgress = FALSE)  # suppress progress messages
#' getSSE(fit1)
#' getSSE(fit2)
#'
#' # Finer resolution grid.
#' fit3 <- fitFMM(FMM_data$y, lengthAlphaGrid = 10, lengthOmegaGrid = 5, numReps = 1)
#' getSSE(fit3)
#'
#' # Two component FMM model with beta and omega restricted
#' restFMM2w_data <- generateFMM(M = 3, A = c(7, 4), alpha = c(0.5, 5), beta = c(rep(3, 2)),
#'                               omega = rep(0.05, 2), from = 0, to = 2*pi, length.out = 100,
#'                               sigmaNoise = 0.3, plot = FALSE)
#' fit2w.rest <- fitFMM(restFMM2w_data$y, nback = 2, maxiter = 1, numReps = 1,
#'                      lengthAlphaGrid = 15, lengthOmegaGrid = 10,
#'                      betaOmegaRestrictions = c(1, 1))
#' plotFMM(fit2w.rest, components = TRUE)
#'
fitFMM <- function(vData, nPeriods = 1, timePoints = NULL,
                   nback = 1, betaOmegaRestrictions = 1:nback,
                   maxiter = nback, stopFunction = alwaysFalse,
                   lengthAlphaGrid = 48, lengthOmegaGrid = 24,
                   numReps = 3, showProgress = TRUE, showTime = FALSE,
                   parallelize = FALSE, restrExactSolution = FALSE){

  alphaGrid <- seq(0, 2*pi, length.out = lengthAlphaGrid)
  omegaMin <- 0.0001
  omegaMax <- 1
  omegaGrid <- exp(seq(log(omegaMin),log(omegaMax),length.out = lengthOmegaGrid))

  betaOmegaRestrictions <- sort(betaOmegaRestrictions)

  if(showTime) time.ini <- Sys.time()

  # minimum number of observations for fitting
  if(length(vData) < 5){
    stop("The minimum number of observations should be 5")
  }

  # If data has more than one period, it must be summarized
  if(nPeriods > 1){
    n <- length(vData)
    if(n %% nPeriods != 0) stop("Data length is not a multiple of nPeriods")
    dataMatrix <- matrix(vData, nrow = nPeriods, ncol = n/nPeriods, byrow = TRUE)
    summarizedData <- apply(dataMatrix, 2, mean)
  } else {
    summarizedData <- vData
  }

  # If data is constant
  if(sd(vData,na.rm=TRUE) == 0){
    stop("Data are constant")
  }

  # Generation of the time points
  if(is.null(timePoints)){
    timePoints<-seqTimes(length(summarizedData))
  } else {
    if(any(timePoints < 0) | any(timePoints > 2*pi)){
      stop("timePoints must be between 0 and 2*pi")
    }
    if(length(timePoints) != length(summarizedData)){
      stop("timePoints must have the same length as one-period data")
    }
  }

  # Used apply function for compute FMM models
  usedApply_Cluster <- getApply(parallelize = parallelize)
  usedApply <- usedApply_Cluster[[1]]

  ### fitFMM_unit
  if(nback == 1){
    fittedFMM <- fitFMM_unit(vData = summarizedData, timePoints = timePoints,
                             lengthAlphaGrid = lengthAlphaGrid, lengthOmegaGrid = lengthOmegaGrid,
                             alphaGrid = alphaGrid, omegaMin = omegaMin, omegaMax = omegaMax,
                             omegaGrid = omegaGrid, numReps = numReps, usedApply = usedApply)
  } else {
    ### fitFMM_back:
    if(length(unique(betaOmegaRestrictions)) == nback){
      fittedFMM <- fitFMM_back(vData = summarizedData, nback = nback, timePoints = timePoints,
                               maxiter = maxiter, stopFunction = stopFunction,
                               lengthAlphaGrid = lengthAlphaGrid, lengthOmegaGrid = lengthOmegaGrid,
                               alphaGrid = alphaGrid, omegaMin = omegaMin, omegaMax = omegaMax,
                               omegaGrid = omegaGrid, numReps = numReps, showProgress = showProgress,
                               usedApply = usedApply)
    ### fitFMM_restr:
    } else {
      #### Exact solution
      if(restrExactSolution){
        fittedFMM <- fitFMM_restr(vData = summarizedData, nback = nback,
                                  betaRestrictions = betaOmegaRestrictions,
                                  omegaRestrictions = betaOmegaRestrictions,
                                  timePoints = timePoints, maxiter = maxiter, stopFunction = stopFunction,
                                  lengthAlphaGrid = lengthAlphaGrid, lengthOmegaGrid = lengthOmegaGrid,
                                  alphaGrid = alphaGrid, omegaMin = omegaMin, omegaMax = omegaMax,
                                  omegaGrid = omegaGrid, numReps = numReps, parallelize = parallelize)
      #### Approximated solution
      } else {
        fittedFMM <- fitFMM_restr_omega_beta(vData = summarizedData, nback = nback,
                                             betaRestrictions = betaOmegaRestrictions,
                                             omegaRestrictions = betaOmegaRestrictions,
                                             timePoints = timePoints, maxiter = maxiter, stopFunction = stopFunction,
                                             lengthAlphaGrid = lengthAlphaGrid, lengthOmegaGrid = lengthOmegaGrid,
                                             alphaGrid = alphaGrid, omegaMin = omegaMin, omegaMax = omegaMax,
                                             omegaGrid = omegaGrid, numReps = numReps, showProgress = showProgress,
                                             parallelize = parallelize)
      }
    }
  }

  cluster <- usedApply_Cluster[[2]]

  if(!is.null(cluster)) parallel::stopCluster(cluster)

  if(showTime){
    time.end <- Sys.time()
    print(time.end-time.ini)
  }

  fittedFMM@nPeriods <- nPeriods
  fittedFMM@data <- vData


  # Restricted algorithm may find models with A<0
  if(any(getA(fittedFMM) < 0)) {
    stop("Invalid solution: check function input parameters.")
  }

  # "Hack" to add show method without hindering paralellized procedure
  addShowMethod()

  return(fittedFMM)
}




