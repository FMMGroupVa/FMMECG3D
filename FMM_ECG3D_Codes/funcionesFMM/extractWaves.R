#' Individual contribution to the fitted values of each FMM wave
#'
#'
#' \code{extractWaves} extracts individual contribution to the fitted values of each FMM wave.
#'
#'
#' @param objFMM Object of class \code{'FMM'}.
#'
#'
#' @return
#' Individual contribution to the fitted values of each FMM wave. It is a \code{list} object with as many elements as FMM components have been fitted.
#'
#'
#' @examples
#' ## Generate example data:
#' fmm2.data <- generateFMM(M = 0, A = rep(1, 2),
#'                          alpha = c(1.5, 3.4), beta = c(0.2, 2.3), omega = c(0.1, 0.2),
#'                          plot = FALSE, outvalues = TRUE,
#'                          sigmaNoise = 0.5) # add a gaussian noise with sigma = 0.5
#'
#' ## Fit the FMM model with nback = 2 components
#' ## fit is an object of S4 class 'FMM'
#' fit <- fitFMM(fmm2.data$y,timePoints = fmm2.data$t,nback = 2,
#'               lengthAlphaGrid = 24,lengthOmegaGrid = 10)
#' ## extracts individual contribution of each FMM wave
#' extractWaves(fit)
#'
extractWaves <- function(objFMM){
   nComponents <- length(getAlpha(objFMM))
   timePoints <- getTimePoints(objFMM)
   firstValue <- getData(objFMM)[1]
   predicted <- list()

   for(i in 1:nComponents){
     predictedComponent <- getA(objFMM)[i]*calculateCosPhi(getAlpha(objFMM)[i], getBeta(objFMM)[i],
                                                       getOmega(objFMM)[i], timePoints)
     predicted[[i]] <- predictedComponent - predictedComponent[1] + firstValue
   }
   return(predicted)
}
