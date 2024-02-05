#'
#' Simulating data from FMM models
#'
#'
#' \code{generateFMM()} simulates data from a FMM model defined by parameters \code{M}, \code{A}, \eqn{\alpha}, \eqn{\beta} and \eqn{\omega}.
#'
#'
#' @param M A numeric vector which contains the value of the intercept parameter \code{M}.
#' @param A A positive numeric vector which contains the value of the FMM wave amplitude parameter \code{A}.
#' @param alpha A numeric vector which contains the value of the FMM wave phase translation parameter \eqn{\alpha}.
#' @param beta A numeric vector which contains the value of the FMM wave skewness parameter \eqn{\beta}.
#' @param omega A numeric vector which contains the value of the FMM wave kurtosis parameter \eqn{\omega}. \code{omega} parameter must be between 0 and 1.
#' @param from A numeric value which contains the initial time point of the simulated data. By default, it is 0.
#' @param to A numeric value which contains the final time point of the simulated data. By default, it is \eqn{2\pi}.
#' @param length.out A non-negative number wich contains the desired length of the simulation. By default, it is 100.
#' @param plot A logical value indicating whether the simulated data should be drawn on a plot. By default, it is \code{TRUE}.
#' @param outvalues A logical value indicating whether the numerical simulation should be return. By default, it is \code{TRUE}.
#' @param sigmaNoise A non-negative number which contains the standard deviation of the gaussian noise to be added. Its default value is zero equivalent to a simulation set-up without noise.
#'
#'
#' @details
#' To simulate a multicomponent FMM model, arguments \code{A}, \code{alpha}, \code{beta} and \code{omega} are vectors of length \eqn{m}, where \eqn{m} represents the number of FMM waves. With different lengths, the smaller vectors will be replicate until thay are the same length as the longest vector.
#'
#' With \code{sigmaNoise = s}, \code{s>0}, the \code{generateFMM} function uses \code{rnorm(length.out, 0, sigmaNoise)} to create the normally distributed noise and adds it to the simulated values.
#'
#'
#' @return
#' When \code{outvalues = TRUE} a list of with the following components is returned:
#' \item{input}{a list with the input parameters \code{M}, \code{A}, \code{alpha}, \code{beta} and \code{omega}.}
#' \item{t}{a numeric vector with the time points at each data is simulated.}
#' \item{y}{a numeric vector with the data simulated.}
#'
#' When \code{plot = TRUE} a scatter plot of \code{y} vs \code{t} is drawn.
#'
#'
#' @references
#' Rueda C, Larriba Y, Peddada SD (2019).
#' Frequency Modulated Moebius Model Accurately Predicts Rhythmic Signals in Biological and Physical Sciences.
#' \emph{Scientific reports}, \bold{9} (1), 18701. \url{https://www.nature.com/articles/s41598-019-54569-1}
#'
#'
#' @examples
#' # Simulate data from a monocomponent FMM model. A plot with the simulated model is shown
#' generateFMM(M = 2,A = 3,alpha = 1.5,beta = 2.3, omega = 0.1, outvalues = FALSE)
#'
#' # Add a gaussian noise with standard deviation 0.3. The numeric results are returned
#' generateFMM(M = 2, A = 3, alpha = 1.5, beta = 2.3, omega = 0.1,
#'             sigmaNoise = 0.3, plot = FALSE, outvalues = TRUE)
#'
#' # Simulate data from a multicomponent FMM model with two FMM waves
#' # both with amplitude parameter = 2
#' generateFMM(M = 0, A = rep(2, 2), alpha = c(1.5, 3.4), beta = c(0.2, 2.3), omega = c(0.1, 0.2))
#'
#'
generateFMM <- function(M, A, alpha, beta, omega, from = 0, to = 2*pi, length.out = 200,
                        plot=TRUE, outvalues = TRUE, sigmaNoise = 0){

  timePoints = seqTimes(length.out)
  nArgs <- max(length(M), length(A), length(alpha), length(beta), length(omega))
  minLength <- min(length(A), length(alpha), length(beta), length(omega))

  if(length(M)>1){
    warning("'M' parameter should be a vector of length 1. Using as 'M' the sum of the specified argument.")
    M <- sum(M)
  }

  if(minLength != nArgs){
    warning("'A', 'alpha', 'beta' and 'omega' parameters have different lengths.")
  }

  A <- rep(A, length.out = nArgs)
  if(sum(A <= 0) > 0) stop("'A' parameter must be positive.")

  alpha <- rep(alpha, length.out = nArgs)
  alpha <- alpha%%(2*pi) # between 0 and 2*pi

  beta <- rep(beta, length.out = nArgs)
  beta <- beta%%(2*pi) # between 0 and 2*pi

  omega <- rep(omega, length.out = nArgs)
  if(any(omega<0) | any(omega>1)) stop("'omega' parameter must be between 0 and 1.")

  y <- as.vector(A%*%t(calculateCosPhi(alpha = alpha, beta = beta, omega = omega,
                                       timePoints = timePoints)) + M)

  if (sigmaNoise > 0) y <- y + rnorm(n = length.out, mean = 0, sd = sigmaNoise)

  if(plot) {
    lineType <- ifelse(sigmaNoise == 0, "l", "p")
    plot(timePoints, y, type = lineType, lwd = 2, col = 2, xlab = "Time",
         ylab = "Response", main = paste("Simulated data from FMM model"))
  }

  if(outvalues) return(list(input=list(M = M, A = A, alpha = alpha, beta = beta,
                                       omega = omega), t = timePoints, y = y))
}
