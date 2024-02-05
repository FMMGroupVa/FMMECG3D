################################################################################
#  FUNCTION:                  DESCRIPTION:
#  summary.FMM                summary method for S4 class 'FMM'
#  coef.FMM                   coef method for S4 class 'FMM'
#  fitted.FMM                 fitted method for S4 class 'FMM'
#  resid.FMM                  resid method for S4 class 'FMM'
#  show.FMM                   show method for S4 class 'FMM'
################################################################################

# ------------------------------------------------------------------------------
setMethod("summary", signature(object="FMM"),function(object,...) {

  # Checks
  if(!is(object, "FMM")){stop("Object must be of class 'FMM'")}

  x <- object
  nComp <- max(c(length(getAlpha(x)),length(getBeta(x)),length(getOmega(x))),na.rm = TRUE)

  # Title:
  cat("\nTitle:\n", "FMM model with ", nComp, " components", "\n", sep = "")

  # Coefficients:
  if (!is.null(getM(x)) & !is.null(getA(x)) & !is.null(getAlpha(x)) & !is.null(getBeta(x)) & !is.null(getOmega(x))) {
    coef.list <- coef(x)
    cat("\nCoefficients:\n")
    cat(paste("M (Intercept): ", round(coef.list$M, digits = 4),"\n", sep=""))
    print(round(coef.list$wave, digits = 4))
    cat("\n")
  }


  # Peak and trough time estimates:
  fp.est <- NULL
  if (!is.null(getM(x)) & !is.null(getA(x)) & !is.null(getAlpha(x)) & !is.null(getBeta(x)) & !is.null(getOmega(x))) {
    cat("Peak and trough times and signals:\n")

    fidPoints <- getFMMPeaks(x, timePointsIn2pi = TRUE)

    fp.est<-cbind(fidPoints$tpeakU,fidPoints$ZU,fidPoints$tpeakL,fidPoints$ZL)
    rownames(fp.est) <- paste("FMM wave ",1:nComp,": ", sep="")
    colnames(fp.est) <- c("t.Upper","Z.Upper","t.Lower","Z.Lower")
    print(round(fp.est, digits = 4))
    cat("\n")
  }

  # Residuals:
  resid <- resid(x)
  cat("Residuals:\n")
  print(summary(resid))


  # R-squared:
  R.sq <- NULL
  if (!is.null(getR2(x))){
    R.sq <- getR2(x)
    if(nComp > 1){
      R.sq <- c(R.sq,sum(getR2(x)))
      names(R.sq) <- c(paste("Wave",1:nComp,sep=" "),"Total")
    }
    cat("\nR-squared:\n")
    print(round(R.sq, digits = 4))
    cat("\n")
  }

  df.sum <- list(coef = coef.list, peak.time = fp.est, resid = resid, R.squared = R.sq)
  invisible(df.sum)
})


# ------------------------------------------------------------------------------
# coef method for S4 class FMM
#
# return:   data.frame with parameters estimates
#
setMethod("coef", "FMM",function(object,...) {
  # Checks
  if(!is(object, "FMM")){stop("Object must be of class 'FMM'")}

  x <- object
  nComp <- max(c(length(getAlpha(x)),length(getBeta(x)),length(getOmega(x))),na.rm = TRUE)

  coef.wave <- NULL
  # Coefficients:
  if (!is.null(getM(x)) & !is.null(getA(x)) & !is.null(getAlpha(x)) & !is.null(getBeta(x)) & !is.null(getOmega(x))) {
    Names <- c("A","alpha","beta","omega")
    coef.wave <- as.data.frame(t(sapply(1:nComp, function(i)
      c(getA(x)[i],getAlpha(x)[i],getBeta(x)[i],getOmega(x)[i]))))
    rownames(coef.wave) <- paste("FMM wave ",1:nComp,": ", sep="")
    colnames(coef.wave) <- Names
    coef.list <- list(M = getM(x), wave = coef.wave)
  }
  return(coef.list)
})


# ------------------------------------------------------------------------------
# fitted method for S4 class FMM
#
# return:   data.frame with two columns: time and fitted,
#           the timePoints and fitted values, respectively
#
setMethod("fitted", "FMM", function(object,...) {
  # Checks
  if(!is(object, "FMM")){stop("Object must be of class 'FMM'")}

  x <- object

  # fitted values:
  fitted.df <- data.frame(timePoints = getTimePoints(x), fitted = getFittedValues(x))

  return(fitted.df)
})


# ------------------------------------------------------------------------------
# resid method for S4 class FMM
#
# return:   residual vector
#
setMethod("resid", "FMM", function(object,...) {
  # Checks
  if(!is(object, "FMM")){stop("Object must be of class 'FMM'")}

  x <- object

  data <-getData(x)
  if(getNPeriods(x) > 1) data <- getSummarizedData(x)

  # residuals:
  res <- NULL
  if (!is.null(data) & !is.null(data - getFittedValues(x))) {
    res <- data - getFittedValues(x)
  }

  return(res)
})

addShowMethod<-function(){
  rlang::env_unlock(env = asNamespace('FMM'))
  rlang::env_binding_unlock(env = asNamespace('FMM'))
  setMethod("show", signature(object="FMM"),function(object) {
    #
    # Checks
    if(!is(object, "FMM")){stop("Object must be of class 'FMM'")}

    text<-""

    x <- object
    nComp <- max(c(length(getAlpha(x)),length(getBeta(x)),length(getOmega(x))),na.rm = TRUE)

    # Title:
    cat("\nTitle:\n", "FMM model with ", nComp, " components", "\n", sep = "")

    # Coefficients:
    if (!is.null(getM(x)) & !is.null(getA(x)) & !is.null(getAlpha(x)) & !is.null(getBeta(x)) & !is.null(getOmega(x))) {
      coef.list <- coef(x)
      cat("\nCoefficients:\n")
      cat(paste("M (Intercept): ", round(coef.list$M, digits = 4),"\n", sep=""))
      print(round(coef.list$wave, digits = 4))
      cat("\n")
    }

    # R-squared:
    R.sq <- NULL
    if (!is.null(getR2(x))){
      R.sq <- getR2(x)
      if(nComp > 1){
        R.sq <- c(R.sq,sum(getR2(x)))
        names(R.sq) <- c(paste("Wave",1:nComp,sep=" "),"Total")
      }
      cat("\nR-squared:\n")
      print(round(R.sq, digits = 4))
      cat("\n")
    }

  })

  rlang::env_binding_lock(env = asNamespace('FMM'))
  rlang::env_lock(asNamespace('FMM'))
}
