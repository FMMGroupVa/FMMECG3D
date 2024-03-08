# S4 object of class FMM
FMM <- setClass("FMM", slots = c(
  timePoints="numeric",
  data="numeric",
  summarizedData="numeric",
  nPeriods="numeric",
  fittedValues="numeric",
  M="numeric",
  A="numeric",
  alpha="numeric",
  beta="numeric",
  omega="numeric",
  SSE="numeric",
  R2="numeric",
  nIter="numeric")
)



# Extracts the slots from the S4 object
getM <- function(objFMM) { objFMM@M }
getA <- function(objFMM) { objFMM@A }
getAlpha <- function(objFMM) { objFMM@alpha }
getBeta <- function(objFMM) { objFMM@beta }
getOmega <- function(objFMM) { objFMM@omega }
getNPeriods <- function(objFMM) { objFMM@nPeriods }
getTimePoints <- function(objFMM) { objFMM@timePoints }
getData <- function(objFMM) { objFMM@data }
getSummarizedData <- function(objFMM) { objFMM@summarizedData }
getFittedValues <- function(objFMM) { objFMM@fittedValues }
getSSE <- function(objFMM) { objFMM@SSE }
getR2 <- function(objFMM) { objFMM@R2 }
getNIter <- function(objFMM) { objFMM@nIter }


