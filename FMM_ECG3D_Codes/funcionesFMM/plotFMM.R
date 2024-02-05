#' Plot fitted FMM models
#'
#' \code{plotFMM()} is used to plot fitted FMM models. The function can either
#' plot the fitted model against the data or each of the components of the model
#' separately. Optionally \code{'ggplot2'} can be used as graphic library.
#'
#' @param objFMM Object of class FMM
#' @param components A logical value indicating if the centered wave components of the model should be separately
#' plotted (case where it is \code{TRUE}). If \code{FALSE}, the default, the fitted FMM model
#' along with the observed data is plotted.
#' @param plotAlongPeriods A logical value indicating if more than one period should be plotted in the plots
#' by default. Its default value is \code{FALSE}.
#' @param use_ggplot2 A logical value. If \code{FALSE}, the default, R base graphics are used. If \code{TRUE},
#' \code{'ggplot2'} library is used as graphics engine.
#' @param legendInComponentsPlot A logical value indicating whether the legend should be plotted in the components
#' plot. By defaults it is \code{TRUE}.
#' @param textExtra A character vector for extra text to be added to the titles of the plots.
#'
#' @details { \code{plotFMM()} can generate two types of plots: the basic plot compares the fitted model against the original data while the components plot represents separately the centered waves of the model (if the argument components is TRUE).
#'
#' The function is also capable of plotting multiple periods if the data has more than one, as is the case in many applications such as chronobiology. In this case, the argument plotAlongPeriods should be TRUE. In the case of components plots the value taken by the latter argument is ignored as they are plotted along just one period.
#'
#' While, by default, plots are created using base R graphics, 'ggplot2' can also be used for more aesthetic and customizable plots. Optional arguments legendInComponentsPlot and textExtra serve to control, respectively, whether a legend to the components plot should be added and adding extra text to the plot's title.
#' }
#'
#' @return None if base R graphics are used, a named \code{ggplot2} list if \code{'ggplot2'} is used.
#'
#' @examples
#'
#' # Simulates an scenario in which an FMM model is suitable,
#' res <- generateFMM(2,3,1.5,2.3,0.1,outvalues = TRUE,sigmaNoise = 0.3, plot=FALSE)
#' # then a FMM model is fitted to the data.
#' fit <- fitFMM(res$y, lengthAlphaGrid=20,lengthOmegaGrid=12)
#' plotFMM(fit)
#'
#' # Components plot of FMM Model fitted to neuronal data with various optional aesthetics
#' data("neuronalSpike")
#' fittedFMM2<-fitFMM(neuronalSpike, nback=2,
#'                    lengthAlphaGrid = 24,lengthOmegaGrid = 10, numReps = 1)
#'
#' plotFMM(fittedFMM2, components = TRUE)
#' plotFMM(fittedFMM2, components = TRUE,
#'         legendInComponentsPlot = FALSE,
#'         textExtra = "Neuronal Data")
#'
#' # With ggplot2, customizable plots can be created,
#' library(ggplot2)
#' # standard plots
#' plotFMM(fittedFMM2, use_ggplot2 = TRUE)
#' # modify x-axis with original timePoints
#' timePoints <- getTimePoints(fittedFMM2)
#' nObs <- length(timePoints)
#' sTimePoints <- round(c(1, nObs*0.25, nObs*0.5, nObs*0.75, nObs))
#' plotFMM(fittedFMM2, use_ggplot2 = TRUE) +
#'   scale_x_continuous(breaks = sTimePoints,
#'                      labels = function(x) round(timePoints[x],2))
#' # and components plots
#' plotFMM(fittedFMM2, components = TRUE, use_ggplot2 = TRUE)
#'
#' # Plot of fitted model to more than one period.
#' data("mouseGeneExp")
#' fittedFMM2<-fitFMM(mouseGeneExp, nPeriods = 2,
#'                    lengthAlphaGrid = 20,lengthOmegaGrid = 10)
#' plotFMM(fittedFMM2, plotAlongPeriods = TRUE)

plotFMM <- function(objFMM, components = FALSE, plotAlongPeriods = FALSE,
                    use_ggplot2 = FALSE, legendInComponentsPlot = TRUE, textExtra = ""){

  nPeriods <- getNPeriods(objFMM)
  if(nPeriods > 1){
    if(plotAlongPeriods & !components){
      vData <- getData(objFMM)
    }else{
      vData <- getSummarizedData(objFMM)
    }
  }else{vData <- getData(objFMM)}
  nObs <- length(vData)

  if(plotAlongPeriods & !components){
    timePoints <- getTimePoints(objFMM)
    timePoints <- rep(timePoints, nPeriods)
  }else{
    timePoints <- getTimePoints(objFMM)
  }

  significantTimePoints <- round(c(1, nObs*0.25, nObs*0.5, nObs*0.75, nObs))

  # Components plot: if there is more than one period, just the data from the first period will be plotted
  if(components){
    title <- ifelse(textExtra != "", paste("Components FMM", textExtra, sep = " - "),"Components FMM")
    nComponents <- length(getAlpha(objFMM))
    # With more than 9 components, the selection of colors must be expanded
    if(nComponents > 9){
      colorsForComponents <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(nComponents)
    }else{
      colorsForComponents <- ifelse(rep(nComponents>3,nComponents),
                                    RColorBrewer::brewer.pal(nComponents, "Set1"),
                                    RColorBrewer::brewer.pal(3, "Set1"))
    }
    componentNames<-paste("Wave ", 1:nComponents, sep = "")

    predicted <- extractWaves(objFMM)

    if(!use_ggplot2){
      yLimits<-c(min(sapply(predicted, min)), max(sapply(predicted, max)))
      plot(1:nObs, vData, ylim = yLimits, xlab = "Time", ylab = "Response",
           main = title, type = "n", xaxt = "n")
      for(i in 1:nComponents){
        points(1:nObs, predicted[[i]], type = "l", lwd = 2, col = colorsForComponents[i])
      }
      axis(1, las = 1, at = significantTimePoints,
           labels = parse(text=paste("t[",significantTimePoints, "]", sep = "")))
      if(legendInComponentsPlot) legend("topright", legend = componentNames, col = colorsForComponents, lty = 1)
    } else {
      requireNamespace("ggplot2", quietly = TRUE)
      requireNamespace("RColorBrewer", quietly = TRUE)

      df <- data.frame("Time" = rep(1:length(timePoints), nComponents),
                       "Response" = unlist(predicted),
                       "Components" = rep(componentNames, each = nObs))

      plot<-ggplot2::ggplot(data = df, ggplot2::aes_(x=~Time, y=~Response, group =~ Components,
                                                     color =~ Components)) +
        ggplot2::geom_line(ggplot2::aes_(color =~ Components),
                           size=1.3,lineend = "round",linejoin = "round") +
        ggplot2::scale_color_manual(values = colorsForComponents) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = ifelse(legendInComponentsPlot,"bottom","none")) +
        ggplot2::labs(title = title) +
        ggplot2::scale_x_continuous(breaks = significantTimePoints,
                                    labels = function(x) parse(text=paste("t[",x,"]")))
      return(plot)
    }

  } else {
    title <- ifelse(textExtra != "", paste("Fitted FMM model",textExtra,sep = " - "),"Fitted FMM model")

    if(!use_ggplot2){
      yLimits<-c(min(vData,getFittedValues(objFMM)), max(vData,getFittedValues(objFMM)))
      plot(1:nObs, vData, xlab = "Time", ylab = "Response", main = title, xaxt = "n",
           ylim = yLimits)
      if(plotAlongPeriods){
        points(1:nObs, rep(getFittedValues(objFMM), nPeriods), type = "l", col = 2, lwd = 2)
      }else{
        points(1:nObs, getFittedValues(objFMM), type = "l", col = 2, lwd = 2)
      }
      axis(1, las = 1, at = significantTimePoints,
           labels = parse(text=paste("t[",significantTimePoints, "]", sep = "")))
    } else {
      requireNamespace("ggplot2", quietly = TRUE)

      if(plotAlongPeriods){
        adjustedModel<-rep(getFittedValues(objFMM),nPeriods)
      }else{
        adjustedModel<-getFittedValues(objFMM)
      }

      fittedData <- data.frame("Time" = 1:nObs, "fitted_FMM" = adjustedModel, "Response" = vData)

      plot <- ggplot2::ggplot(data = fittedData, ggplot2::aes_(x=~Time, y=~Response, color = 1)) +
        ggplot2::geom_point(size = 2, color = "grey65", shape = 21, stroke = 1.1) +
        ggplot2::geom_path(ggplot2::aes_(x=~Time, y=~fitted_FMM, color = "FMM", position = NULL),
                           size=2, lineend = "round", linejoin = "round")+
        ggplot2::labs(title = title) +
        ggplot2::scale_color_manual(values = "red") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::scale_x_continuous(breaks = significantTimePoints,
                                    labels = function(x) parse(text = paste("t[",x,"]")))
      return(plot)
    }
  }
}
