

#' 2-way partial dependence
#'
#' @param object A fitted object of class wbart from the BART R package.
#' @param data The data used to fit the BART model.
#' @param var The name of the exposure of interest in data. The exposure-response function for this variable will be produced
#' @param var2 The name of a second exposure variable. The exposure-response function for var will be produced at specified quantiles of var2.
#' @param qtls The quantiles of var2 to be used.
#' @param L The number of points to estimate the exposure-response function of var at.
#'
#' @return A data.frame with the exposure-response functions
#' @importFrom stats quantile
#' @import BART
#'
#' @export


partialdependence2 <- function(object,
                               data,
                               var,
                               var2,
                               qtls=c(0.25,0.5,0.75),
                               L=50){

  # sample size
  n <- nrow(data)

  # grid to predict over
  x <- seq(min(data[, var]), max(data[, var]), length.out = L)
  x2 <- quantile(data[,var2],qtls)

  # data to predict with
  x_partialplot_temp <- data
  x_partialplot <- NULL
  for(j in 1:L){
    for(j2 in 1:length(qtls)){
      x_partialplot_temp[,var] <- x[j]
      x_partialplot_temp[,var2] <- x2[j2]
      x_partialplot <- rbind(x_partialplot, x_partialplot_temp)
    }
  }

  # predict outcomes
  pred <- BART:::predict.wbart(object, x_partialplot)

  # summarize
  partial <- matrix(nrow = nrow(pred), ncol = L*length(qtls))
  for(j in 1:L) {
    for(j2 in 1:length(qtls)){
      sset <- (j - 1) * n * length(qtls) + (j2 - 1) * n + 1:n
      partial[, (j-1)*length(qtls)+j2] <- rowMeans(pred[, sset])
    }
  }

  # summarize results
  out <- data.frame(
    x = rep(x,each=length(qtls)),
    qtl = rep(qtls, L),
    mean = colMeans(partial),
    lower = apply(partial,2,quantile,0.025),
    upper = apply(partial,2,quantile,0.975)
  )

  return(out)

}
