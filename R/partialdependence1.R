
#' 1-way partial dependence for a single exposure
#'
#' @param object A fitted object of class wbart from the BART R package.
#' @param data The data used to fit the BART model.
#' @param varname Name of exposure  to calculate the partial dependence function for.
#' @param L The number of points to estimate the exposure-response function of var at.
#'
#' @return A data.frame with the estimated exposure-response function.
#'
#' @importFrom stats quantile
#' @import BART

partialdependence1_helper <- function(object,
                               data,
                               varname,
                               L=50){

  # sample size
  n <- nrow(data)

  # grid to predict over
  x <- seq(min(data[, varname]), max(data[, varname]), length.out = L)

  # data to predict with
  x_partialplot_temp <- data
  x_partialplot <- NULL
  for(j in 1:L){
    x_partialplot_temp[,varname] <- x[j]
    x_partialplot <- rbind(x_partialplot, x_partialplot_temp)
  }

  # predict outcomes
  pred <- BART:::predict.wbart(object, x_partialplot)

  # summarize
  partial <- matrix(nrow = nrow(pred), ncol = L)
  for(j in 1:L) {
    sset <- (j - 1) * n + 1:n
    partial[, j] <- rowMeans(pred[, sset])
  }

  # summarize results
  out <- data.frame(
    x = x,
    mean = colMeans(partial),
    lower = apply(partial,2,quantile,0.025),
    upper = apply(partial,2,quantile,0.975)
  )

  return(out)

}



#' 1-way partial dependence
#'
#' @param object A fitted object of class wbart from the BART R package.
#' @param data The data used to fit the BART model.
#' @param exposures Names of exposures to calculate the partial dependence function for. This can be the name of a single exposure or a vector of names.
#' @param L The number of points to estimate the exposure-response function of var at.
#'
#' @return A data.frame with the estimated partial dependence function(s).
#'
#' @export
#'
partialdependence1 <- function(object,
                               data,
                               exposures,
                               L=50){

  all_pd <- data.frame()

  for(exposure  in exposures){

    cat("calculating exposure",which(exposures==exposure),"of",length(exposures))

    temp_pd <- partialdependence1_helper(object,
                                  data=data,
                                  varname = exposure,
                                  L=L)

    temp_pd$exposure <- exposure
    all_pd <- rbind(all_pd,temp_pd)
  }

  return(all_pd)
}
