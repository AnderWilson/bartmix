
#' Total mixture effect
#'
#' @param object A fitted object of class wbart from the BART R package.
#' @param data The data used to fit the BART model.
#' @param exposures Vector with the variable names of exposures in the data.
#' @param qtls A vector of quantiles to estimate the total mixture effect at.
#'
#' @return A data.frame with the total mixture effect.
#'
#' @importFrom stats quantile
#' @import BART
#'
#' @export
#'
totalmixtureeffect <- function(object,
                               data,
                               exposures=colnames(data),
                               qtls = seq(0.2,0.8,0.05)){

  # sample size
  n <- nrow(data)

  # quantiles for mixture effect
  x <- apply(data[,exposures],2,quantile,qtls)
  L <- nrow(x)

  # data to predict with
  x_partialplot_temp <- data
  x_partialplot <- NULL
  for(j in 1:L){
    x_partialplot_temp[,exposures] <- rep(x[j,], each=n)
    x_partialplot <- rbind(x_partialplot, x_partialplot_temp)
  }

  # predict outcomes
  pred <- BART:::predict.wbart(object, x_partialplot)

  # summarize
  # place to store results
  partial <- matrix(nrow = nrow(pred), ncol = L)

  # center on value closest to 0.5
  sset_center <- (which.min(abs(qtls-0.5)) - 1) * n + 1:n
  partial_center <- rowMeans(pred[, sset_center])

  # calculate for all quantiles
  for(j in 1:L) {
    sset <- (j - 1) * n + 1:n
    partial[, j] <- rowMeans(pred[, sset])-partial_center
  }


  # summarize results
  out <- data.frame(
    quantile = qtls,
    mean = colMeans(partial),
    lower = apply(partial,2,quantile,0.025),
    upper = apply(partial,2,quantile,0.975)
  )

  return(out)

}
