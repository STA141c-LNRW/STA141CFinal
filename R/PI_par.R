#' Prediction Intervals and Estimates for New Data - Parallelized
#'
#' This function determines prediction intervals and estimates for new data using the
#' BLB regression coefficient estimates and parallelization.
#' (EXPAND ON THIS)
#'
#' @param lrbs A linear_reg_bs or linear_reg_bs_par object containing BLB regression
#' coefficient estimates.
#' @param x A dataframe of the explanatory variables of unseen observations.
#' @param alpha The significance level. Default value is 0.05.
#' @return The prediction intervals and estimates for the response variable of each
#' unseen observation.
#' @export
PI_par <- function(lrbs, x, alpha = 0.05) {
  coefs <- lrbs$bootstrap_coefficient_estimates
  x1 <- as.matrix(cbind(Intercept = 1, x))
  preds <- future_map(1:length(coefs), function(i) {
    future_map(1:dim(coefs[[i]])[2], function(j) x1 %*% coefs[[i]][, j])
  })
  preds <- future_map(preds, function(sample) matrix(unlist(sample), nrow = dim(x1)[1]))
  PIs <- future_map(preds, function(p) apply(p, 1, quantile,
                                             probs = c((alpha / 2), (1 - (alpha / 2)))))
  fits <- future_map(preds, function(p) apply(p, 1, mean))
  PI <- reduce(PIs, `+`) / length(PIs)
  fit <- reduce(fits, `+`) / length(fits)
  return(cbind(Lower_Bounds = PI[1,], Estimates = fit, Upper_Bounds = PI[2,]))
}
