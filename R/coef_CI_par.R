#' Confidence Intervals and Estimates of Each Regression Coefficient - Parallelized
#'
#' This function determines overall confidence intervals and estimates for each
#' regression coefficient using parallelization.
#' (EXPAND ON THIS)
#'
#' @param lrbs A linear_reg_bs or linear_reg_bs_par object containing BLB regression
#' coefficient estimates.
#' @param alpha The significance level. Default value is 0.05.
#' @return The overall confidence interval for each regression coefficient, along with
#' its overall estimate.
#' @export
coef_CI_par <- function(lrbs, alpha = 0.05) {
  coef <- lrbs$bootstrap_coefficient_estimates
  CIs <- future_map(coef, function(c) apply(c, 1, quantile,
                                            probs = c((alpha / 2), (1 - (alpha / 2)))))
  means <- future_map(coef, function(c) apply(c, 1, mean))
  CI <- reduce(CIs, `+`) / length(CIs)
  beta_hat <- reduce(means, `+`) / length(means)
  return(cbind(Lower_Bounds = CI[1,], Estimates = beta_hat, Upper_Bounds = CI[2,]))
}
