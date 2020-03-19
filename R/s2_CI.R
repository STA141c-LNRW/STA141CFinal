#' Confidence Intervals and Estimates of Sigma-Squared
#'
#' This function determines an overall confidence interval and estimate for
#' sigma-squared.
#' (EXPAND ON THIS)
#'
#' @param lrbs A linear_reg_bs or linear_reg_bs_par object containing BLB sigma-squared
#' estimates.
#' @param alpha The significance level. Default value is 0.05.
#' @return The overall confidence interval for sigma-squared, along with its overall
#' estimate.
#' @export
s2_CI <- function(lrbs, alpha = 0.05) {
  s2 <- lrbs$bootstrap_s2_estimates
  CIs <- map(s2, quantile, probs = c((alpha / 2), (1 - (alpha / 2))))
  means <- map(s2, mean)
  CI <- reduce(CIs, `+`) / length(CIs)
  s2_hat <- reduce(means, `+`) / length(means)
  return(c(Lower_Bound = CI[[1]], Estimate = s2_hat, Upper_Bound = CI[[2]]))
}
