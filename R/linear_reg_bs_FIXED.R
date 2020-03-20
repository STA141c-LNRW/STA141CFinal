

linear_reg_bs <- function(x, y, s = 10, r = 1000) {
  n <- dim(x)[1]
  p <- dim(x)[2] + 1
  print(n)
  x1 <- cbind(Intercept = rep(1, n), x)
  sample_indices <- sample(n)
  samples <- sample(s)
  x_samples <- split(x1[sample_indices,], samples)
  y_samples <- split(y[sample_indices], samples)
  bs_coefs <- list()
  bs_s2 <- list()
  for(i in 1:s) {
    sample_coefs <- NULL
    sample_s2 <- NULL
    for (j in 1:r){
      n_sub <- length(y_samples[[i]])
      freqs <- rmultinom(1, n, rep(1, n_sub))
      subset <- data.frame(x_samples[[i]], y_samples[[i]])
      resamp = subset[rep(seq_len(nrow(subset)), freqs),]
      x_resamp <- as.matrix(resamp[,1:p])
      y_resamp <- as.matrix(resamp[,p+1])
      coefs <- solve(t(x_resamp) %*% x_resamp) %*% t(x_resamp) %*% y_resamp
      fv <- x_resamp %*% coefs
      res <- y_resamp - fv
      s2 <- sum(res^2) / (n_sub - p)
      sample_coefs <- cbind(sample_coefs, coefs)
      sample_s2 <- c(sample_s2, s2)
    }
    bs_coefs[[i]] <- sample_coefs
    bs_s2[[i]] <- sample_s2
  }
  return(list(bootstrap_coefficient_estimates = bs_coefs,
              bootstrap_s2_estimates = bs_s2))
}
