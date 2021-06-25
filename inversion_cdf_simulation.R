cdf <- function(x, w, mean, sd) {
  res <- 0
  for (i in seq_len(length(w))) {
    res <- res + w[i] * pnorm(x, mean[i], sd[i])
  }
  return(res)
}
invcdf <- function(q, w, mean, sd, n=1, n_points=100) {
  y <- seq(min(mean - 3 * sd), max(mean + 3 * sd), length.out = n_points)
  x <- cdf(y, w, mean, sd)
  inv <- approxfun(x, y, yleft = 0, yright = 1, method = "linear")
  u <- runif(n)
  return(inv(u))
}
