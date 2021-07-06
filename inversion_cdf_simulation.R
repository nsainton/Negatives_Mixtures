cdf <- function(x, w, mean, sd) {
  res <- 0
  for (i in seq_len(length(w))) {
    res <- res + w[i] * pnorm(x, mean[i], sd[i])
  }
  return(res)
}

invcdf <- function(w, mean, sd, n=1, n_points=100, tol=0.0001) {
  x <- seq(min(mean - 3 * sd), max(mean + 3 * sd), length.out = n_points)
  y <- cdf(x, w, mean, sd)
  u <- runif(n)
  res <- rep(0, n)
  for (i in 1:n) {
    index_sup <- Position(function(x) x > u[i], y)
    if (is.na(index_sup)) {
      #Search y superior to u[i]
      append(x, qnorm((u[i] + 1) / 2, max(mean), max(sd)))
      append(y, cdf(tail(x, 1), w, mean, sd))
      inf <- x[length(x) - 1]
      sup <- x[length(x)]
      index_sup <- length(x)
    }
    else if (index_sup == 1) {
      #Search y inferior to u[i]
      x <- c(qnorm(u[i] / 2, min(mean), max(sd)), x)
      y <- c(cdf(x[1], w, mean, sd), y)
      inf <- x[1]
      sup <- x[2]
      index_sup <- 2
    }
    else{
      inf <- x[index_sup - 1]
      sup <- x[index_sup]
    }
    x_n <- (inf + sup) / 2
    y_n <- cdf(x_n, w, mean, sd)
    while (abs(y_n - u[i]) > tol) {
      if (y_n < u[i]) {
        inf <- x_n
      }
      else{
        sup <- x_n
      }
      x_n <- (inf  + sup) / 2
      y_n <- cdf(x_n, w, mean, sd)
    }
    res[i] <- x_n
  }
  return(res)
}
