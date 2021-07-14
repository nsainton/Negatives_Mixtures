source("mixture_functions.R")
source("q_min.R")

invcdf <- function(w, means, sd, n=1, n_points=100, tol=0.0001) {
  x <- rep(NA, n + n_points + 2)
  x[1:n_points + 1] <- seq(min(means - 3 * sd), max(means + 3 * sd), length.out = n_points)
  y <- pmixture(x[1:n_points + 1], w, means, sd)
  u <- runif(n)
  y <- c(0, y, u, 1)
  ordered_index <- order(y)
  y <- y[ordered_index]
  x <- x[ordered_index]
  inverse_order_index <- order(ordered_index)
  known <- !is.na(x)
  res <- rep(0, n)
  left_tail <- q_min(min(u), w, means, sd, x[known][2], y[known][2], x[known][3], y[known][3])
  x[1] <- left_tail[1]
  y[1] <- left_tail[2]
  right_tail <- q_max(max(u), w, means, sd, x[known][n_points], y[known][n_points], x[known][n_points - 1], y[known][n_points - 1])
  x[length(x)] <- right_tail[1]
  y[length(y)] <- right_tail[2]
  known[1] <- TRUE
  known[length(known)] <- TRUE
  for (i in 1:n) {
    index_sup <- sum(y < u[i]) + 1
    index_inf <- index_sup - 1
    while (!known[index_inf]) {
      index_inf <- index_inf - 1
    }
    while (!known[index_sup]) {
      index_sup <- index_sup + 1
    }
    inf <- x[index_inf]
    sup <- x[index_sup]
    x_n <- (inf + sup) / 2
    y_n <- pmixture(x_n, w, means, sd)
    while (abs(y_n - u[i]) > tol) {
      if (y_n < u[i]) {
        inf <- x_n
      }
      else{
        sup <- x_n
      }
      x_n <- (inf  + sup) / 2
      y_n <- pmixture(x_n, w, means, sd)
    }
    res[i] <- x_n
    x[inverse_order_index[i + n_points + 1]] <- x_n
    y[inverse_order_index[i + n_points + 1]] <- y_n
    known[inverse_order_index[i + n_points + 1]] <- TRUE
  }
  return(res)
}
