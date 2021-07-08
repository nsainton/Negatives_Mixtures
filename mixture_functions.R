dmixture <- function(x, weights, means, sd) {
  normal_laws <- lapply(seq_len(length(weights)), function(i) {
    weights[i] * dnorm(x, mean = means[i], sd = sd[i])
  })
  Reduce("+", normal_laws)
}

pmixture <- function(x, weights, means, sd) {
  normal_laws <- lapply(seq_len(length(weights)), function(i) {
    weights[i] * pnorm(x, mean = means[i], sd = sd[i])
  })
  Reduce("+", normal_laws)
}