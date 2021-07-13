if (!exists("pmixture", mode = "function")) source("mixture_functions.R")
q_min <- function(q, w, means, sd, left_x, left_y, second_left_x, second_left_y) {
  m <- (left_y - second_left_y) / (left_x - second_left_x)
  b <- left_y - m * left_x
  xres <- -b / m
  yres <- pmixture(xres, w, means, sd)
  while (q < yres) {
    second_left_x <- xres
    second_left_y <- yres
    m <- (left_y - second_left_y) / (left_x - second_left_x)
    b <- left_y - m * left_x
    xres <- -b / m
    yres <- pmixture(xres, w, means, sd)
  }
  return(c(xres, yres))
}
