rm(list = ls())
# We begin by creating the first function we need
# to generate our boxes that will be called boxes.
# This function we'll take the bounds we need as an input,
# the number of boxes, the mean and the
# standard deviation of our normal law and
# the type of boxes we need (over or under the curve)
# And will generate two times the number of boxes to approximate the curve with

boxes <- function(inf, sup, number_of_boxes, mean, sd, coef = 1) {
  bound <- ifelse(abs(mean - inf) > abs(mean - sup), mean + abs(mean - inf), sup)
  # We're just getting sure that the right bound is the more distant
  # to the mean of our normal law between the left and the right bound
  step <- (bound - mean) / number_of_boxes
  if (coef > 0) {
    # If the type is "over" we generate one step further to go one step further on the other side by symmetry
    # We define the intervals accordingly and generate the heights on the giver intervals
    right_interval <- seq(mean, bound + step, step)
    left_interval <- sort(-right_interval) + 2 * mean
    length_right <- length(right_interval) - 1
    left_interval <- left_interval[1:length_right]
    right_interval <- right_interval[1:length_right]
    right_heights <- dnorm(right_interval, mean, sd)
    left_heights <- sort(right_heights)[1:length_right]
  } else if (coef < 0) {
    # If the type is "under" we also generate one step further in order to translate it of one step
    # we generate and cut our intervals at the end of the generation
    right_interval <- seq(mean, bound + step, step)
    right_heights <- dnorm(right_interval, mean, sd)[1:length(right_interval)]
    left_heights <- sort(right_heights)[1:length(right_heights) - 1]
    right_heights <- right_heights[2:length(right_heights)]
    right_interval <- right_interval[1:length(right_interval) - 1]
    left_interval <- sort(-right_interval) + 2 * mean - step
  }
  intervals <- c(left_interval, right_interval)
  heights <- c(left_heights, right_heights)
  heights <- coef * heights
  return(list("intervals" = intervals, "heights" = heights))
}

moyenne <- 12
std <- 5
weights <- c(2, -1)
means <- c(2, 3)
standard_deviations <- c(2, 5)
co <- 4
a <- boxes(inf = -1, sup = 1, mean = moyenne, number_of_boxes = 20, sd = std, coef = co)
x <- a$intervals
y <- a$heights
plot(x, y, type = "s")
lines(x, co * dnorm(x, mean = moyenne, sd = std))
# b=generate_boxes(weights,means,standard_deviations)
