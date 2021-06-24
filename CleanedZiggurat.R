rm(list = ls())
# We begin by creating the first function we need
# to generate our boxes that will be called boxes.
# This function we'll take the bounds we need as an input,
# the number of boxes, the mean and the
# standard deviation of our normal law and
# the type of boxes we need (over or under the curve)
# And will generate two times the number of boxes to approximate the curve with

boxes <- function(inf, sup, number_of_boxes, mean, sd, coef) {
  bound <- ifelse(abs(mean - inf) > abs(mean - sup), mean + abs(mean - inf), sup)
  # We're just getting sure that the right bound is the more distant
  # to the mean of our normal law between the left and the right bound
  step <- (bound - mean) / number_of_boxes
  if (coef > 0) {
    # If the coefficient is positive we generate one step further
    # to go one step further on the other side by symmetry
    # We define the intervals accordingly and
    # generate the heights on the giver intervals
    right_interval <- seq(mean, bound + step, step)
    left_interval <- rev(-right_interval) + 2 * mean
    length_right <- length(right_interval) - 1
    left_interval <- left_interval[1:length_right]
    right_interval <- right_interval[1:length_right]
    right_heights <- dnorm(right_interval, mean, sd)
    left_heights <- rev(right_heights)[1:length_right]
  } else {
    # If the coefficient is negative we also generate one step
    # further in order to translate it of one step
    # we generate and cut our intervals at the end of the generation
    right_interval <- seq(mean, bound + step, step)
    right_heights <- dnorm(right_interval, mean, sd)[seq_len(length(right_interval))]
    left_heights <- rev(right_heights)[seq_len(length(right_heights) - 1)]
    right_heights <- right_heights[seq(2, length(right_heights))]
    right_interval <- right_interval[seq_len(length(right_interval) - 1)]
    left_interval <- rev(-right_interval) + 2 * mean - step
  }
  intervals <- c(left_interval, right_interval)
  heights <- c(left_heights, right_heights)
  heights <- coef * heights
  return(list("dots" = intervals, "heights" = heights))
}

generate_boxes <- function(weights, means, sd, number_of_boxes) {
  lb <- qnorm(0.001, min(means), max(sd))
  rb <- qnorm(0.999, max(means), max(sd))
  size <- length(weights)
  list_of_normal_boxes <- lapply(1:size, function(i) {
    boxes(inf = lb, sup = rb, number_of_boxes = number_of_boxes, mean = means[i], sd = sd[i], coef = weights[i])
  })
  list_of_dots <- lapply(1:size, function(i) list_of_normal_boxes[[i]]$dots)
  list_of_heights <- lapply(1:size, function(i) list_of_normal_boxes[[i]]$heights)
  rank <- rep(0, size)
  dots <- Reduce(c, list_of_dots)
  dots <- dots[order(sapply(dots, function(x) x, simplify = TRUE))]
  number_of_dots <- length(dots)
  heights <- rep(0, number_of_dots)
  #print(length(list_of_dots[[1]]))
  for (i in 1:number_of_dots) {
    for (j in 1:size) {
      if (i > 1 && dots[[i]] != dots[[i - 1]] || i == 1) {
        # if the normal laws are close to each other they may have
        # the same dots. We ensure that there's no problem if this
        # is the case.
        if (dots[[i]] %in% list_of_dots[[j]]) {
          rank[[j]] <- rank[[j]] + 1
        }
      }
      if (rank[[j]] != 0) {
        rank_to_add <- rank[[j]]
        heights[[i]] <- heights[[i]] + list_of_heights[[j]][[rank_to_add]]
      }
    }
    heights[[i]] <- heights[[i]]+0.005
  }
  list("dots" = dots, "heights" = heights)
}

area <- function(dots, heights) {
  size <- length(dots)
  boxes <- lapply(2:size, function(i) {
    heights[i - 1] * (dots[i] - dots[i - 1])
  })
  boxes <- unlist(boxes)
  sum(boxes)
}

cumulative_area <- function(dots, heights) {
  size <- length(dots)
  areas <- rep(0,size-1)
  areas[1] <- heights[1] * (dots[2] - dots[1])
  for(i in seq(2,size-1,1)){
    areas[i] <- areas[i-1] + heights[i - 1] * (dots[i] - dots[i - 1])
  }
  areas
}

simulation <- function(weights, means, sd, number_of_boxes, number_of_variates) {
  mixture <- function(x) {
    normal_laws <- lapply(seq_len(length(weights)), function(i) {
      weights[i] * dnorm(x, mean = means[i], sd = sd[i])
    })
    Reduce("+", normal_laws)
  }
  boxes <- generate_boxes(weights, means, sd, number_of_boxes)
  dots <- boxes$dots
  heights <- boxes$heights
  area <- area(dots, heights)
  cumulative_area <- cumulative_area(dots, heights)
  adjusted_area <- cumulative_area / area
  
  uniform <- runif(number_of_variates)
  #we're selecting the box by stopping right after the last box
  #with a value smaller than our uniform
  sample <- sapply(uniform, function(choice){
    to_sample <- length(adjusted_area[adjusted_area<=choice])
    absc <- 0
    ord <- .Machine$integer.max
    while(mixture(absc)<ord){
      absc <- runif(1, min = dots[to_sample], max = dots[to_sample+1])
      ord <- runif(1, min=0, max = max(heights[to_sample],heights[to_sample+1]))
    }
    ord
  },simplify = TRUE)
}

moyenne <- 11
std <- 4
nb <- 10000
nv <- 10000
weights <- c(2, -1)
means <- c(2, 2)
size <- length(weights)
standard_deviations <- c(2, 1)
co <- 2
a <- boxes(inf = -1, sup = 1, mean = moyenne, number_of_boxes = nb, sd = std, coef = co)
x <- a$dots
y <- a$heights
plot(x, y, type = "s")
t <- co * dnorm(x, mean = moyenne, sd = std)
lines(x, t)
b <- generate_boxes(weights, means, standard_deviations, nb)
x1 <- b$dots
y1 <- b$heights
plot(x1, y1, type = "s")
z <- lapply(1:size, function(i) {
  weights[i] * dnorm(x1, means[i], standard_deviations[i])
})
z <- Reduce("+", z)
lines(x1, z)
#are <- cumulative_area(x1, y1)
a=simulation(weights,means,standard_deviations,nb,nv)