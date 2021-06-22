source("/Users/noahsaintonge/Documents/M1/Mémoire/Juin/Negative_Mixtures/CleanedZiggurat.R")
moyenne <- 11
std <- 4
nb <- 10
weights <- c(2, -1)
means <- c(3, 2)
size <- length(weights)
standard_deviations <- c(2, 2)
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
z1 <- lapply(1:size, function(i) {
  weights[i] * dnorm(x1, means[i], standard_deviations[i])
})

z <- z1[[1]] + z1[[2]]
lines(x1, z)