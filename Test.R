source("https://raw.githubusercontent.com/poulem/Negatives_Mixtures/main/CleanedZiggurat.R")
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
a <- simulation(weights, means, standard_deviations, nb, nv)
hist(main = "Distribution of generated values", xlab = "Values", ylab = "Distribution", a, breaks = 100, freq = FALSE)
