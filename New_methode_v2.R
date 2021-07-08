set.seed(18)


cumu <- function(x, w, mean, sd) {
  normal_laws <- lapply(seq_len(length(w)), function(i) {
    w[i] * pnorm(x, mean = mean[i], sd = sd[i])
  })
  Reduce("+", normal_laws)
}
#Calcul de la fdr
fdr <- function(n=100, w, mean, sd) {
  r_1 <- qnorm(0.001, min(mean), max(sd))
  r_2 <- qnorm(0.999, max(mean), max(sd))
  k <- seq(r_1, r_2, (abs(r_1 - r_2) / n))
  for (i in length(w)){
    k<- c(k, mean[i])
  }
  k<-sort(k)
  cumulate <- matrix(nrow = 1, ncol = length(k))
  cumulate <- lapply(k, function(x) cumu(x, w, mean, sd))
  return(list(cumulate, k))
}

mixture <- function(x) {
  normal_laws <- lapply(seq_len(length(w)), function(i) {
    w[i] * dnorm(x, mean = mean[i], sd = sd[i])
  })
  Reduce("+", normal_laws)
}

simulate <- function(w, mean, sd, n, m) {
  
  liste <- fdr(m, w, mean, sd)[[2]]
  cumulate <- fdr(m, w, mean, sd)[[1]]
  u_1 <- runif(n)
  nb_to_generate <- rep(0, length(liste))
  nbGenere <- 0
  res <- rep(0, n)
  for (i in 1:(length(cumulate) - 1)) {
    nbGenere <- 0
    #Choisis selon quel rectangle de la fdr on va générer
    to_generate <- which(u_1 < cumulate[i + 1] & u_1 >= cumulate[i])
    nb_to_generate[i] <- length(to_generate)
    #Acceptation rejet
    while (nbGenere < nb_to_generate[i]) {
      #Une fois que l'on sait dans quel rectangle tirer Y, C'est une loi uniforme entre 2 points
      Y <- runif((nb_to_generate[i] - nbGenere), min = liste[i], max = liste[i + 1])
      u_2 <- runif((nb_to_generate[i] - nbGenere))
      height <- max(mixture(liste[i]), mixture(liste[i + 1]))
      acceptate <- which((height * u_2) <= mixture(Y))
      if (length(acceptate) > 0) {
        res[to_generate[(nbGenere + 1):min((nbGenere + length(acceptate)), nb_to_generate[i])]]=Y[acceptate[1 : min(nb_to_generate[i], length(acceptate))]]
        nbGenere=nbGenere + length(acceptate)
      }
    }
  }
  return(res)
}

w<- c(2,-1)
mean <- c(2,2)
sd <- c(2,1)
m<-10
n<-100

res <- simulate(w, mean, sd, n, m)
p1 <- hist(main = "Répartition des valeurs générées", xlab = "Valeurs", ylab = "Répartition", res, breaks = 100, freq = FALSE)
curve(2 * dnorm(x, 2, 2) - dnorm(x, 2, 1)+ dnorm(x,3,2), add = TRUE, col = "red")

microbenchmark::microbenchmark(simulate(w, mean, sd, n, m))
