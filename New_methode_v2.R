
set.seed(18)


n <- 100
w <- c(2, -1)
mean <- c(2, 2)
sd <- c(2, 1)
 
cumu <- function(x, w=c(2, -1), mean=c(2, 2), sd=c(2, 1)) {
  cumulative <- 0
  for (i in 1:length(mean)) {
    cumulative <- cumulative + w[i] * pnorm(x, mean[i], sd[i])
  }
  return(cumulative)
}

#Calcul de la fdr
fdr <- function(n=100, w=c(2, -1), mean=c(2, 2), sd=c(2, 1)) {
  r1 <- qnorm(0.001, min(mean), max(sd))
  r2 <- qnorm(0.999, max(mean), max(sd))
  K <- seq(r1, r2, (abs(r1 - r2) / n))
  cumulate <- matrix(nrow = 1, ncol = length(K))
  cumulate <- lapply(K, cumu)
  return(list(cumulate, K))
}

mixture <- function(x) {
  normal_laws <- lapply(seq_len(length(w)), function(i) {
    w[i] * dnorm(x, mean = mean[i], sd = sd[i])
  })
  Reduce("+", normal_laws)
}

choice_h <- function(i) {
  x <- mixture(L[i])
  y <- mixture(L[i + 1])
  if (x > y) {
    return(list(x, y))
  }
  else (return(list(y, x)))
}


L <- fdr(n, w, mean, sd)[[2]]
cumulate <- fdr(n, w, mean, sd)[[1]]

simulate <- function(n = 1000) {
  
  U1 <- runif(n)
  nbAGenerer <- rep(0, length(L))
  nbGenere <- 0
  res <- rep(0, n)
  for (i in 1:(length(cumulate) - 1)) {
    nbGenere <- 0
    #Choisis selon quel rectangle de la fdr on va générer
    AGenerer <- which(U1 < cumulate[i + 1] & U1 >= cumulate[i])
    nbAGenerer[i] <- length(AGenerer)
    #Acceptation rejet
    while(nbGenere < nbAGenerer[i]) {
      #Une fois que l'on sait dans quel rectangle tirer Y, C'est une loi uniforme entre 2 points
      Y <- runif((nbAGenerer[i] - nbGenere), min = L[i], max = L[i + 1])
      U2 <- runif((nbAGenerer[i] - nbGenere))
      height <- choice_h(i)[[1]]
      Acceptes <- which((height * U2) <= mixture(Y))
      if (length(Acceptes) > 0) {
        res[AGenerer[(nbGenere+1):min((nbGenere+length(Acceptes)),nbAGenerer[i])]]=Y[Acceptes[1:min(nbAGenerer[i],length(Acceptes))]]
        nbGenere=nbGenere+length(Acceptes)
      }
    }
  }
  return(res)
}

res <- simulate(100)
p1 <- hist(main = "Répartition des valeurs générées", xlab = "Valeurs", ylab = "Répartition", res, breaks = 100, freq = FALSE)
curve(2 * dnorm(x, 2, 2) - dnorm(x, 2, 1), add = TRUE, col = "red")

microbenchmark::microbenchmark(simulate(1000))
