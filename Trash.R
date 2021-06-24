#Ziggourat simulation.R
rm(list=ls())

#Nombre de boites du ziggourat sur la moiti? droite de chaques lois normales
numberOfBoxes=400
#MAtrice d?finissant le m?lange, la premi?re colonne repr?sente les coefficients, 
#la deuxi?me colonne repr?sente les moyennes et le troisi?me repr?sente les ?carts types
X=matrix(c(2,2,2,-1,2,1),ncol=3,byrow=TRUE)
Ziggourat=function(r1,r2,numberOfBoxes,mu=0,sigma=1){
  const=sigma*sqrt(2*pi)
  f<-function(x){
    return ((1/const)*exp((-(x-mu)^2)/(2*sigma^2)))
  }
  #G?nere le ziggourat sup?rieur ? la courbe sur la partie d?croissante jusqu'au point r qui est le plus distants du mode entre r1 et r2
  if(abs(mu-r1)>abs(mu-r2)){
    r=mu+abs(mu-r1)
  }
  else{
    r=r2
  }
  
  L=seq(mu,r,(r-mu)/numberOfBoxes)
  Hauteurs<-f(L)
  return(list(L,Hauteurs))
}
ZiggouratInscrit=function(r1,r2,numberOfBoxes,mu=0,sigma=1){
  #G?n?re un ziggourat sup?rieur ? la courbe
  result=Ziggourat(r1,r2,numberOfBoxes+1,mu,sigma)
  L=result[1][[1]]
  Hauteurs=result[2][[1]]
  #D?cale le ziggourat vers la gauche le placer sous la courbe
  L=L[2:(length(L)-1)]
  Hauteurs=Hauteurs[2:(length(Hauteurs))]
  
  #Copie le ziggourat par sym?trie sur la partie croissante
  L1=rev(mu-abs(mu-L[1:(length(L))]))
  cutGauche=which(L1>=r1)
  cutDroit=which(L<=r2)
  L=c(L1[cutGauche], L[cutDroit])
  Hauteurs=c(rev(Hauteurs[1:(length(Hauteurs)-1)])[cutGauche],Hauteurs[2:(length(Hauteurs))][cutDroit])
  
  return(list(L,-Hauteurs))
}
ZiggouratExte=function(r1,r2,numberOfBoxes,mu=0,sigma=1){
  result=Ziggourat(r1,r2,numberOfBoxes,mu,sigma)
  
  
  L=result[1][[1]]
  Hauteurs=result[2][[1]]
  #Recopie le ziggourat par sym?trie sur la partie croissante
  L1=rev(mu-abs(mu-L[2:(length(L))]))
  L=L[2:length(L)]
  cutGauche=which(L1>=r1)
  cutDroit=which(L<=r2)
  L=c(L1[cutGauche], L[cutDroit])
  Hauteurs=c(rev(Hauteurs[1:(length(Hauteurs)-1)])[cutGauche],Hauteurs[2:(length(Hauteurs))][cutDroit])
  return(list(L,Hauteurs))
}






tracer=function(X, hauteur){
  #Trace la courbe 2*N(2,2)-N(2,1) et le ziggourat d?finie par X et hauteur
  epsilon=0.00001
  Y=c()
  copie=c()
  Y[1]=X[1]
  for(i in 1:length(X)){
    Y[2*i]=X[i]
    Y[(2*i)+1]=X[i+1]-epsilon
    copie[2*i]=hauteur[i]
    copie[(2*i)+1]=hauteur[i]
    
  }
  plot(Y,copie,type="l")
  print((Y[length(Y)-1]-Y[2])/400)
  Y=seq(Y[2],Y[length(Y)-1], (Y[length(Y)-1]-Y[2])/400)
  lines(Y,2*dnorm(Y,2,2)-dnorm(Y,2,1))
  
}

points<-list()
hauteur=list()


r1=qnorm(0.001,min(X[,2]),max(X[,3]))
r2=qnorm(0.999,max(X[,2]),max(X[,3]))

#G?nere les ziggourats stock?s dans points et hauteurs ou points repr?sente les abcisses de changement de hauteur
#Et hauteur[[1]][i] repr?sente la hauteur du ziggourat entre points[[1]][i] et points[[1]][[i+1]]
for(i in 1:nrow(X)){
  if(X[i,1]>0){
    result=ZiggouratExte(r1,r2,numberOfBoxes,X[i,2],X[i,3])
    points<- c(points,list(result[1][[1]]))
    hauteur<- c(hauteur,list(result[2][[1]]*X[i,1]))
  }
  if(X[i,1]<0){
    result=ZiggouratInscrit(r1,r2,numberOfBoxes,X[i,2],X[i,3])
    points<- c(points,list(result[1][[1]]))
    hauteur<- c(hauteur,list(result[2][[1]]*-1*X[i,1]))
  }
  
}

L=c()
for(i in 1:length(points)){
  L=c(L,points[[i]])
}
L=sort(L)


H=rep(0,length(L))
Debut=rep(1:length(points))

#Additionne toutes les hauteurs des ziggourats entre L[i] et L[i+1] pour obtenir H[i]
for(i in 1:length(L)){
  for(j in 1:length(points)){
    for(z in Debut[j]:(length(hauteur[[j]])-1)){
      if(L[i]>=points[[j]][z] && L[i]<points[[j]][z+1]){
        H[i]=H[i]+hauteur[[j]][z]
        Debut[j]=z
        break
        
      }
      else if(L[i]> points[[j]][length(points[[j]])]){
        H[i]=H[i]+hauteur[[j]][z]
        Debut[j]=length(points[[j]])
        break
      }
    }
  }
}

tracer(L,H)



CalculAire=function(L,H){
  #Calcul l'aire sous le ziggourat
  Aire=0
  for(i in 1:(length(L)-1)){
    Aire=Aire+((L[i+1]-L[i])*H[i])
  }
  return(Aire)
}



CalculAireCumule=function(L,H){
  #Renvoie une liste de l'aire cumul? du ziggourat
  Aire=c(((L[2]-L[1])*H[1]))
  for(i in 2:(length(L)-1)){
    Aire=c(Aire,Aire[length(Aire)]+((L[i+1]-L[i])*H[i]))
  }
  
  return(Aire)
}



simulate=function(L,H,X,r1,r2,n=1000){
  f=function(x){
    y=0
    for(i in 1:nrow(X)){
      y=y+X[i,1]*dnorm(x,X[i,2],X[i,3])
    }
    return(y)
  }
  
  #Calcul l'aire Cumule et la normalise
  AireCumule=CalculAireCumule(L,H)
  Aire=CalculAire(L,H)
  AireCumuleAjuste=AireCumule/Aire
  
  
  U1=runif(n)
  nbAGenerer=rep(0,length(L))
  nbGenere=0
  res=rep(0,n)
  for(i in 1:(length(AireCumuleAjuste)-1)){
    nbGenere=0
    #Choisis selon quel rectangle du ziggourat on va g?n?rer
    AGenerer=which(U1<AireCumuleAjuste[i+1] & U1>=AireCumuleAjuste[i])
    nbAGenerer[i]=length(AGenerer)
    #Acceptation rejet
    while(nbGenere<nbAGenerer[i]){
      #Une fois que l'on sait dans quel rectangle tirer Y, C'est une loi uniforme entre 2 points
      Y=runif((nbAGenerer[i]-nbGenere)*1.1,min=L[i],max=L[i+1])
      U2=runif((nbAGenerer[i]-nbGenere)*1.1)
      Acceptes=which(H[i]*U2<=f(Y))
      if(length(Acceptes)>0){
        res[AGenerer[(nbGenere+1):min((nbGenere+length(Acceptes)),nbAGenerer[i])]]=Y[Acceptes[1:min(nbAGenerer[i],length(Acceptes))]]
        nbGenere=nbGenere+length(Acceptes)
      }
    }
  }
  return(res)
}
res=simulate(L,H,X,r1,r2,100000)
p1=hist(main="R?partition des valeurs g?n?r?es",xlab="Valeurs",ylab="R?partition",res,breaks=100,freq=FALSE)
curve(2*dnorm(x,2,2)-dnorm(x,2,1),add=TRUE,col="red")

#Code en bordel
rm(list=ls())
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
    #left_interval <- c()
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
    right_heights <- dnorm(right_interval, mean, sd)
    left_heights <- rev(right_heights)[seq_len(length(right_heights) - 1)]
    right_heights <- right_heights[seq(2, length(right_heights))]
    right_interval <- right_interval[seq_len(length(right_interval) - 1)]
    left_interval <- rev(-right_interval) + 2 * mean
  }
  intervals <- c(left_interval, right_interval)
  heights <- c(left_heights, right_heights)
  #heights <- coef * heights
  return(list("dots" = intervals, "heights" = heights))
}

moyenne <- 11
std <- 4
nb <- 10
weights <- c(2, -1)
means <- c(3, 2)
size <- length(weights)
standard_deviations <- c(2, 2)
co <- -2
a <- boxes(inf = -1, sup = 1, mean = moyenne, number_of_boxes = nb, sd = std, coef = co)
x <- a$dots
y <- a$heights
plot(x, y, type = "s")
t <- dnorm(x, mean = moyenne, sd = std)
lines(x, t)

to_sample <- length(adjusted_area[adjusted_area<=uniform])
print(to_sample)
absc <- 0
ord <- .Machine$integer.max
for(choice in uniform){
  to_sample <- length(adjusted_area[adjusted_area<=uniform])
  absc <- 0
  ord <- .Machine$integer.max
  while(mixture(absc)<ord){
    absc <- runif(1, min = dots[to_sample], max = dots[to_sample+1])
    ord <- runif(1, min=0, max = max(heights[to_sample],heights[to_sample+1]))
    print("reject")
    print(mixture(absc))
    print(ord)
  }
}

generate_boxes = function(weights,means,sd,numberofboxes){
  lb = qnorm(0.001,min(means),max(sd))
  rb = qnorm(0.999,max(means),max(sd))
  size = length(weights)
  list_of_normal_boxes = lapply(1:size, function(i){
    boxes(inf = lb, sup = rb, number_of_boxes = numberofboxes, mean=means[i], sd=standard_deviations[i],coef = weights[i])
  })
  list_of_dots = lapply(1:size, function(i) list_of_normal_boxes[[i]]$intervals)
  list_of_heights = lapply(1:size, function(i) list_of_normal_boxes[[i]]$heights)
  #In these first lines, we generate the boxes for all the normal laws of the mixture and extract the corresponding heights and dots
  size_list = length(list_of_dots[[1]])
  indexed_dots = lapply(1:size, function(i){
    list_of_dots[[i]] = lapply(1:size_list, function(j){
      c(list_of_dots[[i]][[j]],round(i,0),round(j,0))
    })
  })
  indexed_heights = lapply(1:size, function(i){
    list_of_heights[[i]] = lapply(1:size_list, function(j){
      c(list_of_heights[[i]][[j]],round(i,0),round(j,0))
    })
  })
  #In these two variables we indice the values with two number. The first one represents the position of the normal law in the mixture
  #and the second one represents the position of the value in the list
  dots = Reduce(c,indexed_dots)
  dots = dots[order(sapply(dots,function(x) x[1],simplify = TRUE))]
  #We now put all the indexed dots in a unique list and are going to put all the summed heights in another one
  count = rep(1,size_list)
  heights = lapply(1:length(dots),function(i){
    if(i==1){
      summ=sapply(1:size,function(j) indexed_heights[[i]][[1]],simplify = TRUE)
      sum(summ)
    }else{
      summ=sapply(1:size, function(j){
        if(dots[[i]][[2]]==j) c=1 #empty line just to not forget to continue there
      })
    }
  })
  indexed_heights
}

generate_boxes = function(weights, means, sd, number_of_boxes) {
  lb = qnorm(0.001,min(means),max(sd))
  rb = qnorm(0.999,max(means),max(sd))
  size = length(weights)
  list_of_normal_boxes = lapply(1:size, function(i){
    boxes(inf = lb, sup = rb, number_of_boxes = number_of_boxes, mean = means[i], sd = sd[i], coef = weights[i])
  })
  list_of_dots = lapply(1:size, function(i) list_of_normal_boxes[[i]]$dots)
  list_of_heights = lapply(1:size, function(i) list_of_normal_boxes[[i]]$heights)
  rank = rep(0,size)
  dots = Reduce(c,list_of_dots)
  dots = dots[order(sapply(dots,function(x) x,simplify = TRUE))]
  number_of_dots = length(dots)
  heights = rep(0, number_of_dots)
  for(i in 1:number_of_dots){
    for(j in 1:size){
      if(i>1 && dots[[i]]!=dots[[i-1]] || i==1) {
        if (dots[[i]] %in% list_of_dots[[j]]) {
          rank[[j]] = rank[[j]] + 1
        }
      }
      if (rank[[j]] != 0) {
        rank_to_add = rank[[j]]
        heights[[i]] = heights[[i]] + list_of_heights[[j]][[rank_to_add]]
      }
    }
    if(heights[[i]]<dnorm(dots[i]))
  }
  list("dots" = dots, "heights" = heights)
}

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

#rev_vs_sort.R
test=rep(0,2)
for(i in 1:1000){
  x=1:1000000
  x=dnorm(x)
  start_time = Sys.time()
  d=rev(x)
  end_time = Sys.time()
  
  #print("rev time")
  rev_time=end_time - start_time
  #print(end_time - start_time)
  
  start_time = Sys.time()
  d=sort(x)
  end_time = Sys.time()
  
  #print("sort time")
  sort_time = end_time - start_time
  #print(end_time - start_time)
  
  #print(sort_time<rev_time)
  i = ifelse(sort_time<rev_time, 2,1)
  test[i]=test[i]+1
}
print("sort")
print(test[2])
print("rev")
print(test[1])

#Mélanges finis alternés.R
#Ce programme simule le m?lange d'un m?me nombre de lois normales positives et n?gatives 
#dans le cas ou il est possible de les appareiller 2 ? 2 tel que a_i*f_i >=b_i*g_i


#Matrice de param?tres avec a chaque ligne une nouvelle loi. 
#La colonne 1 repr?sente les moyennes, la 2 les ecarts types et la 3 les coefficients du m?lange. Les n/2 premi?res lois sont celles de coefficients positif et les suivantes de coefficient n?gatif
x<- matrix(c(1,1,0.5,0,1,1.7,0,1,-0.7,1,1,-0.5), ncol=3, byrow=TRUE)
#Nombre de variable ? simuler
n=1000

positive <- function(mu_1, sigma_1, alpha, mu_2, sigma_2, beta){
  
  if (sigma_1 != sigma_2){
    m<- -(sigma_2^2*mu_1 - sigma_1^2*mu_2)/(sigma_2^2-sigma_1^2)
    result<- (((exp((m-mu_1)^2/(2*sigma_1^2))-exp((m-mu_2)^2/(2*sigma_2^2)))<= (sigma_2*alpha)/(sigma_1*(-beta))) & (sigma_1 > sigma_2))
  }
  
  if (sigma_1 == sigma_2){
    result <- (mu_1 == mu_2)
    
    
  }
  
  return(result)
}

costMatrix <- function(mixture){
  
  n <- nrow(mixture)
  accept <- matrix(ncol=(n/2), nrow=(n/2))
  
  for (i in 1:(n/2)){
    
    for (j in ((n/2)+1):n){
      print(mixture)
      if(positive(mixture[i,1], mixture[i,2], mixture[i,3], mixture[j,1], mixture[j,2], mixture[j,3])==TRUE){
        
        accept[i,j-(n/2)] <- 1
      }
      else {
        
        accept[i,j-(n/2)] <- 0
      }
      
      
    }
  }
  
  
  return(accept) 
}

cost=costMatrix(x)
print(cost)



##Implementation de la Methode Hongroise par sergioreyblanco: https://github.com/sergioreyblanco/hungarian_method
solutionChecking <- function(matrix, solution, costMatrix){
  solutionAux <- solution
  foundSolutions <- list()
  flag_found = FALSE
  counter = 0
  previousDuplicate = -1
  duplicate = anyDuplicated(solutionAux[seq(2, 2*nrow(matrix), by = 2)])
  
  while(duplicate != 0){
    
    for(j in 1:ncol(matrix)){
      if(matrix[duplicate, j] == 0){
        counter = counter + 1
      }
      
      if(matrix[duplicate, j] == 0 && j != solutionAux[2*duplicate]){
        
        s<-solutionAux
        s[2*duplicate] = j
        flag_found = FALSE
        for(p in foundSolutions){
          paux <- as.vector(unlist(p))
          if(all(s == paux) == TRUE){
            flag_found = TRUE
            break
          }
        }
        if(flag_found == FALSE){
          solutionAux <- s
          
          foundSolutions[[length(foundSolutions)+1]] <- list(solutionAux)
          
          break
        }
      } 
      
    }
    
    if(anyDuplicated(solutionAux[seq(2, 2*nrow(matrix), by = 2)]) == 0){
      break
    }
    
    if(counter == 1){
      flag_found = TRUE
    }
    counter = 0
    
    if(flag_found == TRUE){
      
      for(i in 1:nrow(matrix)){
        for(j in 1:ncol(matrix)){
          
          if(matrix[i, j] == 0){
            s<-solutionAux
            s[2*(i-1)+1] <- i
            s[2*(i-1)+2] <- j
            flag_found = FALSE
            for(p in foundSolutions){
              paux<-as.vector(unlist(p))
              if(all(s == paux) == TRUE){
                flag_found = TRUE
                break
              }
            }
            if(flag_found == FALSE){
              solutionAux <- s
              break
            }
          }
        }
        if(flag_found == FALSE){
          break
        }
      }
      foundSolutions[[length(foundSolutions)+1]] <- list(solutionAux)
    }
    
    previousDuplicate = duplicate
    
    duplicate = anyDuplicated(solutionAux[seq(2, 2*nrow(matrix), by = 2)])
    
    if(previousDuplicate == duplicate){
      duplicate <- which(duplicated(solutionAux[seq(2, 2*nrow(matrix), by = 2)], fromLast = TRUE)==TRUE)[1]
    }
  }
  
  solutionValue = 0
  for(k in 1:nrow(matrix)){
    if(costMatrix[solutionAux[2*k-1],solutionAux[2*k]] != M){
      solutionValue = solutionValue + costMatrix[solutionAux[2*k-1],solutionAux[2*k]]
    }  
  }
  
  return(list("solution" = solutionAux, "cost" = solutionValue))
}

M <- 10000
Hungarian <- function(costMatrix){
  A <- costMatrix
  
  #Step 1: calculation fo reduced matrix
  
  #the minimum of each row is calculated and they are subtracted
  min = 10*M
  for(i in 1:nrow(A)){
    for(j in 1:ncol(A)){
      if(A[i, j] < min){
        min = A[i, j]
      }
    }
    A[i,] <- A[i,]-min
    
    min = 10*M
  }
  
  #the minimum of each column is calculated and they are subtracted
  min = 10*M
  for(j in 1:ncol(A)){
    for(i in 1:nrow(A)){
      if(A[i, j] < min){
        min = A[i, j]
      }
    }
    A[,j] <- A[,j]-min
    
    min = 10*M
  }
  
  while(TRUE){  
    #Step 2: zero strikethrough cross out  
    assignedRows <- list()
    assignedZeros <- list()
    alreadyCrossout <- list()
    for(i in 1:nrow(A)){
      for(j in 1:ncol(A)){
        if(A[i, j] == 0 && ((i-1)*ncol(A)+j) %in% alreadyCrossout == FALSE){
          assignedZeros <- c(assignedZeros, ((i-1)*ncol(A)+j))
          assignedRows <- c(assignedRows, i)
          if((i+1) < nrow(A)){
            for(k in (i+1):nrow(A)){
              if(A[k, j] == 0  && ((k-1)*ncol(A)+j) %in% alreadyCrossout == FALSE){
                alreadyCrossout<-c(alreadyCrossout, ((k-1)*ncol(A)+j))
              }
            }
          }
          
          break;
        }
      }
    }
    
    rowMarks <- list()
    columnMarks <- list()
    for(i in 1:nrow(A)){
      if(i %in% assignedRows == FALSE){
        rowMarks <- c(rowMarks, i)
      }
    }
    
    for(i in rowMarks){
      for(j in 1:ncol(A)){
        if(A[i, j] == 0 && j %in% columnMarks == FALSE){
          columnMarks <- c(columnMarks, j)
        }
      }
    }
    for(i in columnMarks){
      for(j in 1:nrow(A)){
        if(((j-1)*nrow(A)+i) %in% assignedZeros && j %in% rowMarks == FALSE){
          rowMarks <- c(rowMarks, j)
        }
      }
    }
    
    rowLines <- list()
    columnLines <- list()
    columnLines <- columnMarks
    for(i in 1:nrow(A)){
      if(i %in% rowMarks == FALSE){
        rowLines <- c(rowLines, i)
      }
    }
    
    if(length(rowLines) + length(columnLines) == nrow(A)){
      break;
    }else{
      #Step 3: reduced matrix is updated
      min = M
      for(i in 1:nrow(A)){
        for(j in 1:ncol(A)){
          if(A[i, j] < min && i %in% rowLines == FALSE && j %in% columnLines == FALSE){
            min = A[i, j]
          }
        }
      }
      for(i in 1:nrow(A)){
        for(j in 1:ncol(A)){
          if(i %in% rowLines == FALSE && j %in% columnLines == FALSE){
            A[i, j] <- A[i, j] - min
          } else if(i %in% rowLines == TRUE && j %in% columnLines == TRUE){
            A[i, j] <- A[i, j] + min
          } else if(i %in% rowLines == TRUE && j %in% columnLines == FALSE){
            A[i, j] <- A[i, j]
          } else if(i %in% rowLines == FALSE && j %in% columnLines == TRUE){
            A[i, j] <- A[i, j]
          }
        }
      }
    }
  }  
  
  #Step 4: selection of the solution
  solution <- c()
  for(i in 1:nrow(A)){
    for(j in 1:ncol(A)){
      if(A[i, j] == 0){
        solution <- c(solution, i, j)
        break
      }
    }
  }
  
  list <- solutionChecking(A, solution, costMatrix)
  
  return(list)
}
res=Hungarian(-cost)
assignment=c()
for(i in 1:(length(res$solution)/2)){
  assignment=append(assignment,res$solution[i*2]+nrow(x)/2)
}
print(assignment)
print(res)




rejection <- function( M, g,i,j,n,x) {
  naccepts <- 0
  result.sample <- rep(NA, n)
  
  while (naccepts < n) {
    y <- rnorm(1,x[i,1],x[i,2])
    u <- runif(1)
    
    if ( u <= x[i,3]*dnorm(y,x[i,1],x[i,2])+x[j,3]*dnorm(y,x[j,1],x[j,2]) / (M*dnorm(y,x[i,1],x[i,2])) ) {
      naccepts <- naccepts + 1
      result.sample[naccepts] = y
    }
  }
  
  result.sample
}

negative_mixture_sampling<-function(x,n, assignment){
  coeff=c(0)
  result=vector(,n)
  for(i in 1:(nrow(x)/2)){
    j=assignment[i]
    coeff=append(coeff,coeff[length(coeff)]+x[i,3]+x[j,3])#compute the cumulative sum of alpha-beta for each couple, it will be used as weights
  }
  U=runif(n)
  for(i in 1:(length(coeff)-1)){
    to_simulate=which(U >= coeff[i] & U<coeff[i+1])
    number_simulation=length(to_simulate)
    j=assignment[i]
    
    simulated=rejection(x[i,3],g,i,j,number_simulation,x)
    result[to_simulate]=simulated
  }
  
  return(result)
  
}
negative_mixture_sampling(x,n,assignment)