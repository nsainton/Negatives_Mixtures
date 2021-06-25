rm(list=ls())



set.seed(180)

cumu <- function(x,w=c(2,-1),mean=c(2,2),sd=c(2,1)){
  cumulative <-0
  for(i in 1:length(mean)){
    cumulative <- cumulative + w[i]*pnorm(x,mean[i],sd[i])
  }
  
  return(cumulative)
}

fdr<- function(n=100,w=c(2,-1),mean=c(2,2),sd=c(2,1)){
  
  r1=qnorm(0.001,min(mean),max(sd))
  r2=qnorm(0.999,max(mean),max(sd))
  K=seq(r1,r2,(abs(r1-r2)/n))
  
  cumulate<-matrix(nrow=1,ncol=length(K))
  
  cumulate<-lapply(K,cumu)
  return(list(cumulate,K))
}

#Nombre de boites du ziggourat sur la moitié droite de chaques lois normales
numberOfBoxes=4
#MAtrice définissant le mélange, la première colonne représente les coefficients, 
#la deuxième colonne représente les moyennes et le troisième représente les écarts types
X=matrix(c(2,2,2,-1,2,1),ncol=3,byrow=TRUE)
Ziggourat=function(r1,r2,numberOfBoxes,mu=0,sigma=1){
  const=sigma*sqrt(2*pi)
  f<-function(x){
    return ((1/const)*exp((-(x-mu)^2)/(2*sigma^2)))
  }
  #Génere le ziggourat supérieur à la courbe sur la partie décroissante jusqu'au point r qui est le plus distants du mode entre r1 et r2
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
  #Génère un ziggourat supérieur à la courbe
  result=Ziggourat(r1,r2,numberOfBoxes+1,mu,sigma)
  L=result[1][[1]]
  Hauteurs=result[2][[1]]
  #Décale le ziggourat vers la gauche le placer sous la courbe
  L=L[2:(length(L)-1)]
  Hauteurs=Hauteurs[2:(length(Hauteurs))]
  
  #Copie le ziggourat par symétrie sur la partie croissante
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
  #Recopie le ziggourat par symétrie sur la partie croissante
  L1=rev(mu-abs(mu-L[2:(length(L))]))
  L=L[2:length(L)]
  cutGauche=which(L1>=r1)
  cutDroit=which(L<=r2)
  L=c(L1[cutGauche], L[cutDroit])
  Hauteurs=c(rev(Hauteurs[1:(length(Hauteurs)-1)])[cutGauche],Hauteurs[2:(length(Hauteurs))][cutDroit])
  return(list(L,Hauteurs))
}






tracer=function(X, hauteur){
  #Trace la courbe 2*N(2,2)-N(2,1) et le ziggourat définie par X et hauteur
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


#Génere les ziggourats stockés dans points et hauteurs ou points représente les abcisses de changement de hauteur
#Et hauteur[[1]][i] représente la hauteur du ziggourat entre points[[1]][i] et points[[1]][[i+1]]

generate_sum <-function(r1,r2){
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

return(list(L,H))
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
  #Renvoie une liste de l'aire cumulé du ziggourat
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
  
  res_cumul <-c()
  nbAGenerer=rep(0,1600)
  nbGenere=0
  res=0
  U1=runif(n)
  for(l in 1:n){
    for(i in 1:(length(fdr()[[2]])-2)){
    
    if(U1[l]<fdr()[[1]][i+1] & U1[l]>=fdr()[[1]][i]){
    r1 = fdr()[[2]][i]
    r2 = fdr()[[2]][i+1]
    
  #Calcul l'aire Cumule et la normalise
  Ziggusum <- generate_sum(r1,r2)
  AireCumule=CalculAireCumule(Ziggusum[[1]],Ziggusum[[2]])
  Aire=CalculAire(Ziggusum[[1]],Ziggusum[[2]])
  AireCumuleAjuste=AireCumule/Aire
  U2=runif(1)
  
  for(k in 1:(length(AireCumuleAjuste)-1)){
    nbGenere=0
    #Choisis selon quel rectangle du ziggourat on va générer
    if(U2<AireCumuleAjuste[k+1] & U2>=AireCumuleAjuste[k]){
    #Acceptation rejet
    while(nbGenere<1){
      #Une fois que l'on sait dans quel rectangle tirer Y, C'est une loi uniforme entre 2 points
      Y=runif(1,min=Ziggusum[[1]][k],max=Ziggusum[[1]][k+1])
      U3=runif(1)
      if(Ziggusum[[2]][k]*U3<=f(Y)){
        res=Y
        res_cumul<- c(res_cumul,res)
        nbGenere=nbGenere+1
      }
    }
    }
    }
  }
    }
    }
  print(length(res_cumul))
  return(res_cumul)
}

res=simulate(L,H,X,r1,r2,100)
p1=hist(main="Répartition des valeurs générées",xlab="Valeurs",ylab="Répartition",res,breaks=100,freq=FALSE)
curve(2*dnorm(x,2,2)-dnorm(x,2,1),add=TRUE,col="red")

#microbenchmark::microbenchmark(simulate(L,H,X,r1,r2,10))