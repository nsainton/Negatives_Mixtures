

set.seed(18)

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


