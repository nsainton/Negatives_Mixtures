


set.seed(18)

n <- 100

w=c(2,-1)
mean=c(2,2)
sd=c(2,1)


r1=qnorm(0.001,min(mean),max(sd))
r2=qnorm(0.999,max(mean),max(sd))
L=seq(r1,r2,(abs(r1-r2)/n))


cumu <- function(x){
  cumulative <-0
  for(i in 1:length(mean)){
    cumulative <- cumulative + w[i]*pnorm(x,mean[i],sd[i])
  }
  
  return(cumulative)
}

cumulate<-matrix(nrow=1,ncol=length(L))

cumulate<-lapply(L,cumu)





