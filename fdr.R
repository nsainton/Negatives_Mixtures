


n<-100

f<-function(x,mu,sigma){
  const=sigma*sqrt(2*pi)
  return ((1/const)*exp((-(x-mu)^2)/(2*sigma^2)))
}

mixture_density <- function(x,weight=c(2,-1),mean=c(2,2),variance=c(2,1)){
  mix <-0
  for(i in 1:length(weight)){
      mix<-mix + weight[i]*f(x,mean[i],variance[i])
  }
  
  return(mix)
    
}

sample<- sort(2*rnorm(n,2,2)-rnorm(n,2,1))
cumulative<-c()
for(i in 1:n){
  cumulative<-c(cumulative,integrate(mixture_density,-Inf,sample[i])$value)
}





