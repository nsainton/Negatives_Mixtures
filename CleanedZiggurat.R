rm(list=ls())
#We begin by creating the two first function we need to generate our boxes, generate_boxes
#To do this, we're first going to create a function that generates the boxes on the decreasing
#part of each normal law, the function we'll then copy and paste the boxes on the increasing part
#to have a full approximation with boxes of our density.

boxes=function(left_bound=-1,right_bound=1,number_of_boxes=1,mean=0,sd=1,type="over"){
  const=sd*sqrt(2*pi)
  density=function(x){
    return ((1/const)*exp((-(x-mean)^2)/(2*sd^2)))
  }
  bound = ifelse(abs(mean-left_bound)>abs(mean-right_bound),mean+abs(mean-left_bound),right_bound)
  #We're just getting sure that the right bound is the more distant to the mean of our normal law between the left and the right bound
  step=(bound-mean)/number_of_boxes
  intervals=seq(mean,bound,step)
  if(type=="over"){
    heights=density(intervals)
  }else if(type=="under"){
    heights=density(seq(mean+step,bound+step,step))
  }
  #If we generate boxes under the density we need to go one step further and begin one step further to get sure the boxes stay under the density of our normal law
  return(list(intervals,heights))
}

symetry=function(intervals,heights,mean,bound){
  
}

a<-boxes(number_of_boxes = 10, type = "under")
x<-a[[1]]
y<-a[[2]]
plot(x,y,type='s')
lines(x,dnorm(x))