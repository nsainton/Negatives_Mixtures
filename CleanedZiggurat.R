rm(list=ls())
#We begin by creating the two first function we need to generate our boxes, generate_boxes
#To do this, we're first going to create a function that generates the boxes on the decreasing
#part of each normal law, the function we'll then copy and paste the boxes on the increasing part
#to have a full approximation with boxes of our density.

boxes=function(left_bound=-1,right_bound=1,number_of_boxes=1,mean=0,sd=1){
  const=sd*sqrt(2*pi)
  density=function(x){
    return ((1/const)*exp((-(x-mean)^2)/(2*sd^2)))
  }
  bound = ifelse(abs(mean-left_bound)>abs(mean-right_bound),mean+abs(mean-left_bound),right_bound)
  #We're just getting sure that the right bound is the more distant to the mean of our normal law between the left and the right bound
  intervals=seq(mean,bound,(bound-mean)/numberofboxes)
  heights=density(intervals)
  return(list(intervals,heights))
}

symetry=function(intervals,heights,mean,type="out"){
  
}