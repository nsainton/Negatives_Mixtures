rm(list=ls())
#We begin by creating the first function we need to generate our boxes that will be called boxes.
#This function we'll take the bounds we need as an input, the number of boxes, the mean and the
#standard deviation of our normal law and the type of boxes we need (over or under the curve)
#And will generate two times the number of boxes to approximate the curve with

boxes=function(left_bound=-1,right_bound=1,number_of_boxes=1,mean=0,sd=1,coefficient=1){
  const=sd*sqrt(2*pi)
  density=function(x){
    return ((1/const)*exp((-(x-mean)^2)/(2*sd^2)))
  }
  bound = ifelse(abs(mean-left_bound)>abs(mean-right_bound),mean+abs(mean-left_bound),right_bound)
  #print(bound)
  #We're just getting sure that the right bound is the more distant to the mean of our normal law between the left and the right bound
  step=(bound-mean)/number_of_boxes
  if(coefficient>0){
    #If the type is "over" we generate one step further to go one step further on the other side by symmetry
    #We define the intervals accordingly and generate the heights on the giver intervals
    right_interval=seq(mean,bound+step,step)
    left_interval=sort(-right_interval)+2*mean
    length_right = length(right_interval)-1
    left_interval=left_interval[1:length_right]
    right_interval=right_interval[1:length_right]
    right_heights=density(right_interval)
    left_heights=sort(right_heights)[1:length_right]
  }else if(coefficient<0){
    #If the type is "under" we also generate one step further in order to translate it of one step
    #we generate and cut our intervals at the end of the generation
    right_interval=seq(mean,bound+step,step)
    right_heights=density(right_interval)[1:length(right_interval)]
    left_heights=sort(right_heights)[1:length(right_heights)-1]
    right_heights=right_heights[2:length(right_heights)]
    right_interval=right_interval[1:length(right_interval)-1]
    left_interval=sort(-right_interval)+2*mean-step
  }
  intervals=c(left_interval,right_interval)
  heights=c(left_heights,right_heights)
  return(list(intervals,coefficient*heights))
}

generate_boxes = function(weights,means,standard_deviations,numberofboxes=5){
  lb = qnorm(0.001,min(means),max(standard_deviations))
  rb = qnorm(0.999,max(means),max(standard_deviations))
  list_of_normal_boxes = lapply(1:length(weights), function(i){
    print(i)
    boxes(left_bound = lb, right_bound = rb, number_of_boxes = numberofboxes, mean=means[i], sd=standard_deviations[i],coefficient = weights[i])
  })
  list_of_normal_boxes
}

moyenne=12
std=5
weights=c(2,-1)
means=c(2,3)
standard_deviations=c(2,5)
a=boxes(mean=moyenne,number_of_boxes = 20,sd=std,coefficient=-2)
x=a[[1]]
y=a[[2]]
plot(x,y,type='s')
lines(x,-2*dnorm(x,mean=moyenne,sd=std))
b=generate_boxes(weights,means,standard_deviations)
