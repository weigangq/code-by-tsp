#GOal-to compute the engineering approach statistic for measuring optimization. 
# to run the function you'll need the 4 sets of errors (pol, hydro, vol, iso) from your results. You'll also need 4 sets of random errors (pol, hydro, vol, iso) made from code-shuffle program   

### engineering approach function
engineerApproach(errordata,serrordata)



#Engineering Approach function


##mydata- a list made of your 4 sets of errors.
##randomdata- a list made of 4 sets of random errors

engineerApproach<- function(mydata,randomdata) {
  minerror<-list()
  
  for (i in 1:4) {
    minerror[[i]]<-min(mydata[[i]])
  }
  
  meanserror<-list()
  for (i in 1:4) {
    meanserror[[i]]<-mean(randomdata[[i]])
  }
  
  engineApproach<-list()
  for (i in 1:4) {
    engineApproach[[i]]<-(meanserror[[i]]-sgcerror[[i]])/(meanserror[[i]]-minerror[[i]])
  }
  return(engineApproach)
}
