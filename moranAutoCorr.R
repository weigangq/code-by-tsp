####This file gives the code and the Moran I autocorrelation test for the tRNA tree.


library(ape)
library(ade4)

#####Section 1- tRNA tree and Plot

#t=tRNA tree
t<-read.tree('C:.../tRNA_mdpt.dnd')

#tRNA tree plot.
plot(t,font=1,use.edge.length=TRUE,no.margin=T,edge.width=3, align.tip.label=T)


######Section 2-  amino acids and values

#aa= amino acids file for parsing data
aa<-read.csv('...../aa.csv',sep=",",header = F)
#aa
aa3<-aa[-c(21,22),]  

#aadistances amino acids with physiochemical values. 
aadistances<- read.csv("..../aadistances.csv",header = TRUE)
#aadistances



###Section 3 - Moran's I test for autocorrelation.

##moran_prep- gives a vector consisting of the  measure of polarity, or another characteristic, for each tip label.  
moran_prep<-function(tiplab,measure) {
  tips<-c(rep(1,length(tiplab)))   #this vector will give the polarity values for each elt in t$tip.labels
  colm<-which(colnames(aadistances)==measure)
  
  for (i in 1:nrow(aa3)){
    aavector<-grep(as.character(aa3[i,1]),tiplab) #gives positions in  the tip.labels vector correspond to aa3[i,1]  appears
    initials<-as.character(aa[as.character(aa3[i,1]),'init'])  #gives the initial of the aa3[i,1]
    v<-which(as.character(aadistances[,1])==as.character(aa3[i,2])) #give the position of initial in the aadistances matrix
    meas<-aadistances[v,colm]      #gets polarity value for the set of positions from aavector. 
    for (j in 1:length(aavector)){
      tips[aavector[j]]<-meas
    }
  }
  return (tips)
}
moranI<-function(tree) {
  tiplab<-tree$tip.label     #t$tip.labels
  M_results<-list()          #list of results for each measure
  aad<-aadistances[,2:5]      
  
  for (i in colnames(aad)) {   
    tps<-moran_prep(tiplab, i)    #produces a vector of polarity, or other values,  corresponding to the tip labels
    df.tips<-data.frame(tiplab,tps)
    #w.con a matrix ith don't know how it's made.  Ask Weigang. 
    w.con <- 1/cophenetic.phylo(tree)   #weight matrix
    diag(w.con) <- 0
    w.con[w.con == Inf] <- 0
    m_result<-Moran.I(df.tips[,2],w.con)
    M_results[[length(M_results)+1]]<-m_result
  }
  names(M_results)<-colnames(aad)
  return(M_results)
}

#Moran I autocorrelation test on a tRNA tree, t. 
moranResults<-moranI(t)


###Moran I test results in table form 

read.csv(moran, "......../MoranIAutocorrelation.csv")
