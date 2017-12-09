#Goal: USE Van der Bout algorithm To find an path through the ammino acids using two of more of the following coodinates as coordinates (polarity, hydropathy, volume, isoelectricity).

#The TSP algorithm is found in Miller Bout paper. 


#Sections of code are organized by Roman Numerals.  Subsections of code are organized by Capital Letters


####Section I- Data Wrangling- Construct a distance Matrix Using W Qiu's normalized amino acid measurements.  

###I-A
##Import Dr Qiu's normalized amino acid measurements (aadistances). Make matrix (aadist) that includes columns of aadistances to be optimized, ex (pol, vol).

aadistances<- read.csv("C:/Users/Brian Sulkow/Documents/MonteCarloClub/Neural Nets/TSP/aadistances.csv",header = TRUE)
aadistances

##boxplot of 4 dimensions.  Refer to 12/2/17 for discussion. 
boxplot(aadistances[,-1], xlabs=colnames(aadistances[,-1]), main="Boxplot of Four Types of Amino Acid Measures")

## aadist- new df containing columns of aadistances you want to optimize.
aadist<-aadistances[,-1]
aadist


aadist<-aadist[,-4]
aadist<-aadist[,-2]
colnames(aadist)


###I-B Make the Distance matrix (NormaaMatrix) for TSP hopfield network

## aaMatrix- Use aadist for distances
aaMatrix<-dist(aadist, method="euclidean")


##NormaaMatrix-The matrix used for analysis. max value=1.

NormaaMatrix<-aaMatrix/max(aaMatrix)
NormaaMatrix<-as.matrix(NormaaMatrix)

#STOP Section I

####Section II
###In section II we make matrices that will be used for tsp.  This approach, of making matrices, differs from Olivers approach but is equivalent. 


###Section (II-A)
##Filter cols- This matrix is associated with summation that has the dp parameter in the Exi and energy terms of the algorithm.

#Size of the matrix is 21*21. This is the size of the set of all pairs (amino acid, time). 

filtercols<-matrix(c(rep(0,(21*21)^2)),nrow=21*21)
for (i in 1:21){
  for (j in 1:21){
    for(k in 1:21) {
      if(i!=k) { 
        filtercols[[21*(i-1)+j,21*(k-1)+j]]<-1
      }
    }
  }
}

#check of values of matrix. Should be 1's
filtercols[1,]
for(i in 1:21) {
  print(filtercols[23,21*(i-1)+2])
}



###Section (II-B)

##distpen-This Matrix corresponds to the second summation in the EXi and Energy equations. It is the term where distance appears.   
#(I'll provide a lower dimensional example of this matrix.  for you and Oliver to look at in the next day or two.)



distpen<-matrix(c(rep(0,(21*21)^2)),nrow=21*21)

#the first and last for loops in this construction ensure that you won't have a cyclic path, rather there will be distinct start and end points. 

for(i in 1:21) {
  for(k in 1:21) {
    distpen[[21*(i-1)+1,21*(k-1)+2]]<-NormaaMatrix[i,k]
  }
}

for(i in 1:21) {
  for(l in 2:20) {
    for(k in 1:21) {
      
      distpen[[21*(i-1)+l,21*(k-1)+(l-1)]]<-NormaaMatrix[i,k]
      distpen[[21*(i-1)+l,21*(k-1)+l+1]]<-NormaaMatrix[i,k]
    }
  }
}

for(i in 1:21) {
  for(k in 1:21) {
    distpen[[21*(i-1)+21,21*(k-1)+20]]<-NormaaMatrix[i,k]
  }
}

##STOP Section II


####Section III


###Section III-A
## Initial Conditions

#V-initial condition- The initial condition is a vector V of size 21*21, one entry for each pair (amino acid, time). Hopfield and Tank ('85 pg147) recommend that the initial condition have random values close to .5

V<-c(rep(1,21*21))
for (i in 1:(21*21)) {
  V[i]<-.5+runif(1,-.1,.1)  #random entry around .5
}
length(V)
V



### Section III-B
## Vectors and Dataframes to keep track of Energy, valid paths.  
energytemps<-c()
voltageout<-data.frame (0,0,0)
voltageout
colnames(voltageout)<-c("temperature","dparameter","Valid")
rvoltsls<-list()

###Section III-C
##THis section give a program for the algorithm by Miller, and van der Bout p 131.  ## Some comments on the program will refer to specific lines from the Miller, Van der Bout algorithm. Those  comments will start with the word "LINE".

#Note- there are TWO outer loops. These are not part of the algorithm on pg 131.  The  loop using "m" let's me run the algorithm for various values of the dp parameter (most outerloop). The loop using "k" lets me run the algorithm for various values of T. That's related to the annealing.  Rather than taking the temperature down exponentially as Oliver's algorithm does, the temperature is taken down linearly.  
#Note- Since I found values that worked well, I held dp and T constant in the program below . 
# With dp=.7 and T=.1 I was able to get 52 valid paths from 1500 runs (dim's were pol and volume). If you run this algorithm, try 50 runs.  You'll get 1 or 2 valid paths.

for (m in 0:1500)
{
  for (k in 1:1) {
    T<-.1 #-k*1e-1
    dp<-.7 # .1*(m-1)
    voltvector<-V     #voltvector =initial condition V
    energyresults<-list() #grabs info on energy.
    j<-700    #parameter of while loop given below  
    samplepoint<-c()   # vector to catch the amino acids that were randomly sampled used.
    while(j!=0) {  #LINE- the while loop is used for "do until (a fixed point is found)"
      samplecodon<-sample(c(1:21),1) #LINE-"Select a city at random"
      samplepoint[length(samplepoint)+1]<-samplecodon
      meanfield<-c(rep(1,21))
      sumterms<-0
      for(i in 1:21){  #LINE-for (i in i<=n)..
        meanfield[i]<-dp*filtercols[21*(samplecodon-1)+i,]%*%voltvector +
          distpen[21*(samplecodon-1)+i,]%*%voltvector #LINE- EXi equation.
        sumterms<-sumterms+exp(-meanfield[i]/T)  #LINE- "sum<=sum +e^-EXi/T"
        
      }
      
      for(i in 1:21) { ##LINE=for(i<=1;...) 
        voltvector[21*(samplecodon-1)+i]<- exp(-meanfield[i]/T)/sumterms
      }  #LINE-     "VXi=exp^(-EXi/T)/sum" This loop adjusts the voltvector, changing values corresponding to city=X time=i for all i
    #end loop  
      
      ##sumtest is meant to break the loop if any voltage goes to infinity.
      sumtest<-c(rep(1,21))
      for(i in 1:21){
        sumtest[i]<-voltvector[21*(samplecodon-1)+i]
      }
      
      st<-sum(sumtest)  
      if (is.nan(st)==TRUE) {
        print("bad")
        break}
    #LINE- "E=dp/2...."  THis is the energy equation  
      energy<-(dp/2)*voltvector%*%filtercols%*%voltvector+
        .5*voltvector%*%distpen%*%voltvector
      
      energyresults[[length(energyresults)+1]]<-energy
      j<-j-1
      
    }
    #End While Loop
    
    
    #the next few lines of the program is help guard against non convergence.  
    energytemps[length(energytemps)+1]<-energy
    if(is.nan(st)==TRUE) 
    { voltageout<-rbind(voltageout,c(T,dp,"NaN"))
    
    rvoltsls[[length(rvoltsls)+1]]<-"noconvergence" 
    next }
    volts<-t(matrix(voltvector,nrow=21))  # volts- is the matrix form of  voltvector. Notice we have to take the transpose. 
    
    #rvolts will create our matrix of 1;s and 0;s using volts matrix.  
    rvolts<-volts
    for(p in 1:21) {
      for (q in 1:21) {
        if(rvolts[p,q]==max(volts[p,])) 
        {rvolts[p,q]<-1} else {rvolts[p,q]<-0}
      }}
    
    #Useful output-
    #rvoltsls-a list of rvolts matrices.          #voltage- a dataframe that's used in processing the output. 
    rowvolts<-c()   
    
    
    for (i in 1:21) {
      rowvolts[length(rowvolts)+1]<-sum(rvolts[i,])
    }
    sumrv<-sum(rowvolts)
    
    if(sum(rowvolts)==21)
    {voltageout<-rbind(voltageout,c(T,dp,sumrv))
    rvoltsls[[length(rvoltsls)+1]]<-rvolts
    }
    if(sum(rowvolts)<21)
    {voltageout<-rbind(voltageout,c(T,dp,sumrv))
    rvoltsls[[length(rvoltsls)+1]]<-rvolts
    
    }
  } 
}

#END PROGRAM



















plot(energytemps)
####Section IV
###FINDING ALL VALID PATHS AND GETTING LISTS OF MATRICES OF VALID PATHS

####Section IV-A
##We first need to get column and rowsums.
#rvoltsround- a list of the rounded matrices of rvoltsls.This can be stream lines
rvoltsround<-lapply(rvoltsls,round)
length(rvoltsround)


#Column sums and Row sums of outputs.
sumsrvoltsround<-lapply(rvoltsround,colSums)
sumsrvoltsround
rowsumsrvoltsround<-lapply(rvoltsround,rowSums)
rowsumsrvoltsround

###Section IV-B
##Finding valid paths.
#roundvector- tells you which elements of rvoltsround represent valid paths.
#goodvector-a list of matrices representing good paths.

voltageout<-voltageout[-1,]
voltageout
roundvector<-c()          
#roundvectors find valid paths
for (i in 1:length(sumsrvoltsround)) {
  transformsum<-c()
  for (j in 1:21) {
    transformsum[j]<-sumsrvoltsround[[i]][j]/sumsrvoltsround[[i]][j]
  }
  for (k in 1:21) {
    if (is.nan(transformsum[k])) {
      transformsum[k]<-0
    }
  }
  print(transformsum)
  if (sum(rowsumsrvoltsround[[i]])==21 & sum(transformsum)==21) {
    roundvector[length(roundvector)+1]<- i
  }
}
roundvector
length(roundvector)

goodvolts<-list()
for (i in 1:length(roundvector)) {
  goodvolts[[length(goodvolts)+1]]<- rvoltsls[[roundvector[i]]]
}

roundvolts<-lapply(goodvolts,round)
length(roundvolts)


####Section V
###Output Paths in the form (AminoAcid Names, Paths) 
####NOTE YOU MUST RUN FUNCTION IN SECTION VII BEFORE YOU CAN GO THROUGH THIS SECTION.

AminoMaps<-lapply(roundvolts,maps)
AminoMaps

length(AminoMaps)


aamapsreOrder<-lapply(AminoMaps,reorderAA)
aamapsreOrder



###SECTION VI
##IGNORE, GO TO SECTION VII

#Write 
naa<-c() 

for (i in 1:length(aaMapsFinal)){
  naa[i]<-paste("titvx3",i,sep="_")
  naa[i]<-paste(naa[i],"txt",sep=".")
  naa[i]<-paste("C:/Users/Brian Sulkow/Documents/AANeighbors",naa[i],sep="/")}
naa


for (i in 1:length(aamaps)) {
  write.table(aaMapsFinal[[i]], naa[i], sep="\t")
}


codonnames

results<- read.csv("C:/Users/Brian Sulkow/Documents/AANeighbors/testpol_5.csv",header = FALSE)
summary(results)
results


####SECTION VII
#FUNCTIONS FOR PROCESSING OUTPUT INTO 21 PAIRS (AMINO ACID NAMES, TIME)
#RUN this section before running section V


#maps-Outputs 21 pairs, (Amino Acid Number,Time) .  Notice we get the amino acid number.  The next path with get the pairs.  

maps<-function(data) {
  #the matrix is 21x21  we want to 
  
  
  mapper<-data.frame()
  
  for (l in 1:21){
    for (m in 1:21){ 
      if (data[m,l]!=0 )
      {mapper<-rbind(mapper,c(l,m))}
    }
  }
  colnames(mapper)<-c("AA","codonnum")
  return(mapper)
  
}



#make a AA names vector
 
#change the factor levels of aadistances. THis is a necessary step to avoid confusion. 
newAALevels<-factor(aadistances[,1], levels=(aadistances[,1]))
newAALevels

###reorderAA- Takes the output of maps (aminoacid number, time) and outputs (aminoacid NAME, time)  
reorderAA<-function(data) {
  AAPerm<-c()
  for (i  in 1:21){
  AAPerm[i]<- as.character(newAALevels[data[i,2]])
    
  }
  data[,3]<-AAPerm
  return(data)
}


####SECTION VIII   
###STOP. MOST WHAT"S BELOW IS NOT NEEDED. I WILL TAKE CONTENTS OUT SHORTLY


#aadistances is a frame of aa names and measurements according to pol, hydor, and vol.  sortaadistances reorders the rows with respect to the polarity measurements. 

aadistances<- read.csv("C:/Users/Brian Sulkow/Documents/MonteCarloClub/Neural Nets/TSP/aadistances.csv",header = TRUE)


aadistances
?lapply
sortaadistances<-aadistances[order(aadistances$polarity),]
sortaadistances
###Making an output fit for Weigangs code-stats.pl



#insertAminoAcidCol<- function will be used to add a column to the aamaps data with names of names the amino acid with the data. 
insertAminoAcidCol<-function(data) {
  AminoAcids<-c()
  for (i in 1:64) {
    AminoAcids[i]<-as.character(sortaadistances[data[i,1],1])
  } 
  data[,5]<-AminoAcids
  return(data)
}
aadistances


##capitalLetters function-  makes all letters capital

capitalLetters<- function(data) {
  data.frame(lapply(data, function(v) {
    if (is.character(v)) return(toupper(v))
    else return(v)
  }))
}


#cutCols function to cut two columns  Codon and AA  
# change from cols 4,5 to 3,5 for oldworks
cutCols<- function(data) {
  data<-data.frame(data[,4],data[,5])
  return(data)
  
}

aaMapsFinal

#check that all there are no repeat
checkUnique<-c()
for (i in 1:length(aaMapsFinal)) {
  checkUnique[i]<-length(unique(aaMapsFinal[[1]][,1]))
}
checkUnique

#OldWork- Import the files
oldwork<-list()
for (i in 1:5) {
  oldwork[[i]]<-paste("maps",i, sep="_")
  oldwork[[i]]<-paste(oldwork[[i]],"txt",sep=".")
  oldwork[[i]]<-paste("C:/Users/Brian Sulkow/Documents/MonteCarloClub/rcode",oldwork[[i]],sep="/")
  
}

oldWorkReads<-lapply(oldwork, read.table, header=T)
oldWorkReads

##Computing Distances
computeDistance<- function (data) {
  totaldistance<-0
  for (i in 1:20) {
    totaldistance<- totaldistance+ NormaaMatrix[data[i,2],data[i+1,2]]
  }
  return(totaldistance)
}


aminoDistances<-lapply(aamapsreOrder,computeDistance)
aminoVectorDist<-c()
for ( i in 1:length(aminoDistances)) {
  aminoVectorDist[i]<-aminoDistances[[i]]
}
boxplot(aminoVectorDist)
min(aminoVectorDist)
aminoVectorDist
subset(aminoVectorDist, aminoVectorDist>5)
newAALevels
length(aamapsreOrder)
