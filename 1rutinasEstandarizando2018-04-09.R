#Variable Selection for depth measures

#rm(list=ls(all=TRUE))

install.packages("MASS")
install.packages("stats")
install.packages("graphics")
install.packages("combinat")
install.packages("gdata")    
install.packages("ggplot2")
install.packages("GA")
install.packages("mvtnorm")
install.packages("depth")
install.packages("parallel")
install.packages("doParallel")
install.packages("foreach")
install.packages("iterators")
install.packages("ccaPP")
install.packages("DepthProc")

library(MASS)
library(stats)
library(graphics)
library(combinat)
library(gdata)    
library(ggplot2)
library(GA)

library(depth)
library(parallel)
library(doParallel)
library(foreach)
library(iterators)
library(ccaPP)
library(DepthProc)



######################
# Squared Correlations:
######################
# The function "correlation2" calculates the squared correlation between
# two vectors U and V for different association measures.

correlation2<-function(U,V,correl){
  if(correl=="spearman"){
    cca.res=(corSpearman(U,V,consistent=TRUE))^2
  }
  if(correl=="kendall"){
    cca.res=(corKendall(U,V,consistent=TRUE))^2
  }
  if(correl=="quadrant"){
    cca.res=(corQuadrant(U,V,consistent=TRUE))^2
  }
  if(correl=="M"){
    cca.res=(corM(U,V))^2
  }
  if(correl=="pearson"){
    cca.res=(corPearson(U,V))^2
  }
return(cca.res)
}


###############
# depths:
###############
# Inputs:
# U --> Numerical vector or matrix whose depth is to be calculated. Dimension has to be the same as that of the observations.
# V --> The data as a matrix. If it is a matrix or data frame, then each row is viewed as one multivariate observation.
# prof --> depth measure to be used: "FM" =Fraiman-Muniz, "tukey"= tuckey hyperplane, "proy"=Projection Depth, "oja"- 
depths=function(U,V,prof){
 
  if(dim(U)[2]==1){
    result=rank(U[,1])/length(U[,1])
  }
  else{
    if(prof=="FM"){
      result=fncDepthFM(U,V)
    }
    if(prof=="tukey"){
      result=depthTukey(U,V)
    }
    if(prof=="proy"){
      set.seed(999)
      result=depthProjection(U,V, ndir = 1000, threads = -1)
    }
    if(prof=="oja"){
      result=depth(U,V,method="Oja")
    }
    
  } 
  result
}

#####################
# selvardepthestanga:
#####################
# The function "selvardepthestanga" selects variables using genetic algorithms.
# Is an auxiliary objective function to be maximized inside the genetic algorithm ("GA") function.
# Inputs:
# chromosome --> an initial 0-1 vector, where coordinates with 1 indicate that the variable is taken into account.
# datos --> an nxp data matrix. Each row is an observation.
# crit --> takes value 1 when we minimize MSE and 2 when we maximize correlation.
# lambda --> a penalty parameter on the quantity of variables.
# maxVar --> the maximum quantity of variables admited in the minimization / maximization.
# origDepthYes --> set it 1 when there are previously calculated depths for the analysed data and allows to input the values of such depths and set it as 0 in other case.
# origDepths --> a vector with the n depths of the n observations previously known or set it as 1:n and set origDepthYes=0 otherwise.
# depthFncOrig --> set the depth function to be used to measure the original depths if origDepthYes=0.
# depthFncSearch --> set the depth function to be used to measure depths when minimizing MSE or maximizing correlation.
# indices --> set the indices to be used to calculate the original depths. Set indices=1:p in general, but can fix some particular set of some of the variables to test wether the method is finding the right variables.
# correl --> set the correlation measure to be used when maximizing correlation. Possibilities: "sperman", "kendall", "quadrant", "M" or "pearson".
# outputs:
# fobj2 --> the value of the objective function to be maximized

selvardepthestanga =function(chromosome,datos,crit,lambda,maxVar,origDepthYes,origDepths,depthFncOrig,depthFncSearch,indices,correl){
                                       
  # We set the output fobj2=-10 when for chromosomes with more variables than the maximium admited in order to discard such chromosomes as possible maximizers.
  if(sum(chromosome)>maxVar){
    fobj2=-10
  }
  else{
    # We call "n" the number of observations.
    n=dim(datos)[1]
    # We call "Xdepth" to the depths of the observations that we wish to approximate with less variables. 
    if(origDepthYes==1){
      Xdepth=origDepths
    }
    else{
      # if indices=1:p and origDepthYes=0 then Xdepth are the depths of the observations with all variables.
      Xdepth=depths(datos[,indices],datos[,indices],depthFncOrig)  
    }
    fobj=10 #"fobj2" is the objective function we maximize, and "fobj2" will be "-fobj". We set it like this in order to discard not choosing any variable as a possible maximizer. 
    if (sum(chromosome)>0){
      # "Xnew" is the data matrix restricted to the variables present in the chromosome.
      Xnew=as.matrix(datos[,which(chromosome==1)])
      # "Xnewdepth" gets the depths of the observations if we only consider the variables present in the chromosome.
      Xnewdepth=depths(Xnew,Xnew,depthFncSearch)
    
      if (crit==1){
        #"ecm" is the MSE between the standardized depths of the observations calculated with all variables and with the variable of the chromosome.
        ecm=mean((rank(Xdepth)/n-rank(Xnewdepth)/n)^2) 
        fobj=ecm+lambda*sum(chromosome) # we add the penalty.
        
      }
      if (crit==2){
        # "corr" is the squared correlation between the depths calculated with all variables and with the variables of the chromosome.
        corr=correlation2(Xdepth,Xnewdepth,correl) 
        fobj=(1-corr)+lambda*sum(chromosome) # fobj is a function one wish to minimize under both criteria.
      } 
      
    }
    fobj2=-fobj 
    # we multiply by -1 because the fitness function that take the Genetic Algorithm "GA" is a function to be maximized.
    
  }
  return(fobj2)  
}




########################
#selvardepthestanexahus:
########################
# The function "selvardepthestanexahus" excecutes an exhaustive search of variables maximizing/minimizing.
# INPUTS:
# X --> an nxp data matrix. Each row is an observation.
# origDepthYes --> set it 1 when there are previously calculated depths for the analysed data and allows to input the values of such depths and set it as 0 in other case.
# origDepths --> a vector with the n depths of the n observations previously known or set it as 1:n and set origDepthYes=0 otherwise.
# indices --> set the indices to be used to calculate the original depths. Set indices=1:p in general, but can fix some particular set of some of the variables to test wether the method is finding the right variables.
# lambda --> a penalty parameter on the quantity of variables.
# n1 --> maximum number of variables to be selected.
# solini --> vector with 0's and 1's, 1's indicates the presence of the variable for the analysis.
# crit --> takes value 1 when we minimize MSE and 2 when we maximize correlation.
# depthFncOrig --> set the depth function to be used to measure the original depths if origDepthYes=0.
# depthFncSearch --> set the depth function to be used to measure depths when minimizing MSE or maximizing correlation.
# correl --> set the correlation measure to be used when maximizing correlation. Possibilities: "sperman", "kendall", "quadrant", "M" or "pearson".
# OUTPUTS:
# exhaustivo --> matrix where each row represents the result for a set of variables. In the last column we observe the objective function and we observe ones in the columns of the variables which belong to the set of variables and zeros in the others.
# bestfobj --> matrix where each row shows for each number of variables (1,2,3,...) the optimum value of the objective function.
# bc --> matrix where each the i-th row shows the variables of the optimum set of i variables.

selvardepthestanexahus=function(X,origDepthYes,origDepths,indices,lambda,n1,solini,crit, depthFncOrig,depthFncSearch,correl){  
  
  n=dim(X)[1]
  # "np" is the total number of variables where we make the exhaustive search.
  np=sum(solini==1)
  # "combinaciones" will be a matrix where each row represents a possible set of variables.
  combinaciones=t(rep(0,dim(X)[2])) #the first row wont be taken into account.
  
  fobj=0 # we set the objective function equal to 0.
  
  # We call "Xdepth" to the depths of the observations that we wish to approximate with less variables. 
  if(origDepthYes==1){
    Xdepth=origDepths
  }
  if(origDepthYes==0){
    Xdepth=depths(X[,indices],X[,indices],depthFncOrig)  
  }
  
  
  # the sets of variables for which we will maximize/minimize will have no mor than min(np,n1) variables.
  for (i in 1:min(np,n1)){
   fcomb=combn(which(solini==1),i) #each column is a possible set of i variables
   # the matrix "aux" is a 1-0 matrix where each row represents a possible set of i variables. The "1" indicates the presence of the variable that corresponds to the column.
   aux=matrix(0,nrow=dim(as.matrix(fcomb))[2],ncol=dim(X)[2])
    # we insert a "1" in the desired places.
    for (jj in 1:nrow(aux)) {
      aux[jj,as.matrix(fcomb)[,jj]] = 1
    }
    combinaciones=rbind(combinaciones,aux)
    szfcomb=dim(as.matrix(fcomb))[2] #quantity of sets of i variables of the possible variables of search (solini==1).
    
    for (j in 1:szfcomb){
      # By varying "j", "xnew" varies between the data matrix restricted to the possible sets of i variables.
      xnew=as.matrix(X[,which(aux[j,]==1)])
      # "Xnewdepth" gets the depths of the observations of "xnew".
      Xnewdepth=depths(xnew,xnew,depthFncSearch)  
      
      if (crit==1){
        #"ecm" is the MSE between the standardized depths of the observations calculated with all variables and with the variable of xnew.
        ecm=mean((rank(Xdepth)/n-rank(Xnewdepth)/n)^2) 
        # in "fobj" we save the values of the penalized objective function.
        fobj=cbind(fobj,ecm+lambda*sum(aux[j,]))
      }
      else{
        corr2=correlation2(Xdepth,Xnewdepth,correl) # "corr2" is the squared correlation between the depths calculated with all the variables and with the variables of xnew
        fobj=cbind(fobj,(1-corr2)+lambda*sum(aux[j,]))} 
    }
   
}
  # "exhaustivo" will save all possible sets of variables and its respective value for the penalized objective function.
  exhaustivo=cbind(combinaciones,t(fobj))
  
  
  nrovar=1:min(n1,np)  #all the possible number of variables.
  # in "bestfobj" we save the optimum values of the penalized objective function for the differente quantity of variables.
  bestfobj=rep(0,max(nrovar))
  # in "bestcomb" we save the optimum variables of the penalized objective function for the differente quantity of variables.
  bestcomb=matrix(nrow=max(nrovar),ncol=dim(combinaciones)[2])
  
  for (i in nrovar){
    bestfobj[i]=min(fobj[which(apply(combinaciones,1,sum)==i)])
    if(bestfobj[i] != 0){
      #we only save the first set of variables which optimizes if "fobj" is different from 0.
      bestcomb[i,]=combinaciones[which(fobj==bestfobj[i])[1],]       
    }
    if(bestfobj[i] == 0){
      #we only save the first set of variables which optimizes if "fobj" is different from 0 discarding the first row which is equal to 0.
      bestcomb[i,]=combinaciones[which(fobj==bestfobj[i])[2],]      
    }
    
  
  }
  #bc will be a list of the selected variables (best combinations) for each possible quantity of variables.
  bc=list()
  for (i in nrovar){ bc[[i]]=which(bestcomb[i,]>0)}
  #bestGlobalFobj is the optimum value of the objective function independent of the quantity of variables.
  bestGlobalFobj=min(bestfobj)
  #bestGlobalComb are the variable where the optimum is obtained.
  bestGlobalComb=bc[[which(bestfobj==bestGlobalFobj)[1]]]
  return(list(exhaustivo=exhaustivo,bestfobj=bestfobj,bc=bc,bestglobalfobj=bestGlobalFobj,bestglobalcomb=bestGlobalComb))
#  return(list(exhaustivo=exhaustivo,bestfobj=bestfobj,bc=bc,bestglobalfobj=bestGlobalFobj)).
  
}


#####################
# cross validation
#####################
# selvardepthestanCV:
#####################
# Inputs:
# X --> an nxp data matrix. Each row is an observation.
# origDepthYes --> set it 1 when there are previously calculated depths for the analysed data and allows to input the values of such depths and set it as 0 in other case.
# origDepths --> a vector with the n depths of the n observations previously known or set it as 1:n and set origDepthYes=0 otherwise.
# depthFncOrig --> set the depth function to be used to measure the original depths if origDepthYes=0.
# depthFncSearch --> set the depth function to be used to measure depths when minimizing MSE or maximizing correlation.
# maxIter --> the maximum quantity of iterations when using Genetic algorithm.
# indices --> set the indices to be used to calculate the original depths. Set indices=1:p in general, but can fix some particular set of some of the variables to test wether the method is finding the right variables.
# method --> the method of search. The options are: "GA" (genetic Algorithm) or "EX" (the exhaustive search).
# Lambdas --> the set of possible values of lambda to use as penalty parameters in the cross validation.
# kfold --> the total of groups in which we divide to run the cross validation.
# n1 --> maximum number of variables to be selected when using the exhaustive method.
# solini --> vector with 0's and 1's, 1's indicates the presence of the variable for the analysis when using the exhaustive method.
# crit --> takes value 1 when we minimize MSE and 2 when we maximize correlation.
# maxVar --> the maximum quantity of variables admited in the minimization / maximization.
# correl --> set the correlation measure to be used when maximizing correlation. Possibilities: "sperman", "kendall", "quadrant", "M" or "pearson".
##########################################################
selvardepthestanCV=function(X,origDepthYes,origDepths,depthFncOrig,depthFncSearch,maxIter,indices,method,Lambdas,kfold,n1,solini, crit,maxVar,correl){
                                                                  
  nn=dim(X)[1]   # the quantity of observations.
  cantxgrupo=nn%/%kfold  # the quantity of observations in each group.
  # since "nn" might not be divisible by "kfold" we need to aggregate 1 observation to some of the groups.
  resto=nn%%kfold
  sobrevive=1:nn # from this vector we will sample the indices of the observations of the groups.
  grupo<-list()  #will be the list of the sampled groups.
  
  for(gg in 1:kfold){
    if(gg<=resto){
      grupo[[gg]]=sample(sobrevive,size=cantxgrupo+1)
      sobrevive=sobrevive[-grupo[[gg]]]
    }
    else{
      grupo[[gg]]=sample(sobrevive,size=cantxgrupo)
    }
  }
  #We will calculate the objective function for each group and each lambda .
  gruposCant=length(grupo)  # the quantity of groups.
  lambdasCant=length(Lambdas) # the quantity of lambdas.
  # for each group and each lambda, we calculate the optimum indices  excuding the observations of the group and penalizing with lambda.
  # for such indices we calculate the objective function using the obtained indices in the observations wich were excluded previously.
  error2Prom=vector(,length=lambdasCant) # the average of the objective function between the groups for each lambda minimizing MSE.
  corr2Prom=vector(,length=lambdasCant) # the average of the objective function between the groups for each lambda maximizing correlation.
  for(ll in 1:lambdasCant){
    lambda=Lambdas[ll]
    error2=vector(,length=gruposCant)
    corr2=vector(,length=gruposCant)
    for(gg in 1:gruposCant){
      todosMenos=list()
      indicesOpt=list()
      #"todosmenos[[gg]]" will be the indices of the observations wich do not belong to group gg.
      todosMenos[[gg]]=(1:nn)[-grupo[[gg]]]
      if(method=="GA"){
        
        varsel=ga(type="binary",fitness=selvardepthestanga,datos=X[todosMenos[[gg]],],crit,lambda,maxVar,origDepthYes,origDepths,depthFncOrig,depthFncSearch,indices,correl,nBits=ncol(X[todosMenos[[gg]],]),maxiter=maxIter)
                                                                                                                                    
        indicesOpt[[gg]]=as.vector(which(varsel@solution[1,]==1))
      }
      if(method=="EX"){
        varsel=selvardepthestanexahus(X[todosMenos[[gg]],],origDepthYes,origDepths[todosMenos[[gg]]],indices,lambda,n1,solini,crit,depthFncOrig,depthFncSearch,correl)
                                                                                                                                  
        indicesOpt[[gg]]=as.vector(varsel$bestglobalcomb)
        
        
      }
      
      indexN=length(indicesOpt[[gg]])
      #We calculate the depths of the observations of group gg both with all the variables and witth the optimum variables obtained without using gg.
      XdepthEst=depths(as.matrix(X[grupo[[gg]],indicesOpt[[gg]]]),as.matrix(X[grupo[[gg]],indicesOpt[[gg]]]),depthFncSearch) 
      XdepthOrig=depths(as.matrix(X[grupo[[gg]],]),as.matrix(X[grupo[[gg]],]),depthFncSearch)
      if(crit==1){
        if(gg<=resto){
          error2[gg]=sum((rank(XdepthEst)/(cantxgrupo+1)-rank(XdepthOrig)/(cantxgrupo+1))^2) 
        }
        else{
          error2[gg]=sum((rank(XdepthEst)/(cantxgrupo)-rank(XdepthOrig)/(cantxgrupo))^2)
        }
        
      }
      if(crit==2){
        
        corr2[gg]=correlation2(XdepthEst,XdepthOrig,correl)
        
      }
    }
    if(crit==1){error2Prom[ll]=mean(error2)}
    if(crit==2){corr2Prom[ll]=mean(corr2)}  
  }
  if(crit==1){lambdaOpt=Lambdas[which(error2Prom==min(error2Prom))[1]]}  
  if(crit==2){lambdaOpt=Lambdas[which(corr2Prom==max(corr2Prom))[1]]}
  if(method=="GA"){
    varselFin=ga(type="binary",fitness=selvardepthestanga,datos=X,crit,lambda=lambdaOpt,maxVar,origDepthYes,origDepths,depthFncOrig,depthFncSearch,indices,correl,nBits=ncol(X),maxiter=maxIter)
    fobj=varselFin@fitnessValue
    combinaciones=as.vector(varselFin@solution[1,])
    lambdaopt=lambdaOpt
    
  }  
  if(method=="EX"){
    varselFin=selvardepthestanexahus(X,origDepthYes,origDepths,indices,lambdaOpt,n1,solini,crit,depthFncOrig,depthFncSearch,correl)
                                    
    
    fobj=varselFin$bestglobalfobj
    combinaciones=varselFin$bestglobalcomb
    lambdaopt=lambdaOpt
  }
    list(fobj=fobj,combinaciones=combinaciones,lambdaopt=lambdaopt)
}
  
