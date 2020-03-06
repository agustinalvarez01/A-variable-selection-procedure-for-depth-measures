rm(list=ls(all=TRUE))


library("readxl")
library("mvtnorm")
source("1programa_simul_2018-04-09.R")

################
# simulacionaux:
################
# The function "simulacionaux" is used to perform the simulation.
# Inputs:
# nrep --> the quantity of replications of the simulation study.
# nmuestra --> the size of the sample in each replication.
# nx --> the dimension of the X vector.
# nz --> the dimension of the Z vector.
# nw --> the dimension of the W vector.
# sigy --> the standard deviation of the Y variables.
# sigeps --> the standard deviation of the Epsilon vector. 
# sigw --> the standard deviation of the W variables. 
# rho --> the correlation between distinct X variables.
# the other input variables correspond to the auxiliary functions used in the simulation.

simulacionaux<-function(origDepthYes,origDepths,depthFncOrig,depthFncSearch,indices,nrep,nmuestra,nx,nz,nw,sigy,sigeps,sigw,rho,method,n1,solini, crit,maxIter,correl){
  
  indi=toString(indices)
  # we define the name of the files to be saved after simulation
  fobj.global.arch<-paste("fobj_global_metodo_",method,"_ODY_",origDepthYes,"_DFO_",depthFncOrig,"_DFS_",depthFncSearch,"_indices_",indi,"_nrep_",nrep,"_nmuestra_",nmuestra,"_nx_",nx,"_nz_",nz,"_nw_",nw,"_sigy_",sigy,"_sigeps_",sigeps,"_sigw_",sigw,"_rho_",rho,"_n1_",n1,"_crit_",crit,"_maxiter_",maxIter,"_cor_",correl,".txt",sep="")    
  combi.arch<-paste("combi_metodo_",method,"_ODY_",origDepthYes,"_DFO_",depthFncOrig,"_DFS_",depthFncSearch,"_indices_",indi,"_nrep_",nrep,"_nmuestra_",nmuestra,"_nx_",nx,"_nz_",nz,"_nw_",nw,"_sigy_",sigy,"_sigeps_",sigeps,"_sigw_",sigw,"_rho_",rho,"_n1_",n1,"_crit_",crit,"_maxiter_",maxIter,"_cor_",correl,".txt",sep="")
  lambda.opt.arch<-paste("lambda_opt_metodo_",method,"_ODY_",origDepthYes,"_DFO_",depthFncOrig,"_DFS_",depthFncSearch,"_indices_",indi,"_nrep_",nrep,"_nmuestra_",nmuestra,"_nx_",nx,"_nz_",nz,"_nw_",nw,"_sigy_",sigy,"_sigeps_",sigeps,"_sigw_",sigw,"_rho_",rho,"_n1_",n1,"_crit_",crit,"_maxiter_",maxIter,"_cor_",correl,".txt",sep="")
  Sigma=matrix(rep(rho,nx*nx),nx,nx)+diag(1-rho,nx)
  # we start the simulation
  for(m in 1:nrep){
    print(m)
    print(Sys.time())
    set.seed(999+m)
    # we define the vectors X, Z and W.
    x=rmvnorm(nmuestra,mean=rep(0,nx),sigma=Sigma)
    y=sigy*rnorm(nmuestra)
    eps=rmvnorm(nmuestra,mean=rep(0,nz),sigma=diag(rep(sigeps,nz)))
    z=matrix(,nmuestra,nz)
    for(nn in 1:nz){
      z[,nn]=(x[,1]+nn*y+eps[nn])/2       
    }
    w=rmvnorm(nmuestra, sigma=diag(rep(sigw,nw)))
    datos=cbind(x,z,w)
    origDepths=depths(datos[,indices],datos[,indices],depthFncOrig)
    # we perform the variable selection procedure.
    salida=selvardepthestanexahus(datos,origDepthYes,origDepths,indices,0.01,n1,solini,crit,depthFncOrig,depthFncSearch,correl)
    
    # we save the outputs and write the results to the files .
    guardo.fobj=c(m,salida$bestfobj)
    guardo.comb=list()
    write(t(guardo.fobj),file=fobj.global.arch,ncolumns=n1+1,append=T)
    for(i in 1:n1){
      guardo.comb[[i]]=c(m,salida$bc[[i]])  #ya no lo entiendo
      write(t(guardo.comb[[i]]),file=combi.arch,ncolumns=dim(datos)[2]+1,append=T)
    }
    
  }
  
}


#### we show an example

nrep=5
nmuestra=50
nx=6
nz=6
nw=2
sigy=2
sigeps=0.01
sigw=0.1
rho=0.90

method="EX"
n1=nx+nz+nw # 
solini=rep(1,n1) # We use all the variables in the analysis
crit= 1 # we minimize MSE
depthFncSearch ="FM"
depthFncOrig="FM"
origDepthYes=0
origDepths=1:nmuestra
indices=1:n1
maxVar=4
maxIter=3
kfold=3
correl="pearson"

simulacionaux(origDepthYes,origDepths,depthFncOrig,depthFncSearch,indices,nrep,nmuestra,nx,nz,nw,sigy,sigeps,sigw,rho,method,n1,solini, crit,maxIter,correl)
