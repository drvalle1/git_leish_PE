rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(1)

#read relevant functions
setwd('U:\\GIT_models\\git_leish_PE')
source('gibbs PE leish functions.R')
sourceCpp('aux_PE_leish.cpp')

#get incidence
setwd('U:\\GIT_models\\git_leish_PE\\simulated data')
dat=data.matrix(read.csv('simulated data.csv',as.is=T))
nloc=nrow(dat)
nanos=ncol(dat)

#get distance
setwd('U:\\anaia\\derived data')
dist=read.csv('matriz distancia.csv',as.is=T)
rownames(dist)=dist$X
ind=which(colnames(dist)=='X')
dist=dist[,-ind]
OneOverDist=1/data.matrix(dist)
diag(OneOverDist)=0
SomaOneOverDist=rowSums(OneOverDist)

#get covariates
setwd('U:\\anaia\\derived data\\covs')
z=list.files()
ind=grep('a1_',z); z1=z[ind]; z1
ncov=2#length(z1)
wmat=list()
for (i in 1:ncov){
  tmp=read.csv(z1[i],as.is=T)
  wmat[[i]]=data.matrix(tmp)
}
nparam=length(wmat)

#get covariates for detection probability
xmat=matrix(0,nloc,nanos)
colnames(xmat)=paste0('x',7:17)
xmat[,paste0('x',15:17)]=1
xmat=cbind(1,as.numeric(xmat)) #add 1's for intercept

#initial values for parameters
alpha=rep(0,length(wmat)+1+1) #for intercept and for invasion pressure
betas=rep(0,2)
sd.alpha=c(sqrt(10),rep(1,length(alpha)-1))

#get latent z's
z1=matrix(NA,nloc,nanos)
for (i in 1:nloc){
  tmp=dat[i,]
  ind=which(tmp==1)
  tmp1=rep(0,nanos); 
  if (length(ind)> 0) {
    ind1=min(ind)
    tmp1[ind1:nanos]=1
  }
  z1[i,]=tmp1
}

#useful stuff
Identifiers=matrix(1:(nloc*nanos),nloc,nanos)
IP=CalcInvasionPressure(z=z1, OneOverDist=OneOverDist, 
                        nanos=nanos, nlocs=nloc,
                        SomaOneOverDist=SomaOneOverDist)
IP[,nanos]=NA

medias=get.medias(wmat=wmat,alpha=alpha,nloc=nloc,nanos=nanos,IP=IP)
ind.random=list()
for (i in 1:nloc){
  tmp=which(dat[i,]==1)
  tmp1=1:(nanos+1) #to account for the possibility of remaining uninvaded
  if (dat[i,1]==1) tmp1=numeric()
  if (length(tmp) != 0 & dat[i,1]!=1) tmp1=1:(min(tmp)-1)
  ind.random[[i]]=tmp1
}

#stuff for gibbs sampler
ngibbs=1000
nburn=ngibbs/2
vec.z=matrix(NA,ngibbs,nloc*nanos)
vec.betas=matrix(NA,ngibbs,ncol(xmat))
vec.alpha=matrix(NA,ngibbs,length(wmat)+1+1) #intercept + invasion pressure
jump1=list(betas=rep(1,ncol(xmat)),alpha=rep(1,length(alpha)))
accept1=list(betas=rep(0,ncol(xmat)),alpha=rep(0,length(alpha)))
accept.output=50
for (i in 1:ngibbs){
  print(i)

  z1=sample.z(nloc=nloc,nanos=nanos,ind.random=ind.random,z1=z1,
              OneOverDist=OneOverDist,SomaOneOverDist=SomaOneOverDist,
              wmat=wmat,alpha=alpha,xmat=xmat,betas=betas,dat=dat)
  # z1=z1.true
  IP=CalcInvasionPressure(z=z1, OneOverDist=OneOverDist, 
                          nanos=nanos, nlocs=nloc,
                          SomaOneOverDist=SomaOneOverDist)

  medias=get.medias(wmat=wmat,alpha=alpha,nloc=nloc,nanos=nanos,IP=IP)

  tmp=sample.betas(betas=betas,jump=jump1$betas,xmat=xmat,dat=dat,z1=z1)
  accept1$betas=accept1$betas+tmp$accept
  betas=tmp$betas
  # betas=betas.true
  
  tmp=sample.alpha(z1=z1,IP=IP,nanos=nanos,sd.alpha=sd.alpha,alpha=alpha,jump=jump1$alpha,
                   wmat=wmat,nloc=nloc)
  alpha=tmp$alpha
  accept1$alpha=accept1$alpha+tmp$accept
  # alpha=alpha.true
  medias=get.medias(wmat=wmat,alpha=alpha,nloc=nloc,nanos=nanos,IP=IP)
  
  #adaptation MH algorithm
  if (i<nburn & i%%accept.output==0){
    k=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output)
    accept1=k$accept1
    jump1=k$jump1
  }

  #store results
  vec.z[i,]=z
  vec.betas[i,]=betas
  vec.alpha[i,]=alpha
}
z1.estim=z1
betas.estim=betas
alpha.estim=alpha