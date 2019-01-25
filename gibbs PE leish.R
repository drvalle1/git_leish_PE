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

#get dat complete
dat.complete=dat
for (i in 2:nanos){
  cond=dat.complete[,i-1]==1
  dat.complete[cond,i]=1
}

#get distance
setwd('U:\\anaia\\derived data')
dist=read.csv('matriz distancia.csv',as.is=T)
rownames(dist)=dist$X
ind=which(colnames(dist)=='X')
dist=dist[,-ind]
OneOverDist=1/data.matrix(dist)
diag(OneOverDist)=0
soma=rowSums(OneOverDist)

#get invasion pressure
IP=matrix(NA,nloc,nanos)
for (i in 2:nanos){
  tmp=OneOverDist*matrix(dat.complete[,i-1],nloc,nloc,byrow=T)
  IP[,i]=apply(tmp,1,sum,na.rm=T)/soma  
}

#get covariates
setwd('U:\\anaia\\derived data\\covs')
z=list.files()
ind=grep('a1_',z); z1=z[ind]; z1
ncov=3#length(z1)
wmat=list()
for (i in 1:ncov){
  tmp=read.csv(z1[i],as.is=T)
  wmat[[i]]=data.matrix(tmp)
}
n=length(wmat)
wmat[[n+1]]=IP
nparam=n+1

#useful stuff
Identifiers=matrix(1:(nloc*nanos),nloc,nanos)

#get covariates for detection probability
xmat=matrix(0,nloc,nanos)
colnames(xmat)=paste0('x',7:17)
xmat[,paste0('x',15:17)]=1
xmat=cbind(1,as.numeric(xmat)) #add 1's for intercept

#initial values for parameters
alpha=rep(0,length(wmat)+1) #for intercept
betas=rep(0,2)
sd.alpha=c(sqrt(10),rep(1,length(alpha)-1))
z1=dat.complete

medias=get.medias(wmat=wmat,alpha=alpha,nloc=nloc,nanos=nanos,nparam=nparam)
ind.random=list()
for (i in 1:nloc){
  tmp=which(dat[i,]==1)
  tmp1=0:nanos
  if (dat[i,1]==1) tmp1=numeric()
  if (length(tmp) != 0 & dat[i,1]!=1) tmp1=0:(min(tmp)-1)
  ind.random[[i]]=tmp1
}

#stuff for gibbs sampler
ngibbs=1000
nburn=ngibbs/2
vec.z=matrix(NA,ngibbs,nloc*nanos)
vec.betas=matrix(NA,ngibbs,ncol(xmat))
vec.alpha=matrix(NA,ngibbs,length(wmat)+1) #+1 because of intercept
jump1=list(betas=rep(1,ncol(xmat)),alpha=rep(1,length(alpha)))
accept1=list(betas=rep(0,ncol(xmat)),alpha=rep(0,length(alpha)))
accept.output=50
for (i in 1:ngibbs){
  print(i)

  z1=sample.z(nloc=nloc,nanos=nanos,ind.random=ind.random,z1=z1,
              medias=medias,xmat=xmat,betas=betas,dat=dat)
  # z1=z1.true

  tmp=sample.alpha(z1=z1,IP=IP,nanos=nanos,sd.alpha=sd.alpha,alpha=alpha,jump=jump1$alpha,
                   wmat=wmat,nloc=nloc,nparam=nparam)
  alpha=tmp$alpha
  accept1$alpha=accept1$alpha+tmp$accept
  # alpha=alpha.true
  medias=get.medias(wmat=wmat,alpha=alpha,nloc=nloc,nanos=nanos,nparam=nparam)
  
  tmp=sample.betas(betas=betas,jump=jump1$betas,xmat=xmat,dat=dat,z1=z1)
  accept1$betas=accept1$betas+tmp$accept
  betas=tmp$betas
  # betas=betas.true
  
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