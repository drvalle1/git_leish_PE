rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(4)

#read relevant functions
setwd('U:\\GIT_models\\git_leish_PE')
source('gibbs PE leish functions.R')
source('gibbs PE leish.R')
sourceCpp('aux_PE_leish.cpp')

#get incidence
setwd('U:\\GIT_models\\git_leish_PE\\simulated data')
dat=data.matrix(read.csv('simulated data.csv',as.is=T))
nanos=ncol(dat)
nloc=nrow(dat)

#get dat complete
dat.complete=dat
for (i in 2:nanos){
  cond=dat.complete[,i-1]==1
  dat.complete[cond,i]=1
}
image(dat.complete)

#get distance
setwd('U:\\anaia\\derived data')
tmp=read.csv('matriz One_over_dist.csv',as.is=T)
rownames(tmp)=tmp$X
ind=which(colnames(tmp)=='X')
One.over.dist=tmp[,-ind]
soma=rowSums(One.over.dist)

#get invasion pressure
IP=matrix(NA,nloc,nanos)
for (i in 2:nanos){
  tmp=One.over.dist*matrix(dat.complete[,i-1],nloc,nloc,byrow=T)
  IP[,i]=apply(tmp,1,sum,na.rm=T)/soma  
}

#get covariates
setwd('U:\\anaia\\derived data\\covs')
z=list.files()
ind=grep('a1_',z); z1=z[ind]; z1
ncov=5#length(z1)
wmat=1
for (i in 1:ncov){
  tmp=read.csv(z1[i],as.is=T)
  wmat=cbind(wmat,as.numeric(unlist(tmp)))
}
nomes=z1[1:ncov]
nomes=gsub('a1_','',nomes); nomes=gsub('.csv','',nomes)
colnames(wmat)=c('interc',nomes)
wmat=data.matrix(wmat)

#useful stuff
Identifiers=matrix(1:(nloc*nanos),nloc,nanos)

#get covariates for detection probability
setwd('U:\\anaia\\derived data')
accesso=read.csv('derived assessoramento.csv',as.is=T)
ind=which(colnames(accesso)=='X'); accesso=accesso[,-ind]
capacit=read.csv('derived capacitacao.csv',as.is=T)
ind=which(colnames(capacit)=='X'); capacit=capacit[,-ind]
xmat=cbind(1,unlist(accesso),unlist(capacit)) #add 1's for intercept

#prior for alphas
sd.alpha=c(sqrt(10),rep(0.25,length(wmat)))
sd.psi=1

ngibbs=1000
nburn=ngibbs/2
res=gibbs.leish(wmat=wmat,xmat=xmat,dat.complete=dat.complete,
                dat=dat,ngibbs=ngibbs,nburn=nburn,IP=IP,
                sd.alpha=sd.alpha,sd.psi=sd.psi)

par(mfrow=c(2,1))
plot(res$gamma,type='l')
plot(res$psi,type='l')

z1.estim=res$z[nrow(res$z),]
betas.estim=colMeans(res$betas)
alpha.estim=colMeans(res$alpha)

