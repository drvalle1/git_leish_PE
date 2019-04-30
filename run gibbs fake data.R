rm(list=ls(all=TRUE))
library('Rcpp')
library('mvtnorm')
set.seed(4)

#read relevant functions
setwd('U:\\GIT_models\\git_leish_PE')
source('gibbs PE leish functions.R')
source('gibbs PE leish.R')
sourceCpp('aux_PE_leish.cpp')

#get incidence
setwd('U:\\GIT_models\\git_leish_PE\\simulated data')
dat0=read.csv('simulated data.csv',as.is=T)
ind=which(colnames(dat0)=='X')
dat=data.matrix(dat0[,-ind])
nanos=ncol(dat)
nloc=nrow(dat)

#get dat complete
dat.complete=dat
for (i in 2:nanos){
  cond=dat.complete[,i-1]==1
  dat.complete[cond,i]=1
}
image(dat.complete)

#get covariates
setwd('U:\\anaia\\derived data\\covs')
tmp=read.csv('covs combo standard.csv',as.is=T)
sum(rep(floor(dat0[,'X']/10),times=11)!=tmp$ID.IBGE) #basic checking
covs.invasion=c('Indice.Saneamento','Indice.Inst.Sanita','Cultivos.total','Formação.florestal')
wmat=data.matrix(cbind(1,tmp[,covs.invasion]))

#useful stuff
Identifiers=matrix(1:(nloc*nanos),nloc,nanos)

#get covariates for detection probability
assesso=tmp[,'assessoramento']
capacit=tmp[,'capacitacao']
xmat=cbind(1,assesso,capacit) #add 1's for intercept

#prior for alphas
sd.alpha=c(sqrt(10),rep(sqrt(10),length(wmat)))

ngibbs=1000
nburn=ngibbs/2
res=gibbs.leish(wmat=wmat,xmat=xmat,dat.complete=dat.complete,
                dat=dat,ngibbs=ngibbs,nburn=nburn,
                sd.alpha=sd.alpha)