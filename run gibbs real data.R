rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(4)

#read relevant functions
setwd('U:\\GIT_models\\git_leish_PE')
source('gibbs PE leish functions.R')
source('gibbs PE leish.R')
sourceCpp('aux_PE_leish.cpp')

#get incidence
setwd('U:\\anaia\\derived data')
dat0=data.matrix(read.csv('derived incidence.csv',as.is=T))
locs=dat0[,'X']
ind=which(colnames(dat0)=='X')
dat=dat0[,-ind]
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
covs.invasion=c('Indice.Saneamento','Indice.Inst.Sanita','Cultivos.total','Forma��o.florestal')
wmat=data.matrix(cbind(1,tmp[,covs.invasion]))

#useful stuff
Identifiers=matrix(1:(nloc*nanos),nloc,nanos)

#get covariates for detection probability
assesso=tmp[,'assessoramento']
capacit=tmp[,'capacitacao']
xmat=cbind(1,assesso,capacit) #add 1's for intercept

#prior for alphas
sd.alpha=c(sqrt(10),rep(sqrt(10),length(wmat)))

ngibbs=10000
nburn=ngibbs/2
res=gibbs.leish(wmat=wmat,xmat=xmat,dat.complete=dat.complete,
                dat=dat,ngibbs=ngibbs,nburn=nburn,
                sd.alpha=sd.alpha)

setwd('U:\\anaia\\results')
betas=res$betas; colnames(betas)=colnames(xmat)
alpha=res$alpha; colnames(alpha)=colnames(wmat)
write.csv(betas,'betas.csv',row.names=F)
write.csv(alpha,'alpha.csv',row.names=F)
write.csv(res$gamma,'gamma.csv',row.names=F)
seq1=seq(from=1,to=nrow(res$z),length.out=5000)
write.csv(res$z[seq1,],'z.csv',row.names=F)
