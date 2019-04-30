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
dat=data.matrix(read.csv('derived incidence.csv',as.is=T))
locs=dat[,'X']
ind=which(colnames(dat)=='X')
dat=dat[,-ind]
nanos=ncol(dat)
nloc=nrow(dat)

#get dat complete
dat.complete=dat
for (i in 2:nanos){
  cond=dat.complete[,i-1]==1
  dat.complete[cond,i]=1
}
image(dat.complete)

#get IP
setwd('U:\\anaia\\derived data')
tmp=read.csv('derived IP.csv',as.is=T)
sum(locs!=tmp$X) #check if order of rows is correct
rownames(tmp)=tmp$X
ind=which(colnames(tmp)=='X')
IP=data.matrix(tmp[,-ind])

#get covariates
setwd('U:\\anaia\\derived data\\covs')
z1=c('a1_Acessib_median.csv','a1_Formação florestal.csv',
     'a1_IDHM_2010.csv','a1_Urbana.csv')
ncov=length(z1)
wmat=1
for (i in 1:ncov){
  tmp=read.csv(z1[i],as.is=T)
  print(c(i,sum(locs!=tmp$X))) #check if order of rows is correct
  ind=which(colnames(tmp)=='X')
  tmp=tmp[,-ind]
  wmat=cbind(wmat,as.numeric(unlist(tmp)))
}
nomes=z1[1:ncov]
nomes=gsub('a1_','',nomes); nomes=gsub('.csv','',nomes)
colnames(wmat)=c('interc',nomes)
wmat=data.matrix(wmat)

cor(wmat[,-1])
apply(wmat,2,mean)
apply(wmat,2,sd)

#useful stuff
Identifiers=matrix(1:(nloc*nanos),nloc,nanos)

#get covariates for detection probability
setwd('U:\\anaia\\derived data')
accesso=read.csv('derived assessoramento.csv',as.is=T)
sum(locs!=accesso$X) #check if order of rows is correct
ind=which(colnames(accesso)=='X'); accesso=accesso[,-ind]
capacit=read.csv('derived capacitacao.csv',as.is=T)
sum(locs!=capacit$X) #check if order of rows is correct
ind=which(colnames(capacit)=='X'); capacit=capacit[,-ind]
notif=read.csv('derived notificacao.csv',as.is=T)
sum(floor(locs/10)!=notif$InputID) #check if order of rows is correct
notif1=matrix(notif$dist.notif,nloc,nanos)
#add distance to notificacao
xmat=cbind(1,unlist(accesso),unlist(capacit),matrix(notif1,nanos*nloc,1)) #add 1's for intercept
colnames(xmat)=c('interc','accesso','capacit','dist.notif')
xmat=data.matrix(xmat)

apply(xmat,2,range)
apply(xmat,2,mean)
apply(xmat,2,sd)

#prior for alphas
sd.alpha=c(1,rep(0.25,length(wmat)))
sd.psi=1

ngibbs=50000
nburn=ngibbs/2
res=gibbs.leish(wmat=wmat,xmat=xmat,dat.complete=dat.complete,
                dat=dat,ngibbs=ngibbs,nburn=nburn,IP=IP,
                sd.alpha=sd.alpha,sd.psi=sd.psi)

setwd('U:\\anaia\\results')
betas=res$betas; colnames(betas)=colnames(xmat)
alpha=res$alpha; colnames(alpha)=colnames(wmat)
write.csv(betas,'betas.csv',row.names=F)
write.csv(alpha,'alpha.csv',row.names=F)
write.csv(res$gamma,'gamma.csv',row.names=F)
write.csv(res$psi,'psi.csv',row.names=F)
seq1=seq(from=1,to=nrow(res$z),length.out=5000)
write.csv(res$z[seq1,],'z.csv',row.names=F)
