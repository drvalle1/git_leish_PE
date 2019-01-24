rm(list=ls(all=TRUE))
library('Rcpp')

#read relevant functions
setwd('U:\\GIT_models\\git_leish_PE')
sourceCpp('aux_PE_leish.cpp')

#get incidence
setwd('U:\\anaia\\simulated data')
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

IP=CalcInvasionPressure(z=dat, OneOverDist=OneOverDist, 
                        nanos=nanos, nlocs=nloc,
                        SomaOneOverDist=SomaOneOverDist)

IP1=matrix(NA,nloc,nanos)
for (i in 1:nanos){
  z1=matrix(dat[,i],nloc,nloc,byrow=T)
  tmp=rowSums(OneOverDist*z1)
  tmp1=tmp/SomaOneOverDist
  IP1[,i]=tmp1
}

plot(IP,IP1)
rango=c(-10,10)
lines(rango,rango,col='red')