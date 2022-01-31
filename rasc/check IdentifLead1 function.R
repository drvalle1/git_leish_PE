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

Identifiers=matrix(1:(nloc*nanos),nloc,nanos)
res=IdentifLead1(z=dat, nanos=nanos, nlocs=nloc,Identifiers=Identifiers) 

res1=rep(0,nloc)
for (i in 1:nloc){
  for (j in 2:nanos){
    if (dat[i,j]==1 & dat[i,j-1]==0) res1[i]=Identifiers[i,j]    
  }
}
plot(res,res1)
unique(res-res1)