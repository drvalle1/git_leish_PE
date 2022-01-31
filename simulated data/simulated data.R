rm(list=ls(all=TRUE))
set.seed(110)

setwd('U:\\anaia\\derived data')
dat=read.csv('derived incidence.csv',as.is=T)
ind=which(colnames(dat)=='X')
dat1=dat[,-ind]
locs=dat[,ind]
nloc=nrow(dat1)
nanos=ncol(dat1)

#get covariates
setwd('U:\\anaia\\derived data\\covs')
tmp=read.csv('covs combo standard.csv',as.is=T)
sum(rep(floor(dat$X/10),times=11)!=tmp$ID.IBGE) #basic checking
covs.invasion=c('Indice.Saneamento','Indice.Inst.Sanita',
                'Cultivos.total','Formação.florestal') 
wmat=data.matrix(cbind(1,tmp[,covs.invasion]))
covs.detect=c('assessoramento','capacitacao','dist_sede')
# covs.detect='dist_sede'
xmat=data.matrix(cbind(1,tmp[,covs.detect]))

#get regression parameters
alpha=rep(0,ncol(wmat))
alpha[1]=-1; alpha[2]=0.5;alpha[3]=-0.5
alpha.true=alpha
betas=c(-0.5,0,0.4,0.5)
betas.true=betas

#media
pmedia=matrix(wmat%*%alpha,nloc,nanos)
prob.invasion=exp(pmedia)/(1+exp(pmedia)); 
hist(prob.invasion)

#get latent z's
y1=z1=u1=matrix(NA,nloc,nanos)
prob.z1=gamma1.true=0.23;prob.z1
z1[,1]=rbinom(nloc,size=1,prob=prob.z1) #true invasion status

#simulate forward
for (i in 2:nanos){
  z1[,i]=rbinom(nloc,size=1,prob=prob.invasion[,i-1])
  
  #once invaded, always invaded
  cond=z1[,i-1]==1
  z1[cond,i]=1
}
image(z1)

#detection probability
logit.p=xmat%*%betas
p=exp(logit.p)/(1+exp(logit.p)); 
prob.detect=matrix(p,nloc,nanos); hist(prob.detect)

#get observation
tmp=rbinom(nloc*nanos,size=1,prob=z1*prob.detect)
y1=matrix(tmp,nloc,nanos)
image(y1)

#look at simulated results
image(z1)
image(y1)
z1.true=z1

#export results
X=dat$X
colnames(y1)=paste0('inc',7:17)
setwd('U:\\GIT_models\\git_leish_PE\\simulated data')
write.csv(cbind(X,y1),'simulated data.csv',row.names=F)


