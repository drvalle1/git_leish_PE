rm(list=ls(all=TRUE))
set.seed(10)

setwd('U:\\anaia\\derived data\\covs')
tmp=read.csv('matriz One_over_dist.csv',as.is=T)
rownames(tmp)=tmp$X
ind=which(colnames(tmp)=='X')
One.over.dist=tmp[,-ind]

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
covs.invasion=c('Indice.Saneamento','Indice.Inst.Sanita','Cultivos.total','Formação.florestal') 
wmat=data.matrix(cbind(1,tmp[,covs.invasion]))

#get regression parameters
alpha=rep(0,ncol(wmat))
alpha[1]=-2; alpha[2]=0.5;alpha[3]=-0.5
alpha.true=alpha
betas=betas.true=c(-1,0.5,1)

#media
pmedia=matrix(wmat%*%alpha,nloc,nanos)
prob.invasion=exp(pmedia)/(1+exp(pmedia)); 
hist(prob.invasion)

#detection probability
assesso1=tmp$assessoramento;
capacit1=tmp$capacitacao
logit.p=matrix(betas[1]+betas[2]*assesso1+betas[3]*capacit1,nloc,nanos)
p=exp(logit.p)/(1+exp(logit.p)); 
p1=data.matrix(p); hist(p1)

#get latent z's
y1=z1=u1=matrix(NA,nloc,nanos)
gamma=gamma.true=0.1 #probability of occupancy in year 1
z1[,1]=rbinom(nloc,size=1,prob=gamma) #true invasion status
y1[,1]=rbinom(nloc,size=1,prob=z1[,1]*p1[,1]) #observed invasion status
y.complete=y1

#simulate forward
for (i in 2:nanos){
  z1[,i]=rbinom(nloc,size=1,prob=prob.invasion[,i-1])
  
  #once invaded, always invaded
  cond=z1[,i-1]==1
  z1[cond,i]=1
  
  #get observation
  y1[,i]=rbinom(nloc,size=1,prob=z1[,i]*p1[,i])
}

#look at simulated results
image(z1)
image(y1)
z1.true=z1

#export results
X=dat$X
colnames(y1)=paste0('inc',7:17)
setwd('U:\\GIT_models\\git_leish_PE\\simulated data')
write.csv(cbind(X,y1),'simulated data.csv',row.names=F)


