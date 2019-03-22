rm(list=ls(all=TRUE))
set.seed(10)

setwd('U:\\anaia\\derived data')
tmp=read.csv('matriz One_over_dist.csv',as.is=T)
rownames(tmp)=tmp$X
ind=which(colnames(tmp)=='X')
One.over.dist=tmp[,-ind]

dat=read.csv('derived incidence.csv',as.is=T)
ind=which(colnames(dat)=='X')
dat1=dat[,-ind]
locs=dat[,ind]
nloc=nrow(dat1)
nanos=ncol(dat1)

#get capacitacao
capacit=read.csv('derived capacitacao.csv',as.is=T)
sum(locs!=capacit$X)
ind=which(colnames(capacit)=='X')
capacit1=capacit[,-ind]

#get assessoramento
assesso=read.csv('derived assessoramento.csv',as.is=T)
sum(locs!=assesso$X)
ind=which(colnames(assesso)=='X')
assesso1=assesso[,-ind]

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

#get regression parameters
alpha=rep(0,ncol(wmat))
alpha[1]=-1; alpha[2]=0.5;
alpha.true=alpha
betas=betas.true=c(-0.8,0,2)
psi=psi.true=1

#parte media
pmedia=matrix(wmat%*%alpha,nloc,nanos)

#detection probability
logit.p=betas[1]+betas[2]*assesso1+betas[3]*capacit1
p=exp(logit.p)/(1+exp(logit.p)); 
p1=data.matrix(p); hist(p1)

#get latent z's
y1=z1=u1=matrix(NA,nloc,nanos)
gamma=gamma.true=0.1
z1[,1]=rbinom(nloc,size=1,prob=gamma)
y1[,1]=rbinom(nloc,size=1,prob=z1[,1]*p1[,1])
y.complete=y1
soma=apply(One.over.dist,1,sum,na.rm=T)
IP.true=media.true=matrix(NA,nloc,nanos)
for (i in 2:nanos){
  #get invasion pressure
  tmp=One.over.dist*matrix(y.complete[,i-1],nloc,nloc,byrow=T)
  IP.true[,i]=apply(tmp,1,sum,na.rm=T)/soma

  #get mean
  media.true[,i]=pmedia[,i]+psi*IP.true[,i]
  prob=exp(media.true[,i])/(1+exp(media.true[,i])); hist(prob)
  z1[,i]=rbinom(nloc,size=1,prob=prob)
  
  #once invaded, always invaded
  cond=z1[,i-1]==1
  z1[cond,i]=1
  
  #get observation
  y.complete[,i]=y1[,i]=rbinom(nloc,size=1,prob=z1[,i]*p1[,i])
  cond=y.complete[,i-1]==1
  y.complete[cond,i]=1
}
range(IP.true,na.rm=T)

#look at simulated results
image(z1)
image(y.complete)
image(y1)
z1.true=z1

#export results
setwd('U:\\GIT_models\\git_leish_PE\\simulated data')
write.csv(y1,'simulated data.csv',row.names=F)
write.csv(z1,'simulated status.csv',row.names=F)


