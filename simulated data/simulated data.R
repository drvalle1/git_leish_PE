rm(list=ls(all=TRUE))
set.seed(1)

setwd('U:\\anaia\\derived data')
dist=read.csv('matriz distancia.csv',as.is=T)
rownames(dist)=dist$X
ind=which(colnames(dist)=='X')
dist=dist[,-ind]
diag(dist)=NA

dat=read.csv('derived incidence.csv',as.is=T)
ind=which(colnames(dat)=='X')
dat1=dat[,-ind]
nloc=nrow(dat1)
nanos=ncol(dat1)

#get covariates
setwd('U:\\anaia\\derived data\\covs')
z=list.files()
ind=grep('a1_',z); z1=z[ind]; z1
ncov=length(z1)
wmat=list()
for (i in 1:ncov){
  tmp=read.csv(z1[i],as.is=T)
  wmat[[i]]=tmp
}

#get covariates for detection probability
xmat=matrix(0,nloc,nanos)
colnames(xmat)=paste0('x',7:17)
xmat[,paste0('x',15:17)]=1
xmat=cbind(1,as.numeric(xmat)) #add 1's for intercept

#get regression parameters
alpha=rep(0,length(wmat)+1+1) #for intercept and for invasion pressure
alpha[1]=-0.5; alpha[3]=-1; alpha[length(alpha)]=0.5
alpha.true=alpha
betas=betas.true=c(-0.4,1)

#parte media
pmedia=alpha[1]
for (j in 1:ncov){
  pmedia=pmedia+alpha[j+1]*wmat[[j]]
}

#detection probability
logit.p=xmat%*%betas
p=exp(logit.p)/(1+exp(logit.p)); hist(p)
p1=matrix(p,nloc,nanos)

#get latent z's
y1=z1=u1=matrix(NA,nloc,nanos)
z1[,1]=rbinom(nloc,size=1,prob=0.1)
y1[,1]=rbinom(nloc,size=1,prob=z1[,1]*p1[,1])
y.complete=y1
soma=apply(1/dist,1,sum,na.rm=T)
for (i in 2:nanos){
  #get invasion pressure
  IP=(1/dist)*matrix(y.complete[,i-1],nloc,nloc,byrow=T)
  IP1=apply(IP,1,sum,na.rm=T)/soma
  
  #get mean
  media=pmedia[,i]+alpha[length(alpha)]*IP1
  prob=exp(media)/(1+exp(media)); hist(prob)
  z1[,i]=rbinom(nloc,size=1,prob=prob)
  
  #once invaded, always invaded
  cond=z1[,i-1]==1
  z1[cond,i]=1
  
  #get observation
  y.complete[,i]=y1[,i]=rbinom(nloc,size=1,prob=z1[,i]*p1[,i])
  cond=y.complete[,i-1]==1
  y.complete[cond,i]=1
}

#look at simulated results
image(z1)
image(y.complete)
image(y1)

#export results
setwd('U:\\GIT_models\\git_leish_PE\\simulated data')
write.csv(y1,'simulated data.csv',row.names=F)
write.csv(z1,'simulated status.csv',row.names=F)


