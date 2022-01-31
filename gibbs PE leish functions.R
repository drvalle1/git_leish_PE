tnorm <- function(n,lo,hi,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo[z == -Inf]
  z[z == Inf]   <- hi[z == Inf]
  z
}
#-------------------------------
acceptMH <- function(p0,p1,x0,x1,BLOCK){   #accept for M, M-H
  # if BLOCK, then accept as a block,
  # otherwise, accept individually
  
  nz           <- length(x0)  #no. to accept
  if(BLOCK) nz <- 1
  
  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a)
  
  if(BLOCK & length(keep) > 0) x0 <- x1
  if(!BLOCK)                   x0[keep] <- x1[keep]           
  accept <- length(keep)        
  
  list(x = x0, accept = accept)
}
#----------------------------------
print.adapt = function(accept1z,jump1z,accept.output){
  accept1=accept1z; jump1=jump1z; 
  
  for (k in 1:length(accept1)){
    z=accept1[[k]]/accept.output
    print(names(accept1)[k])
    print(mean(z)); print(mean(jump1[[k]]))
  }
  
  for (k in 1:length(jump1)){
    cond=(accept1[[k]]/accept.output)>0.4 & jump1[[k]]<100
    jump1[[k]][cond] = jump1[[k]][cond]*2       
    cond=(accept1[[k]]/accept.output)<0.2 & jump1[[k]]>0.01
    jump1[[k]][cond] = jump1[[k]][cond]*0.5
    accept1[[k]][]=0
  }
  
  return(list(jump1=jump1,accept1=accept1))
}
#----------------------------------
sample.z=function(nloc,nanos,ind.random.loc,ind.random.zeroes,z1,medias,xmat,betas,dat,gamma1){
  z1.old=z1.new=z1
  for (i in ind.random.loc){ #if dat[i,1]=1 then don't do anything as there are no random z's in that location
    #propose new set of z's for location i
    tmp=ind.random.zeroes[[i]]
    k=sample(tmp,size=1)
    if (k==0)        z1.new[i,tmp[-1]]=0 #everything is zero
    if (k!=0)        z1.new[i,k:nanos]=1
    if (k> 1)        z1.new[i,1:(k-1)]=0 
  }
  # rowSums(z1.old!=z1.new)
  
  #calculate loglikelihood
  llk.old=get.llk(xmat=xmat,betas=betas,dat=dat,z1=z1.old,nloc=nloc,nanos=nanos)
  llk.new=get.llk(xmat=xmat,betas=betas,dat=dat,z1=z1.new,nloc=nloc,nanos=nanos)
      
  #calculate process probability
  lpprob.old=get.lprocess.prob(nanos=nanos,z1=z1.old,nloc=nloc,media=medias,gamma1=gamma1)
  lpprob.new=get.lprocess.prob(nanos=nanos,z1=z1.new,nloc=nloc,media=medias,gamma1=gamma1)
      
  #calculate final probability of changing to new state
  lp.old=lpprob.old+llk.old
  lp.new=lpprob.new+llk.new
  max1=apply(cbind(lp.old,lp.new),1,max)
  lp.old1=lp.old-max1
  lp.new1=lp.new-max1
  p.old=exp(lp.old1)
  p.new=exp(lp.new1)
  prob=p.new/(p.old+p.new)
      
  #accept or reject move
  cond=runif(nloc)<prob
  z1.old[cond,]=z1.new[cond,]
  z1.old
}
#----------------------------------
get.lprocess.prob=function(nanos,z1,nloc,media,gamma1){
  lprob=rep(0,nloc)
  
  #for first time step
  cond=z1[,1]==1
  lprob[cond]=log(gamma1)
  lprob[!cond]=log(1-gamma1)
    
  #for other time steps
  prob=exp(media)/(1+exp(media))
  for (i in 2:nanos){
    cond1=z1[,i-1]==0 & z1[,i]==1
    lprob[cond1]=lprob[cond1]+log(prob[cond1,i-1])
    cond0=z1[,i-1]==0 & z1[,i]==0
    lprob[cond0]=lprob[cond0]+log(1-prob[cond0,i-1])
  }
  lprob
}
#----------------------------------
get.llk=function(xmat,betas,dat,z1,nloc,nanos){
  media=xmat%*%betas
  tmp=exp(media); 
  prob=matrix(tmp/(1+tmp),nloc,nanos)
  rowSums(dbinom(dat,size=1,prob=z1*prob,log=T))
}
#----------------------------------
sample.betas=function(betas,jump,xmat,dat,z1,sd.betas){
  betas.old=betas.new=betas
  nparam=length(betas)
  z2=as.numeric(z1)
  sd1=sd.betas
  for (i in 1:nparam){
    betas.new=betas.old
    betas.new[i]=rnorm(1,mean=betas.old[i],sd=jump[i])
    llk.old=get.llk(xmat,betas.old,dat,z1,nloc,nanos)
    llk.new=get.llk(xmat,betas.new,dat,z1,nloc,nanos)
    prior.old=dnorm(betas.old[i],mean=0,sd=sd1[i],log=T)
    prior.new=dnorm(betas.new[i],mean=0,sd=sd1[i],log=T)
    k=acceptMH(p0=sum(llk.old)+prior.old,
               p1=sum(llk.new)+prior.new,
               x0=betas.old[i],
               x1=betas.new[i],BLOCK=F)
    betas.old[i]=k$x
  }
  list(betas=betas.old,accept=betas.old!=betas)
}
#----------------------------------
sample.alpha=function(z1,nanos,sd.alpha,alpha,jump,wmat,nloc,nparam,gamma1){
  alpha.old=alpha
  for (i in 1:length(alpha)){
    alpha.new=alpha.old
    alpha.new[i]=rnorm(1,mean=alpha.old[i],sd=jump[i])
    
    medias.old=get.medias(wmat=wmat,alpha=alpha.old,nloc=nloc,nanos=nanos)
    medias.new=get.medias(wmat=wmat,alpha=alpha.new,nloc=nloc,nanos=nanos)
    lpprob.old=get.lprocess.prob(nanos=nanos,z1=z1,nloc=nloc,media=medias.old,gamma1=gamma1)
    lpprob.new=get.lprocess.prob(nanos=nanos,z1=z1,nloc=nloc,media=medias.new,gamma1=gamma1)
    
    lprior.old=dnorm(alpha.old[i],mean=0,sd=sd.alpha[i],log=T)
    lprior.new=dnorm(alpha.new[i],mean=0,sd=sd.alpha[i],log=T)
    
    k=acceptMH(p0=sum(lpprob.old)+lprior.old,
               p1=sum(lpprob.new)+lprior.new,
               x0=alpha.old[i],x1=alpha.new[i],BLOCK=F)
    alpha.old[i]=k$x
  }
  list(alpha=alpha.old,accept=alpha!=alpha.old)
}
#----------------------------------
sample.gamma=function(z1,nloc){
  soma=sum(z1[,1])
  a1=soma+1
  b1=nloc-soma+1
  rbeta(1,a1,b1)
}
#----------------------------------
get.medias=function(wmat,alpha,nloc,nanos){
  matrix(wmat%*%alpha,nloc,nanos)
}