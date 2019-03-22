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
sample.z=function(nloc,nanos,ind.random,z1,medias,xmat,betas,dat,gamma){
  z1.old=z1.new=z1
  for (i in 1:nloc){
    z1.new=z1.old
    if (length(ind.random[[i]])!=0){ #if dat[i,] starts with 1, then ind.random[[i]]=numeric(0) and there are no random z's in that location
      #propose new set of z's for location i
      k=sample(ind.random[[i]],size=1)
      if (k==0)        z1.new[i,ind.random[[i]]]=0 
      if (k!=0)        z1.new[i,k:nanos]=1
      if (k> 1)        z1.new[i,1:(k-1)]=0 

      #calculate loglikelihood
      llk.old=get.llk(xmat=xmat,betas=betas,dat=dat,z1=z1.old,nloc=nloc,nanos=nanos)
      llk.new=get.llk(xmat=xmat,betas=betas,dat=dat,z1=z1.new,nloc=nloc,nanos=nanos)
      
      #calculate process probability
      lpprob.old=get.lprocess.prob(nanos=nanos,z1=z1.old,nloc=nloc,media=medias,gamma=gamma)
      lpprob.new=get.lprocess.prob(nanos=nanos,z1=z1.new,nloc=nloc,media=medias,gamma=gamma)
      
      #calculate final probability of changing to new state
      lp.old=lpprob.old+llk.old
      lp.new=lpprob.new+llk.new
      max1=max(c(lp.old,lp.new))
      lp.old1=lp.old-max1
      lp.new1=lp.new-max1
      p.old=exp(lp.old1)
      p.new=exp(lp.new1)
      prob=p.new/(p.old+p.new)
      
      #accept or reject move
      k1=rbinom(1,size=1,prob=prob)
      if (k1==1) z1.old[i,]=z1.new[i,]
    }
  }
  z1.old
}
#----------------------------------
get.lprocess.prob=function(nanos,z1,nloc,media,gamma){
  lprob=matrix(NA,nloc,nanos)
  
  #for first time step
  cond=z1[,1]==1
  lprob[cond,1]=log(gamma)
  lprob[!cond,1]=log(1-gamma)
    
  #for other time steps
  prob=exp(media)/(1+exp(media))
  for (i in 2:nanos){
    cond1=z1[,i-1]==0 & z1[,i]==1
    lprob[cond1,i]=log(prob[cond1,i])
    cond0=z1[,i-1]==0 & z1[,i]==0
    lprob[cond0,i]=log(1-prob[cond0,i])
  }
  sum(lprob,na.rm=T)
}
#----------------------------------
get.llk=function(xmat,betas,dat,z1,nloc,nanos){
  media=xmat%*%betas
  tmp=exp(media); 
  prob=matrix(tmp/(1+tmp),nloc,nanos)
  sum(dbinom(dat,size=1,prob=z1*prob,log=T))
}
#----------------------------------
sample.betas=function(betas,jump,xmat,dat,z1){
  betas.old=betas.new=betas
  nparam=length(betas)
  z2=as.numeric(z1)
  sd1=c(sqrt(10),rep(1,nparam-1))
  for (i in 1:nparam){
    betas.new=betas.old
    betas.new[i]=rnorm(1,mean=betas.old[i],sd=jump[i])
    llk.old=get.llk(xmat,betas.old,dat,z1,nloc,nanos)
    llk.new=get.llk(xmat,betas.new,dat,z1,nloc,nanos)
    prior.old=dnorm(betas.old[i],mean=0,sd=sd1[i],log=T)
    prior.new=dnorm(betas.new[i],mean=0,sd=sd1[i],log=T)
    k=acceptMH(p0=llk.old+prior.old,
               p1=llk.new+prior.new,
               x0=betas.old[i],
               x1=betas.new[i],BLOCK=F)
    betas.old[i]=k$x
  }
  list(betas=betas.old,accept=betas.old!=betas)
}
#----------------------------------
sample.alpha=function(z1,IP,nanos,sd.alpha,alpha,jump,wmat,nloc,nparam,gamma,psi){
  alpha.old=alpha
  for (i in 1:length(alpha)){
    alpha.new=alpha.old
    alpha.new[i]=rnorm(1,mean=alpha.old[i],sd=jump[i])
    medias.old=get.medias(wmat=wmat,alpha=alpha.old,nloc=nloc,nanos=nanos,nparam=nparam,
                          psi=psi,IP=IP)
    medias.new=get.medias(wmat=wmat,alpha=alpha.new,nloc=nloc,nanos=nanos,nparam=nparam,
                          psi=psi,IP=IP)
    
    lpprob.old=get.lprocess.prob(nanos=nanos,z1=z1,nloc=nloc,media=medias.old,gamma=gamma)
    lpprob.new=get.lprocess.prob(nanos=nanos,z1=z1,nloc=nloc,media=medias.new,gamma=gamma)
    
    lprior.old=dnorm(alpha.old[i],mean=0,sd=sd.alpha[i],log=T)
    lprior.new=dnorm(alpha.new[i],mean=0,sd=sd.alpha[i],log=T)
    k=acceptMH(p0=lpprob.old+lprior.old,p1=lpprob.new+lprior.new,
               x0=alpha.old[i],x1=alpha.new[i],BLOCK=F)
    alpha.old[i]=k$x
  }
  list(alpha=alpha.old,accept=alpha!=alpha.old)
}
#----------------------------------
get.medias=function(wmat,alpha,nloc,nanos,nparam,psi,IP){
  matrix(wmat%*%alpha,nloc,nanos)+psi*IP
}
#----------------------------------
sample.gamma=function(z1,nloc){
  soma=sum(z1[,1])
  a1=soma+1
  b1=nloc-soma+1
  rbeta(1,a1,b1)
}
#----------------------------------
sample.psi=function(z1,IP,nanos,alpha,jump,wmat,nloc,nparam,gamma,psi,sd.psi){
  psi.old=psi
  psi.new=rnorm(1,mean=psi.old,sd=jump)
  medias.old=get.medias(wmat=wmat,alpha=alpha,nloc=nloc,nanos=nanos,nparam=nparam,
                          psi=psi.old,IP=IP)
  medias.new=get.medias(wmat=wmat,alpha=alpha,nloc=nloc,nanos=nanos,nparam=nparam,
                          psi=psi.new,IP=IP)
    
  lpprob.old=get.lprocess.prob(nanos=nanos,z1=z1,nloc=nloc,media=medias.old,gamma=gamma)
  lpprob.new=get.lprocess.prob(nanos=nanos,z1=z1,nloc=nloc,media=medias.new,gamma=gamma)
    
  lprior.old=dnorm(psi.old,mean=0,sd=sd.psi,log=T)
  lprior.new=dnorm(psi.new,mean=0,sd=sd.psi,log=T)
  k=acceptMH(p0=lpprob.old+lprior.old,p1=lpprob.new+lprior.new,
             x0=psi.old,x1=psi.new,BLOCK=F)
  list(psi=k$x,accept=k$x!=psi.old)
}