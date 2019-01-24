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
sample.u=function(z1,nanos,nloc,Identifiers,medias){
  ind0=which(z1==0)
  ind0=ind0[!(ind0%in%Identifiers[,1])] #we are not interested in u's at time 1
  ind1=IdentifLead1(z=z1, nanos=nanos, nlocs=nloc,Identifiers=Identifiers) 
  ind1=ind1[ind1!=0] #eliminate locations that have zeroes for all years
  ind1=ind1[!(ind1%in%Identifiers[,1])] #we are not interested in u's at time 1
  u=matrix(NA,nloc,nanos)
  u[ind0]=tnorm(length(ind0),lo=-Inf,hi=0,mu=medias[ind0],sig=1)
  u[ind1]=tnorm(length(ind1),lo=0 ,hi=Inf,mu=medias[ind1],sig=1)
  u
}
#----------------------------------
sample.z=function(nloc,nanos,ind.random,z1,OneOverDist,SomaOneOverDist,wmat,alpha,
                  xmat,betas,dat){
  z1.old=z1.new=z1
  for (i in 1:nloc){
    z1.new=z1.old
    if (length(ind.random[[i]])!=0){ #if dat[i,] starts with 1, then ind.random[[i]]=numeric(0) and there are no random z's in that location
      #propose new set of z's for location i
      k=sample(ind.random[[i]],size=1)
      if (k!=nanos+1) z1.new[i,k:nanos]=1
      if (k!=1) z1.new[i,1:(k-1)]=0 #notice that if k=anos+1, then this allows for the possibility that the site was never invaded

      #calculate loglikelihood
      llk.old=get.llk(xmat=xmat,betas=betas,dat=dat,z1=z1.old,nloc=nloc,nanos=nanos)
      llk.new=get.llk(xmat=xmat,betas=betas,dat=dat,z1=z1.new,nloc=nloc,nanos=nanos)
      
      #calculate process probability
      IP.old=CalcInvasionPressure(z=z1.old, OneOverDist=OneOverDist, nanos=nanos, nlocs=nloc,
                                  SomaOneOverDist=SomaOneOverDist)
      
      IP.new=CalcInvasionPressure(z=z1.new, OneOverDist=OneOverDist, nanos=nanos, nlocs=nloc,
                                  SomaOneOverDist=SomaOneOverDist)

      medias.old=get.medias(wmat=wmat,alpha=alpha,nloc=nloc,nanos=nanos,IP=IP.old)
      medias.new=get.medias(wmat=wmat,alpha=alpha,nloc=nloc,nanos=nanos,IP=IP.new)
      
      lpprob.old=get.lprocess.prob(nanos=nanos,z1=z1.old,nloc=nloc,media=medias.old)
      lpprob.new=get.lprocess.prob(nanos=nanos,z1=z1.new,nloc=nloc,media=medias.new)
      
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
get.lprocess.prob=function(nanos,z1,nloc,media){
  lprob=matrix(NA,nloc,nanos)
  
  #for first time step
  cond=z1[,1]==1
  lprob[cond,1]=log(0.1)
  lprob[!cond,1]=log(0.9)
    
  #for other time steps
  for (i in 2:nanos){
    cond1=z1[,i-1]==0 & z1[,i]==1
    lprob[cond1,i]=pnorm(media[cond1,i],log=T)
    cond0=z1[,i-1]==0 & z1[,i]==0
    lprob[cond0,i]=pnorm(media[cond0,i],log=T,lower.tail=F)
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
sample.alpha=function(u1,wmat.design,IP,nanos,invT.alpha){
  #create design matrix
  IP1=cbind(NA,IP[,-nanos]) #mean at year 2 uses invasion pressure for year 1
  wmat1=cbind(1,wmat.design,as.numeric(IP1))
  
  #eliminate missing u1 (either at year 1 or after an invasion has already occurred)
  ind=which(!is.na(u1))
  wmat2=wmat1[ind,]
  u2=u1[ind]
  
  #sample alpha
  prec=t(wmat2)%*%wmat2+invT.alpha
  var1=solve(prec)
  pmedia=t(wmat2)%*%u2
  t(rmvnorm(1,mean=var1%*%pmedia,sigma=var1))
}
#----------------------------------
get.medias=function(wmat,alpha,nloc,nanos,IP){
  medias=matrix(alpha[1],nloc,nanos)
  for (i in 1:length(wmat)){
    medias=medias+wmat[[i]]*alpha[i+1]
  }
  medias+alpha[length(alpha)]*cbind(NA,IP[,-nanos]) #mean for 2nd year depends on IP from year 1
}