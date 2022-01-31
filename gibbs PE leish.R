gibbs.leish=function(wmat,xmat,dat.complete,dat,ngibbs,nburn,sd.alpha,sd.betas){
  #useful values
  nloc=nrow(dat)
  nanos=ncol(dat)
  nparam=length(wmat)
  
  #initial values for parameters
  alpha=rep(0,ncol(wmat))
  betas=rep(0,ncol(xmat))
  gamma1=0.1
  z1=dat.complete

  #which locations were not already infected in year 1?
  ind.random.loc=which(dat[,1]==0)
  
  #get indicators for z's that are truly latent
  ind.random.zeroes=list()
  for (i in ind.random.loc){
    cond=sum(dat[i,]==0)==nanos #all zeroes
    if (cond) ind.random.zeroes[[i]]=c(0,1:nanos)
    if (!cond) {
      min1=min(which(dat[i,]==1))
      ind.random.zeroes[[i]]=c(0,1:(min1-1))
    }
  }

  #stuff for gibbs sampler
  vec.z=matrix(NA,ngibbs,nloc*nanos)
  vec.betas=matrix(NA,ngibbs,ncol(xmat))
  vec.alpha=matrix(NA,ngibbs,ncol(wmat))
  vec.gamma1=matrix(NA,ngibbs,1)
  jump1=list(betas=rep(1,ncol(xmat)),alpha=rep(0.1,length(alpha)))
  accept1=list(betas=rep(0,ncol(xmat)),alpha=rep(0,length(alpha)))
  accept.output=50

  #start gibbs sampler
  for (i in 1:ngibbs){
    print(i)
  
    medias=get.medias(wmat=wmat,alpha=alpha,nloc=nloc,nanos=nanos)
    z1=sample.z(nloc=nloc,nanos=nanos,z1=z1,
                ind.random.zeroes=ind.random.zeroes,ind.random.loc=ind.random.loc,
                medias=medias,xmat=xmat,betas=betas,dat=dat,gamma1=gamma1)
    # z1=z1.true
  
    tmp=sample.alpha(z1=z1,nanos=nanos,sd.alpha=sd.alpha,alpha=alpha,jump=jump1$alpha,
                     wmat=wmat,nloc=nloc,nparam=nparam,gamma1=gamma1)
    alpha=tmp$alpha
    accept1$alpha=accept1$alpha+tmp$accept
    # alpha=alpha.true
    
    tmp=sample.betas(betas=betas,jump=jump1$betas,xmat=xmat,dat=dat,z1=z1,
                     sd.betas=sd.betas)
    accept1$betas=accept1$betas+tmp$accept
    betas=tmp$betas
    # betas=betas.true
    
    gamma1=sample.gamma(z1=z1,nloc=nloc)
    # gamma1=gamma1.true
    
    #adaptation MH algorithm
    if (i<nburn & i%%accept.output==0){
      k=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output)
      accept1=k$accept1
      jump1=k$jump1
    }
  
    #store results
    vec.z[i,]=z1
    vec.betas[i,]=betas
    vec.alpha[i,]=alpha
    vec.gamma1[i]=gamma1
  }

  seq1=nburn:ngibbs
  list(z=vec.z[seq1,],
       betas=vec.betas[seq1,],
       alpha=vec.alpha[seq1,],
       gamma1=vec.gamma1[seq1])
}