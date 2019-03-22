gibbs.leish=function(wmat,xmat,dat.complete,dat,ngibbs,nburn,IP,sd.alpha,sd.psi){
  #useful values
  nloc=nrow(dat)
  nanos=ncol(dat)
  nparam=length(wmat)
  
  #initial values for parameters
  psi=0
  alpha=rep(0,ncol(wmat))
  betas=rep(0,ncol(xmat))
  z1=dat.complete
  medias=get.medias(wmat=wmat,alpha=alpha,nloc=nloc,nanos=nanos,nparam=nparam,psi=psi,IP=IP)
  gamma=0.5

  #get indicators for z's that are truly latent
  ind.random=list()
  for (i in 1:nloc){
    tmp=which(dat[i,]==1)
    tmp1=0:nanos
    if (dat[i,1]==1) tmp1=numeric()
    if (length(tmp) != 0 & dat[i,1]!=1) tmp1=0:(min(tmp)-1)
    ind.random[[i]]=tmp1
  }

  #stuff for gibbs sampler
  vec.z=matrix(NA,ngibbs,nloc*nanos)
  vec.betas=matrix(NA,ngibbs,ncol(xmat))
  vec.alpha=matrix(NA,ngibbs,ncol(wmat))
  vec.gamma=vec.psi=matrix(NA,ngibbs,1)
  jump1=list(betas=rep(1,ncol(xmat)),alpha=rep(1,length(alpha)),psi=0.2)
  accept1=list(betas=rep(0,ncol(xmat)),alpha=rep(0,length(alpha)),psi=0)
  accept.output=50

  #start gibbs sampler
  for (i in 1:ngibbs){
    print(i)
  
    z1=sample.z(nloc=nloc,nanos=nanos,ind.random=ind.random,z1=z1,
                medias=medias,xmat=xmat,betas=betas,dat=dat,gamma=gamma)
    # z1=z1.true
  
    tmp=sample.alpha(z1=z1,IP=IP,nanos=nanos,sd.alpha=sd.alpha,alpha=alpha,jump=jump1$alpha,
                     wmat=wmat,nloc=nloc,nparam=nparam,gamma=gamma,psi=psi)
    alpha=tmp$alpha
    accept1$alpha=accept1$alpha+tmp$accept
    # alpha=alpha.true
    medias=get.medias(wmat=wmat,alpha=alpha,nloc=nloc,nanos=nanos,nparam=nparam,
                      psi=psi,IP=IP)
    
    tmp=sample.betas(betas=betas,jump=jump1$betas,xmat=xmat,dat=dat,z1=z1)
    accept1$betas=accept1$betas+tmp$accept
    betas=tmp$betas
    # betas=betas.true
    
    gamma=sample.gamma(z1=z1,nloc=nloc)
    
    tmp=sample.psi(z1=z1,IP=IP,nanos=nanos,alpha=alpha,jump=jump1$psi,
                   wmat=wmat,nloc=nloc,nparam=nparam,gamma=gamma,
                   psi=psi,sd.psi=sd.psi)    
    accept1$psi=accept1$psi+tmp$accept
    psi=tmp$psi
    
    #adaptation MH algorithm
    if (i<nburn & i%%accept.output==0){
      k=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output)
      accept1=k$accept1
      jump1=k$jump1
    }
  
    #store results
    vec.z[i,]=z
    vec.betas[i,]=betas
    vec.alpha[i,]=alpha
    vec.gamma[i]=gamma
    vec.psi[i]=psi
  }

  seq1=nburn:ngibbs
  list(z=vec.z[seq1,],
       betas=vec.betas[seq1,],
       alpha=vec.alpha[seq1,],
       gamma=vec.gamma[seq1],
       psi=vec.psi[seq1])
}