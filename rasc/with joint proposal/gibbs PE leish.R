gibbs.leish=function(wmat,xmat,dat.complete,dat,ngibbs,nburn,sd.alpha){
  #useful values
  nloc=nrow(dat)
  nanos=ncol(dat)
  nparam=length(wmat)
  
  #initial values for parameters
  alpha=rep(0,ncol(wmat))
  betas=rep(0,ncol(xmat))
  z1=dat.complete
  medias=get.medias(wmat=wmat,alpha=alpha,nloc=nloc,nanos=nanos)
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
  vec.gamma=matrix(NA,ngibbs,1)
  jump1=list(betas=rep(1,ncol(xmat)),alpha=0.01)
  accept1=list(betas=rep(0,ncol(xmat)),alpha=0)
  accept.output=50
  cov1=diag(1,length(alpha))

  #start gibbs sampler
  for (i in 1:ngibbs){
    print(i)
  
    z1=sample.z(nloc=nloc,nanos=nanos,ind.random=ind.random,z1=z1,
                medias=medias,xmat=xmat,betas=betas,dat=dat,gamma=gamma)
    # z1=z1.true
  
    tmp=sample.alpha(z1=z1,IP=IP,nanos=nanos,sd.alpha=sd.alpha,alpha=alpha,jump=jump1$alpha,
                     wmat=wmat,nloc=nloc,nparam=nparam,gamma=gamma,cov1=cov1)
    alpha=tmp$alpha
    accept1$alpha=accept1$alpha+tmp$accept
    # alpha=alpha.true
    medias=get.medias(wmat=wmat,alpha=alpha,nloc=nloc,nanos=nanos)
    
    tmp=sample.betas(betas=betas,jump=jump1$betas,xmat=xmat,dat=dat,z1=z1)
    accept1$betas=accept1$betas+tmp$accept
    betas=tmp$betas
    # betas=betas.true
    
    gamma=sample.gamma(z1=z1,nloc=nloc)
    
    #adaptation MH algorithm
    if (i<nburn & i%%accept.output==0){
      cov1=var(vec.alpha[(i-accept.output+1):(i-1),])
      k=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output)
      accept1=k$accept1
      jump1=k$jump1
      jump1$alpha=1
    }
  
    #store results
    vec.z[i,]=z1
    vec.betas[i,]=betas
    vec.alpha[i,]=alpha
    vec.gamma[i]=gamma
  }

  seq1=nburn:ngibbs
  list(z=vec.z[seq1,],
       betas=vec.betas[seq1,],
       alpha=vec.alpha[seq1,],
       gamma=vec.gamma[seq1])
}