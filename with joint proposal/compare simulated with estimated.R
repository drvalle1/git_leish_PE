k=res$alpha
cor(k)

par(mfrow=c(1,1))
plot(res$gamma,type='l')

z1.estim=res$z[nrow(res$z),]
betas.estim=colMeans(res$betas)
alpha.estim=colMeans(res$alpha)

alpha.true1=c(alpha.true,psi.true)
rango=range(c(alpha.true1,alpha.estim))
plot(alpha.true1,alpha.estim,xlim=rango,ylim=rango)
lines(rango,rango,col='red')

rango=range(c(betas.true,betas.estim))
plot(betas.true,betas.estim,xlim=rango,ylim=rango)
lines(rango,rango,col='red')

nsim=1000
res=rep(NA,nsim)
for (i in 1:nsim){
  tmp=sample.psi(z1=z1.true,IP=IP.true,nanos=nanos,alpha=alpha.true,jump=0.1,
                 wmat=wmat,nloc=nloc,nparam=nparam,gamma=gamma.true,
                 psi=psi,sd.psi=100)  
  psi=tmp$psi
  res[i]=psi
}
plot(res,type='l')