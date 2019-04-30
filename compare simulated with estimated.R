k=res$alpha
z=cor(k)
z[z > -0.5 & z < 0.5]=NA
diag(z)=NA; z

par(mfrow=c(1,1))
plot(res$gamma,type='l')

z1.estim=res$z[nrow(res$z),]
betas.estim=colMeans(res$betas)
alpha.estim=colMeans(res$alpha)

alpha.true1=alpha.true
# alpha.estim[1]=alpha.estim[1]-(alpha.estim[length(alpha.estim)]*media.ip/sd.ip)
# alpha.estim[length(alpha.estim)]=alpha.estim[length(alpha.estim)]/sd.ip
rango=range(c(alpha.true1,alpha.estim))
plot(alpha.true1,alpha.estim,xlim=rango,ylim=rango)
lines(rango,rango,col='red')

rango=range(c(betas.true,betas.estim))
plot(betas.true,betas.estim,xlim=rango,ylim=rango)
lines(rango,rango,col='red')

#look at z's
ztmp=data.frame(z.estim=apply(res$z,2,mean),
                z.true=matrix(z1.true,nloc*nanos,1))
ind=which(dat[,-1]==1)
ztmp1=ztmp[-ind,]
boxplot(z.estim~z.true,data=ztmp1)
