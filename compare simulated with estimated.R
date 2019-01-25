rango=range(c(alpha.true,alpha.estim))
plot(alpha.true,alpha.estim,xlim=rango,ylim=rango)
lines(rango,rango,col='red')

rango=range(c(betas.true,betas.estim))
plot(betas.true,betas.estim,xlim=rango,ylim=rango)
lines(rango,rango,col='red')

ind=unlist(ind.random)
dat=data.frame(estim=as.numeric(z1.estim[ind]),true=as.numeric(z1.true[ind]))
table(dat)