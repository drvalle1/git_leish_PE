rm(list=ls())

setwd('U:\\GIT_models\\git_leish_PE\\results')
betas.estim=read.csv('betas.csv')
alpha.estim=read.csv('alpha.csv')
gamma.estim=read.csv('gamma1.csv')

# seq1=nburn:ngibbs
# betas.estim=vec.betas[seq1,]
# alpha.estim=vec.alpha[seq1,]
# gamma.estim=vec.gamma1[seq1,]

#get regression parameters
alpha=rep(0,ncol(alpha.estim))
alpha[1]=-1; alpha[2]=0.5;alpha[3]=-0.5
alpha.true=alpha

betas=c(-0.5,0,0.4,0.5)
betas.true=betas
gamma1.true=0.23

#compare
plot(gamma.estim$x,type='l')
# plot(gamma.estim,type='l')
abline(h=gamma1.true,col='red')

#invasion process param
boxplot(alpha.estim,ylim=range(c(alpha.true,alpha.estim)))
points(1:ncol(alpha.estim),alpha.true,col='blue',pch=19)

#observation param
boxplot(betas.estim,ylim=range(c(betas.true,betas.estim)))
points(1:ncol(betas.estim),betas.true,col='blue',pch=19)

for (i in 1:ncol(betas.estim)) plot(betas.estim[,i],type='l')
cor(cbind(gamma.estim$x,betas.estim))