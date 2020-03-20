library(ncdf4)
library(viridis)
library(fields)
library(rstan)
options(mc.cores = parallel::detectCores())

############################################################################
## DATA ####################################################################
############################################################################
setwd('d:/dropbox/working/bayesian-matrixmodel/')

ndays            <- 4.0
limit_to_numdays <- 4.0
stride_t_obs     <- 10
data             <- list()
data$dt          <- 10

source('data_processing.r')

##############################################################################
## STAN MODELING #############################################################
##############################################################################

data$i_test = c(rep(0,round(dim(data$obs)[2]*(2/3))),rep(1,round(dim(data$obs)[2]*(1/3))))

mod_free    <- stan_model('d:/dropbox/working/bayesian-matrixmodel/matrixmodel_freedelta_CV.stan')
mod_orig    <- stan_model('d:/dropbox/working/bayesian-matrixmodel/matrixmodel_orig_CV.stan')
mod_sig     <- stan_model('d:/dropbox/working/bayesian-matrixmodel/matrixmodel_sigmoidaldelta_CV.stan')
mod_free_np <- stan_model('d:/dropbox/working/bayesian-matrixmodel/matrixmodel_freedelta_no_prior_CV.stan')

mcmc_free    <- sampling(mod_free,    data=data, open_progress=TRUE,chains=4)
mcmc_orig    <- sampling(mod_orig,    data=data, open_progress=TRUE,chains=4)
mcmc_sig     <- sampling(mod_sig,     data=data, open_progress=TRUE,chains=4)
mcmc_free_np <- sampling(mod_free_np, data=data, open_progress=TRUE,chains=4)

post_free    <- extract(mcmc_free)
post_orig    <- extract(mcmc_orig)
post_sig     <- extract(mcmc_sig)
post_free_np <- extract(mcmc_free_np) 

free_mean    <- apply(post_free$mod_obspos,c(2,3),mean)
orig_mean    <- apply(post_orig$mod_obspos,c(2,3),mean)
sig_mean     <- apply(post_sig$mod_obspos, c(2,3),mean)
free_np_mean <- apply(post_free_np$mod_obspos,c(2,3),mean) 

free_resid    <- free_mean - data$obs
orig_resid    <- orig_mean - data$obs
sig_resid     <- sig_mean - data$obs
free_np_resid <- free_np_mean - data$obs

xin <- seq(1,data$nt_obs*30,length.out=175)

par(mfrow=c(2,2),mar=c(3,4,2,3),oma=c(2,2,2,2))
image.plot(x=xin,y=v,z=t(free_resid),xlab='',ylab=''); box(); mtext('free with prior')
image.plot(x=xin,y=v,z=t(orig_resid),xlab='',ylab=''); box(); mtext('allometric')
image.plot(x=xin,y=v,z=t(sig_resid),xlab='',ylab=''); box(); mtext('sigmoidal')
image.plot(x=xin,y=v,t(free_np_resid),xlab='',ylab=''); box(); mtext('free no prior')
	mtext(side=1,outer=TRUE,'Time (minutes)')
	mtext(side=2,outer=TRUE,expression('Size ('*mu*'m)'))
	
par(mfrow=c(5,3),mar=c(2,2,1,2))
for(i in 15:1){
	acf(sig_res[i,],lag.max=157,main='',ylim=c(-1,1),ci=0,xaxt='n')
	axis(side=1,at=c(0,50,100,150),labels=c(0,50,100,150)*30)
	mtext(round(v[i],3),cex=0.7)
}

par(mfrow=c(5,3),mar=c(2,2,2,2))
for(i in 15:1) {
	plot(data$obs[i,],ylim=c(0,0.25),cex=0.3,pch=2)
	lines(sig_mean[i,])
	mtext(round(v[i],3),cex=0.7)
}
