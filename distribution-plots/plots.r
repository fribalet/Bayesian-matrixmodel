##############################################
## ANALYZE MOMENTS ###########################
##############################################
library(fields)
library(e1071)
library(viridis)
library(Weighted.Desc.Stat)

dat <- read.csv('d:/dropbox/working/bayesian_matrix_model/github/data/SeaFlow_PSD_hourlyCOUNT_v3.csv')
st  <- as.numeric(substr(as.character(dat[,1]),12,13))
s   <- as.matrix(dat[,7:ncol(dat)])
sv  <- as.numeric(substr(colnames(s),2,7)
sv  <- sv/max(sv)
image.plot(s)

plot(s[2,])

zin <- read.csv('d:/dropbox/working/bayesian_matrix_model/github/data/FSC_Array5_logtransform.csv')
z   <- as.matrix(zin[,2:ncol(zin)])
zt  <- zin[,1]
zv  <- as.numeric(substr(colnames(z),2,7))
zv  <- zv/max(zv)
plot(z[5,])

smom=zmom <- data.frame(mean=numeric(),sd=numeric(),skew=numeric())

iss <- 24:72
for(i in iss){
#for(i in 1:nrow(s)){
	meann <- w.mean(sv,s[i,])
	sdd   <- w.sd(sv,s[i,])
	skeww <- w.skewness(sv,s[i,])
	smom  <- rbind(smom,data.frame(mean=meann,sd=sdd,skew=skeww))
}
for(i in 1:nrow(z)){
	meann <- w.mean(zv,z[i,])
	sdd   <- w.sd(zv,z[i,])
	skeww <- w.skewness(zv,z[i,])
	zmom  <- rbind(zmom,data.frame(mean=meann,sd=sdd,skew=skeww))
}

par(mfrow=c(1,2))
cols <- viridis(30)
is <- 2:26
plot(s[1,]/sum(s[1,]),type='l',ylim=c(0,0.1),col=cols[1])
for(i in is) lines(s[i,]/sum(s[i,]),col=cols[i])
legend('topright',legend=st,col=cols[1:26],lty=1,cex=0.5,bty='n')

is <- 2:25
plot(z[1,]/sum(z[1,]),type='l',col=cols[1])
for(i in is) lines(z[i,]/sum(z[i,]),col=cols[i])
legend('topright',legend=zt,col=cols[c(1,is)],lty=1,cex=0.5)


par(mfrow=c(3,2),mar=c(2,2,2,2),oma=c(2,2,2,2))
plot(smom$mean,type='l',ylim=c(range(c(smom$mean,zmom$mean))),xaxt='n')
	axis(1,labels=st[iss],at=seq(1,length(st[iss])))
	plot(zmom$mean,type='l',ylim=c(range(c(smom$mean,zmom$mean))),xaxt='n')
	axis(1,labels=zt,at=seq(1,length(zt)))
plot(smom$sd,type='l',ylim=c(range(c(smom$sd,zmom$sd))),xaxt='n',xaxt='n')
	axis(1,labels=st[iss],at=seq(1,length(st[iss])))
	plot(zmom$sd,type='l',ylim=c(range(c(smom$sd,zmom$sd))),xaxt='n')
	axis(1,labels=zt,at=seq(1,length(zt)))
plot(smom$skew,type='l',ylim=)
	plot(zmom$skew,type='l')

par(mfrow=c(2,2))
cols <- viridis(ncol(z))
plot(z[,1]/sum(z[,1]),type='l',ylim=c(0,0.3),col=cols[1])
for(i in 2:ncol(z)) lines(z[,i]/sum(z[,i]),col=cols[i])
plot(z[,1],type='l',ylim=c(0,5000000),col=cols[1],xlim=c(0,30))
for(i in 2:ncol(z)) lines(z[,i],col=cols[i])
legend('topright',legend=round(zv[seq(1,199,length.out=20)],4),lty=1,col=cols[seq(1,199,length.out=20)],cex=0.5,bty='n')

cols <- viridis(ncol(s))
plot(s[24:72,1]/sum(s[24:72,1]),type='l',ylim=c(0,0.1),col=cols[1])
for(i in 2:ncol(s)) lines(s[24:72,i]/sum(s[24:72,i]),col=cols[i])
plot(s[24:72,1],type='l',ylim=c(0,13),col=cols[1],xlim=c(0,58))
for(i in 2:ncol(s)) lines(s[24:72,i],col=cols[i])
legend('topright',legend=round(sv[seq(1,48,length.out=20)],4),lty=1,col=cols[seq(1,48,length.out=20)],cex=0.5,bty='n')




