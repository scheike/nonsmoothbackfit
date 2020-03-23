
library(mets)
###library(prodlim)
data(TRACE)

source("my-functions-backfit.r")

boot.back <-function(At,Aa,res2,xx,xt,boot=100,sd.only=0,wild=0,bands=1)
{# {{{

xxx <- At$cox.prep
idjump <- xxx$id[xxx$jumps+1]+1
ttimes <- At$cumhaz[,1]
dAt <- diff(c(0,At$cumhaz[,2]))
At <- Cpred(At$cumhaz,xt)[,2]
###
atimes <- Aa$cumhaz[,1]
xxx <- Aa$cox.prep
idjumpa <- xxx$id[xxx$jumps+1]+1
dAa <- diff(c(0,Aa$cumhaz[,2]))
Aa <- Cpred(Aa$cumhaz,xx)[,2]

res3<-backfitEaEt(res2$E1,res2$E2,Aa,At,xt,xx,it=40)
At0 <- res3$At
Aa0 <- res3$Ba

wherea <- sindex.prodlim(c(0,atimes),xx)
wheret <- sindex.prodlim(c(0,ttimes),xt)

AtB <- BaB <- matrix(0,length(xx),boot)

for (j in 1:boot) {
g <- rnorm(nrow(ct))+wild
At <- cumsum(dAt*g[idjump])
Aa <- cumsum(dAa*g[idjumpa])
###
At <- c(0,At)[wheret]
Aa <- c(0,Aa)[wherea]
###
res3<-backfitEaEt(res2$E1,res2$E2,Aa,At,xt,xx,it=40)
###res3$At
###res3$Ba
AtB[,j] <- res3$At   ### -tail(res3$At,1)*(xt/mt)
BaB[,j] <- res3$Ba   ### +tail(res3$At,1)*(xx)
}

sd.At=apply(AtB,1,sd);sd.Ba=apply(BaB,1,sd) 

if (bands==1) {
ppA <- percen(apply(abs(AtB/sd.At),1,max),0.95)
bandAt <- cbind(At0-ppA*sd.At,At0+ppA*sd.At)
ppB <- percen(apply(abs(BaB/sd.Ba),1,max),0.95)
bandBa <- cbind(Aa0-ppB*sd.Ba,Aa0+ppB*sd.Ba)
confAt <- cbind(At0-1.96*sd.At,At0+1.96*sd.At)
confBa <- cbind(Aa0-1.96*sd.Ba,Aa0+1.96*sd.Ba)
} else {
	ppA <- ppB <- bandAt <- bandBa <- NULL
	confAt <- confBa <- NULL
}

if (sd.only==0) return(
		       list(atimes=xx,ttimes=xt,AtB=AtB,BaB=BaB,
			    confAt=confAt,confBa=confBa,
			    At=At0,Ba=Aa0,ppA=ppA,ppB=ppB,
			    bandAt=bandAt,bandBa=bandBa,sd.At=sd.At,sd.Ba=sd.Ba)) 
else 
return(list(atimes=xx,ttimes=xt,At=At0,Ba=Aa0,ppA=ppA,ppB=ppB,
			    confAt=confAt,confBa=confBa,
	    bandAt=bandAt,bandBa=bandBa,sd.At=sd.At,sd.Ba=sd.Ba, 
	    mean.At=apply(AtB,1,mean),mean.Ba=apply(BaB,1,mean))) 
}# }}}

##################################################
### Running the trace data 
##################################################

### setting up the data  
data(TRACE)
dsort(TRACE) <- ~id
set.seed(100)
## make sure we do not have ties 
TRACE$time <- TRACE$time+runif(nrow(TRACE))*0.01
TRACE$age0 <- TRACE$age
TRACE$age <- TRACE$age+TRACE$time
TRACE$time0 <- 0
## start age and start time=0
TRACE$age.null <- TRACE$age0
TRACE$time.null <- TRACE$time0
head(TRACE)
head(TRACE)
dcut(TRACE) <- ~age
TRACE$vff <- factor(TRACE$vf)

dtable(TRACE,~vf+chf+diabetes+sex,level=2)

## decide where to look at model here box 
t0 <- 0
mt <- 5
a0 <- 40
ma <- 90
## ct is data that is only using relevant box where t \in [0,5] a \in [40,72]
## model \beta(a) + \alpha(t) 
ct <- timescale.box(t0,mt,a0,ma,TRACE,gplot=1)
ct$statusd <- ct$status!=0

### fit marginal two models  on time and age
at <- aalen(Surv(time0,time,status!=0)~+1,data=ct,robust=0,start.time=t0,max.time=mt)
aa <- aalen(Surv(age0,age,status!=0)~+1,data=ct,robust=0,start.time=a0,max.time=ma)

t0=0; mt=5   # interval time direction
a0=40; ma=90;   # interval age direction
m <- 100
xx <- seq(a0,ma,length=m)
xt <- seq(0,mt,length=m)

## overall 100
At <- Cpred(at$cum,xt)[,2]
Aa <- Cpred(aa$cum,xx)[,2]
res2<-getEaEt2.new(ct,t0,mt,a0,ma,xt,xx)
res3<-backfitEaEt(res2$E1,res2$E2,Aa,At,xt,xx)
At<-res3$At
Ba<-res3$Ba
par(mfrow=c(2,2))
plot(xx,Ba,type="l")
plot(xt,At,type="l")


################################################################
############ Bootstrapping the SE's ############################
################################################################


pat <- phreg(Surv(time0,time,status!=0)~cluster(id),data=ct)
paa <- phreg(Surv(age0,age,status!=0)~cluster(id),data=ct)
sdtrace <- boot.back(pat,paa,res2,xx,xt,boot=500,sd.only=0,wild=0)

papir=0
if (papir==1) pdf("fig1-marg-backfit.pdf")
par(mfrow=c(1,2))
plotl(aa$cum,xlab="Age (years)\n (a)",ylab="Cumulative hazard")
lines(xx,Ba,col=2,lty=2)
title(main="Age time-scale")
legend("topleft",c("marginal","two-time-scale model"),col=1:2,lty=1:2)
mets:::plotConfRegion(xx,sdtrace$bandBa,col=2)
mets:::plotConfRegion(xx,sdtrace$confBa,polygon=FALSE,col=2,lty=3)
###
plotl(at$cum,xlab="Duration (years)\n (b)",ylab="Cumulative hazard")
lines(xt,At,col=2,lty=2)
abline(v=xt[ which(At==max(At))],lty=3)
title(main="Duration time-scale")
legend("topleft",c("marginal","two-time-scale model"),col=1:2,lty=1:2)
mets:::plotConfRegion(xt,sdtrace$bandAt,col=2)
mets:::plotConfRegion(xt,sdtrace$confAt,polygon=FALSE,col=2,lty=3)
if (papir==1) dev.off()

xt[ which(At==max(At))]*365.25


### survival predictions for different starting ages 
survP <- function(sdtrace,a1)
{# {{{
BaB50 <- sdtrace$BaB[sdtrace$atimes>=a1,]
###
tt <- sdtrace$ttimes
atimes <- a1+tt
wherea <- sindex.prodlim(c(0,sdtrace$atimes),atimes)
###
BaB50 <-Cpred(cbind(sdtrace$atimes,sdtrace$BaB),atimes)[,-1]
BaB50 <- t(t(BaB50)- BaB50[1,])
AtB <- 	sdtrace$AtB
cumHazB <- BaB50+AtB
sd.cumH=apply(cumHazB,1,sd);
###
Ba <- Cpred(cbind(sdtrace$atime,sdtrace$Ba),atimes)[,-1]
Ba <- Ba-Ba[1]
cumhaz <- sdtrace$At+Ba
###
ppA <- percen(apply(abs(cumHazB/sd.cumH),1,max),0.95)
bandAt <- cbind(cumhaz-ppA*sd.cumH,cumhaz+ppA*sd.cumH)

return(list(band=bandAt,atimes=atimes,
	    times=atimes-a1,a1=a1,cumhaz=cumhaz))
}# }}}

S60 <- survP(sdtrace,60)
S70 <- survP(sdtrace,70)
S80 <- survP(sdtrace,80)

## marginal age survival
survmargP <- function(aa,a1)
{# {{{

cum <- aa$cum[aa$cum[,1]>=a1,]
cumhaz <- cum[,2]
cumhaz <- cumhaz-cumhaz[1]
atimes=cum[,1]
return(list( atimes=atimes,
	    times=atimes-a1,a1=a1,cumhaz=cumhaz))
}# }}}

mS60 <- survmargP(aa,60)
mS70 <- survmargP(aa,70)
mS80 <- survmargP(aa,80)

papir=0
if (papir==1) pdf("fig2-survival-pred.pdf")
par(mfrow=c(1,1))
with(S60,
plot(times,exp(-cumhaz),type="l",ylim=c(0,1),xlab="time (years)",ylab="Survival probability"))
with(S60,mets:::plotConfRegion(times,exp(-band),col=1))
with(S70,lines(times,exp(-cumhaz),col=2,lwd=2,lty=1))
with(S70,mets:::plotConfRegion(times,exp(-band),col=2))
with(S80,lines(times,exp(-cumhaz),col=3,lwd=2,lty=1))
with(S80,mets:::plotConfRegion(times,exp(-band),col=3))
###
with(mS60,lines(times,exp(-cumhaz),lty=2,col=1))
with(mS70,lines(times,exp(-cumhaz),lty=2,col=2))
with(mS80,lines(times,exp(-cumhaz),lty=2,col=3))
legend("bottomleft",
   c("two-time scale","age only","duration time only","60, 70, 80"),lty=c(1:3,0))
if (papir==1) dev.off()

