library('ProjectTemplate')
load.project()


s=stack('data/recon_v3.1/recondata_2000_v3.1.nc')

nwtind=grep('NIWOT',snotellocs$Site_Name)
nwtloc=snotellocs[nwtind,]
nwt=extract(s,nwtloc)
ind=which(nwt<1 & nwt>0)
nwt2=nwt[ind[1]:(ind[length(ind)]+1)]

a=min(nwt2)
b=max(nwt2)
z=(nwt2-a)/(b-a)
hist(z)
fitdistr(z,'beta',start=list(shape1=6,shape2=.3))


## Calculate mean and variance of the data.
temp=as.numeric(nwt2/(1-.7))
a=min(temp)
b=max(temp)
z=(temp-a)/(b-a)
hist(z)
mu <- mean(z)
var <- var(z)

## Define function to estimate parameters of a beta distribution.
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}


## Apply function to your data.
bparams=estBetaParams(mu, var)
alpha=bparams$alpha
beta=bparams$beta

z2=z[z<1 & z>0]
fitdistr(z2,'beta',start=list(shape1=.68,shape2=.3))
plot(1:length(z2),z2)
z2
dbeta()
dbeta(1:length(z2),alpha, beta)
predict.beta=function(alpha,beta)


library(betareg)
trial=data.frame(y=nwt[1,]/(1-0.7),x=1:length(nwt))
trial2=subset(trial,y>0 & y < 1)
mdl=betareg(y~x,data=trial2)
plot(trial2$y,predict(mdl,type='response'))
abline()


sbeta=calc(scrop,function(x){
	dF=data.frame(y=x/(1-0.7),x=1:length(x))
	dF=subset(dF,y>0 & y<1)
	mdl=betareg(y~x,data=dF)
	tryCatch(predict(mdl),error=function(e) rep(NA,nrow(dF)))
})