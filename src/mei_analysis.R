setwd('~/Documents/snotel-regression_project')
library(ProjectTemplate)
load.project()


fn='~/Documents/snotel-regression_project/diagnostics/rswe_v3.1/covrangeidp1/snotelscale/fullpreds/r2/netcdf/fullpreds-phvrcn_2001_blend-fsca.nc'
system(paste('gdalinfo',fn))

s=stack(fn)

z=s[50000]
x=z[z<200]
calc(s,function(x){
	x=x[x<200]
	qx=quantile(x)
	qf=rep(NA, length(x))
	for(p in c('75%','50%','25%')){
		qi=grep(p,attr(qx,'names'))	
		qf[x<=qx[qi]]=p			
	}
qf	
F=ecdf(x)
F(x)
quantile(F)
as.numeric(x)
})

mei=read.table('data/climate/wolter_mei.txt',sep='\t',header=T,fill=T)
mei=subset(mei,YEAR>2000)

str(snotelrec)
mth_avg_swe=dcast(snotelrec,wy+yr+mth~.,value.var='swe',mean)
colnames(mth_avg_swe)=c('wy','yr','mth','swe')
mei.m=melt(mei,'YEAR')
mei.m$mth=ifelse(grepl('DECJAN',mei.m$variable),1,
	ifelse(grepl('JANFEB',mei.m$variable),2,
	ifelse(grepl('FEBMAR',mei.m$variable),3,
	ifelse(grepl('MARAPR',mei.m$variable),4,
	ifelse(grepl('APRMAY',mei.m$variable),5,
	ifelse(grepl('MAYJUN',mei.m$variable),6,
	ifelse(grepl('JUNJUL',mei.m$variable),7,
	ifelse(grepl('JULAUG',mei.m$variable),8,
	ifelse(grepl('AUGSEP',mei.m$variable),9,
	ifelse(grepl('SEPOCT',mei.m$variable),10,
	ifelse(grepl('OCTNOV',mei.m$variable),11,
	ifelse(grepl('NOVDEC',mei.m$variable),12,NA))))))))))))
colnames(mei.m)=c('yr','mnths','mei','mth')

mei_swe=merge(mth_avg_swe,mei.m,by=c('yr','mth'))
mei_swe2=ddply(mei_swe,.(wy),function(dF){
	dF$swe_sc=scale(dF$swe,scale=F)
	return(dF)
})
with(mei_swe2,cor(swe_sc,mei))

pdo=read.table('data/climate/pdo.txt',header=F,sep='\t')
colnames(pdo)=c('yr','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
pdo.m=melt(pdo,'yr',variable.name='mnth',value.name='pdo')
pdo.m$mth=as.numeric(pdo.m$mnth)

master=merge(mei_swe2,pdo.m,c('yr','mth'))
master.sub=subset(master,mth<8 | mth>10)
master.sub$swe[master.sub$swe==0]=0.001

mdl=glm(swe~mei+pdo,data=master.sub,family=gaussian)

          myglm= function(formula,data) {
               glm(formula,data,family=gaussian())
          }
          mypredict.glm <- function(object, newdata){
               predict(object, newdata,type='response')
          }

hist(master.sub$mei)
hist(master.sub$pdo)
hist(master.sub$pdo*master.sub$mei)
mean(master.sub$swe)
errorest(formula=swe~mei+pdo,data=master.sub, model=myglm, predict=mypredict.glm, estimator='cv',est.para=control.errorest(k=nrow(master.sub), predictions=T))

plot(master.sub$swe,predict(mdl,type='response'))
plot(master.sub$swe,resid(mdl,type='response'))
