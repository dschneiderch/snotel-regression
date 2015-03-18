library('ProjectTemplate')
setwd('~/GoogleDrive/snotel-regression_project')
#setwd('~/Documents/snotel-regression_project')
library(doMC)
load.project()

# select recon version ----------------------------------------------------
recon.version='v3.2'

# run model -------------------
cost='r2'#cor, r2, mae, rmse
style='real-time'#real-time'#'real-time','reanalysis'
spatialblend='blend'#flag to do geostatistical blending or not (prediction stage only; always does in the CV stage). blending takes long time for the entire domain..
output='points'#points' #'surface' #just predict at snotel pixels #for 'points' spatialblend must also be 'blend'
covrange='300km'

registerDoMC(3)

fullpreds=read.table(paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/','fullpreds/xval/',style,'_snotel_xval_bestpreds_',cost,'.txt'),header=T)
#
snotellocs=snotellocs[snotellocs$Station_ID %in% fullpreds$Station_ID,]#this is incase we do another year such as fassnacht years. keep in mind that only years availble 2000-2012 will be included regardless if there wre others.
snotellocs.usgs=spTransform(snotellocs,CRS('+init=epsg:5070'))

## --- calc Global Moran on residuals
moran.df=ddply(fullpreds,.(yrdoy),calc_GMoran,snotellocs.usgs,recon.version,covrange,.parallel=3)
write.table(moran.df,paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/',style,'_moran_info_for_recondate_selection_',cost,'.txt'),sep='\t',row.names=F,quote=F)     


## ---- calc Local Moran of SWE obs vs SWE model and compare
lmoran.df=ddply(fullpreds,.(yrdoy),calc_localMoran,snotellocs.usgs,recon.version,covrange,.parallel=3)
lmoran.df$date=as.POSIXct(lmoran.df$date)

plttitle2=paste0('\nIDW p=0.5, binary neighborhood')#from calc_localMoran.R

metric='mae'
if(metric=='mae'){ 
	costfun<-function(y, yhat) mean(abs(yhat-y),na.rm=T) #mae
} else if (metric=='rmse') {
	costfun<-function(y,yhat) sqrt(mean((yhat-y)^2,na.rm=T)) #rmse
}


snotellm=ddply(lmoran.df,.(date),summarise,
	min=min(snotel.I,na.rm=T),
	max=max(snotel.I,na.rm=T),
	avg=mean(snotel.I,na.rm=T),
	med=median(snotel.I,na.rm=T),
	yr=strftime(unique(date),'%Y'))


lmoran.metric=ddply(lmoran.df,.(yrdoy),function(dF){
	date=as.POSIXct(unique(dF$date))
	yr=strftime(date,'%Y')
	phv=costfun(dF$snotel.I,dF$phv.I)
	phvrel=phv/mean(dF$snotel.I,na.rm=T)
	phvrcn=costfun(dF$snotel.I,dF$phvrcn.I)
	phvrcnrel=phvrcn/mean(dF$snotel.I,na.rm=T)
	phv.full=costfun(dF$snotel.I,dF$phv.full)
	phv.fullrel=phv.full/mean(dF$snotel.I,na.rm=T)
	phvrcn.full=costfun(dF$snotel.I,dF$phvrcn.full)
	phvrcn.fullrel=phvrcn.full/mean(dF$snotel.I,na.rm=T)
	return(
		data.frame(	date,yr,phv,phvrel,phvrcn,phvrcnrel,phv.full,phv.fullrel,phvrcn.full,phvrcn.fullrel)
			)
})

cols=colnames(lmoran.metric)[!grepl('rel',colnames(lmoran.metric))]
lmoran=melt(lmoran.metric[,cols],.(date,yr,yrdoy))
g=ggplot(lmoran)+
	geom_linerange(data=snotellm,aes(x=date,ymin=min,ymax=max))+
	geom_point(data=snotellm,aes(x=date,y=avg))+
	geom_point(aes(x=date,y=value, colour=variable))+
	labs(x='Date',y=paste0(toupper(metric),' [unitless]'),title=paste0(toupper(metric),' of Local Moran I\'s',plttitle2))+
	coord_cartesian(ylim=c(-1,1))+
	theme_minimal()+
    theme(strip.text=element_text(face='bold',size=14),
        panel.background=element_rect(colour='grey80'))+
	facet_wrap(~yr,scale='free_x')

ggsave(plot=g,filename=paste0('graphs/rswe_',recon.version,'/covrange',covrange,'/scatterplot_SWE-localMoranI-',toupper(metric),'_',cost,'_',style,'_byyear_blended&unblended.pdf'))

cols=colnames(lmoran.metric)[grep('rel',colnames(lmoran.metric))]
lmoran.rel=melt(lmoran.metric[,c('date','yr',cols)],.(date,yr))
grel=ggplot(lmoran.rel)+
	geom_linerange(data=snotellm,aes(x=date,ymin=min,ymax=max))+
	geom_point(data=snotellm,aes(x=date,y=avg))+
	geom_point(aes(x=date,y=value, colour=variable))+
	labs(x='Date',y=paste0('Relative ',toupper(metric),' [/100]'),title=paste0('Relative ',toupper(metric),' of Local Moran I\'s',plttitle2))+
	coord_cartesian(ylim=c(-20,20))+
	theme_minimal()+
    theme(strip.text=element_text(face='bold',size=14),
        panel.background=element_rect(colour='grey80'))+
	facet_wrap(~yr,scale='free')

ggsave(plot=grel,filename=paste0('scatterplot_SWE-localMoranI-rel',toupper(metric),'_',cost,'_',style,'_byyear_blended&unblended.pdf'))