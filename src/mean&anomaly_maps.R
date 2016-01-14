library('ProjectTemplate')
setwd('~/Documents/snotel-regression_project')
load.project()
# library(raster) 
# library(rgdal)
# library(plyr)
# library(reshape2)
# library(gtable)
# library(gridExtra)

recon.version='v3.1'
covrange='idp1'
cost='r2'
product='phvrcn'
fscaMatch='fsca'
resid='full'
config='bic-v2.5-removed_NWbdiff'#'bic-v2.2-nosnoteltransform'
yr=2010
mth=4
dte=as.POSIXct('1900-04-08',tz='MST')
for(snotelscale in c('scale')){
for(product in c('phvrcn')){#
for(fscaMatch in c('fsca')){
for(resid in c('full')){
for(yr in seq(2001,2012)){
	for(mth in strftime(dte,'%m')){
		month=strftime(dte,'%B')
		for(dy in strftime(dte,'%d')){
			basepath=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',snotelscale,'/fullpreds/',cost,'/netcdf/',config,'/')
			fns=list.files(path=basepath,pattern=glob2rx(paste0(resid,'preds-',product,'_',dy,month,'*',fscaMatch,'.tif')),full.names=T)
			b=stack(fns)
			projection(b)='+proj=longlat +datum=WGS84'
			b[b==253]=NA
			bavg=mean(b,na.rm=T)
			bavg.fn=paste0(basepath,resid,'preds-',product,'_',dy,month,'_mean-',fscaMatch,'.tif')
			writeRaster(bavg,filename=bavg.fn)
			#
			banom=b-bavg
			banom.fn=paste0(basepath,resid,'preds-',product,'_',dy,month,'_anom-',fscaMatch,'.tif')
			writeRaster(banom,filename=banom.fn,bylayer=T)
			#
			banom.percent=banom/bavg
			banom.per.fn=paste0(basepath,resid,'preds-',product,'_',dy,month,'_anompercent-',fscaMatch,'.tif')
			writeRaster(banom.percent,filename=banom.per.fn,bylayer=T)
		}
	}
}
}
}
}}



