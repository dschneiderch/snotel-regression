library('ProjectTemplate')
setwd('~/GoogleDrive/snotel-regression_project')
#setwd('~/Documents/snotel-regression_project')
library(doMC)
reload.project()

# select recon version ----------------------------------------------------
recon.version='v3.2'
ryrs=seq(2000,2012)
snotelrecon=get_snotelrecondata(recon.version,ryrs)

# generate swedata dataframe
swedata=generate_swedata(snotelrecon)
#subset snotel locations to those we have.
snotellocs=snotellocs[snotellocs$Station_ID %in% unique(swedata$Station_ID),]

# get modelling data --------------------------------------------------------
mdldata=swedata[swedata$mth<6,]#if changing mth<6, need to change dates2model too
mdldata$recondate=mdldata$date

# select recon dates to evaluate model with
#recondata is used to iterate through recondates. 
#these dates can be different than olddata/mdldata if wanted but must include at least 1 date before each mdl date if used in real-time mode
recondata=subset(mdldata, (dy==1 | dy==15) ) 

## model options
cost='rmse'#cor, r2, mae, rmse
style='real-time'#real-time'#'real-time','reanalysis'
spatialblend='noblend'#flag to do geostatistical blending or not (prediction stage only; always does in the CV stage). blending takes long time for the entire domain..
output='surface'#points' #'surface' #just predict at snotel pixels #for 'points' spatialblend must also be 'blend'
covrange='300km'

##### To continue from previous run
which_recon_date=read.table(paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/',style,'_recondate_selection_',cost,'.txt'),sep='\t',header=T,stringsAsFactors=F)
which_recon_date$date=as.POSIXct(strptime(which_recon_date$date,'%Y-%m-%d',tz='MST'))
which_recon_date$phvrcn_recondate=as.POSIXct(strptime(which_recon_date$phvrcn_recondate,'%Y-%m-%d',tz='MST'))
which_recon_date$recon_costdate=as.POSIXct(strptime(which_recon_date$recon_costdate,'%Y-%m-%d',tz='MST'))
str(which_recon_date)
#####


#define domain for prediction -------
snotellocs.usgs=spTransform(snotellocs,CRS('+init=epsg:5070'))
newdata=as.data.frame(scale(ucophv))#ucophv is automatically loaded with project and  contains newdata for prediction to domain
newdatalocs=SpatialPoints(ucophv[,c('Long','Lat')])
proj4string(newdatalocs)='+proj=longlat +datum=WGS84'
newdatalocs.usgs=spTransform(newdatalocs,CRS('+init=epsg:5070'))
newdatalocs.agg=SpatialPoints(raster('data/gmted_combined_uco6km.tif'))
proj4string(newdatalocs.agg)='+proj=longlat +datum=WGS84'
newdatalocs.agg.usgs=spTransform(newdatalocs.agg,CRS('+init=epsg:5070'))
#
#geope=projectExtent(ucophv.stack[[1]],crs=projection(ucophv.stack))
#
#rskel=ucophv.stack[[1]]
#values(rskel)=F
#names(rskel)=NULL
#rskel=projectRaster(rskel,crs=CRS('+init=epsg:5070'))

# which_recon_date=which_recon_date[c(1,2,32,33),]
registerDoMC()
# save predicted surfaces to netcdf by year.
d_ply(which_recon_date,.(yr),predict_wrapper,style,newdata,newdatalocs.usgs,newdatalocs.agg.usgs,spatialblend,.inform=T,.parallel=F,.paropts=list(.export=ls(), .packages=.packages(all=T)))
#d_ply(wrd,.(yr),predict_wrapper,style,newdata,newdatalocs.usgs,newdatalocs.agg.usgs,spatialblend,.inform=T,.parallel=F,.paropts=list(.export=ls(), .packages=.packages(all=T)))

#quit(save='no')