library('ProjectTemplate')
#setwd('~/GoogleDrive/snotel-regression_project')
library(doMC)
load.project()

# select recon version ----------------------------------------------------
recon.version='v3.2'
ryrs=seq(2000,2011)
snotelrecon=get_snotelrecondata(recon.version,ryrs)

# generate swedata dataframe
swedata=generate_swedata(snotelrecon)

# setup parallel processing -----------
getnodes <- function() {
       f <- Sys.getenv('PBS_NODEFILE')
       x <- if (nzchar(f)) readLines(f) else rep('localhost', 2)
       as.data.frame(table(x), stringsAsFactors=FALSE)
      }
nodes <- getnodes()
cl <- makeSOCKcluster(nodes$x)

print(cl)
print(nodes)

registerDoSNOW(cl)

setcores <- function(cl, nodes) {
     cores=nodes$Freq
   f <- function(cores) assign('allocated.cores', cores, envir=.GlobalEnv)
   clusterApply(cl, nodes$Freq, f)
 }
setcores(cl, nodes)


# get modelling data --------------------------------------------------------
mdldata=swedata[swedata$mth<6,]#if changing mth<6, need to change dates2model too
mdldata$recondate=mdldata$date

# select dates to model ----------
# Option A will only model dates for dates selected from modscag images.
# Option B will model selected dates from modscag for months prior to March and then 1st, 8th, 15th, 22nd of March, April, May
dateselect=dates2model(opt='B')
ind=mdldata$date %in% dateselect
doidata=mdldata[ind,]

# select recon dates to evaluate model with
#recondata is used to iterate through recondates. 
#these dates can be different than olddata/mdldata if wanted but must include at least 1 date before each mdl date if used in real-time mode
recondata=subset(mdldata, (dy==1 | dy==15) & yr!=2012 ) 

# run model -------------------
cost='r2'#cor, r2, mae, rmse
style='real-time'#'real-time','reanalysis'
spatialblend='blend'#flag to do geostatistical blending or not (prediction stage only; always does in the CV stage). blending takes long time for the entire domain..
output='points' #'surface' #just predict at snotel pixels #for 'points' spatialblend must also be 'blend'

if(style=='real-time'){
     doidata=subset(doidata,yr>2000)
}
#  dte=strptime('20010524','%Y%m%d','MST')
#  ldte=strptime('20020302','%Y%m%d','MST')
#  doidata=subset(doidata,date>=dte)
#  doidata=subset(doidata,date<ldte)
#find which recondate gives best estimate.
#output dF of GLobal Moran I and objective functions for best model estimates for each yrdoy
# doylist=dlply(doidata,.(yrdoy),doDOYfit,cost,style,.parallel=F, .drop=F,.inform=F)
# cleandF=function(dF){
#      mutate(dF,
#             yrdoy=attr(doylist,'split_labels')$yrdoy,
#             date=as.POSIXct(strptime(yrdoy,'%Y%j','MST')),
#             yr=strftime(date,'%Y'))
# }
# which_recon_date=cleandF(ldply(doylist,'[[',1))
# moran.df=cleandF(ldply(doylist,'[[',2))


# write.table(which_recon_date,paste0('diagnostics/rswe_',recon.version,'/',style,'_recondate_selection_',cost,'.txt'),sep='\t',row.names=F,quote=F)     
# write.table(moran.df,paste0('diagnostics/rswe_',recon.version,'/',style,'_moran_info_for_recondate_selection_',cost,'.txt'),sep='\t',row.names=F,quote=F)     


##### To continue from previous run
which_recon_date=read.table(paste0('diagnostics/rswe_',recon.version,'/',style,'_recondate_selection_',cost,'.txt'),sep='\t',header=T,stringsAsFactors=F)
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
d_ply(which_recon_date,.(yr),predict_wrapper,style,newdata,newdatalocs.usgs,spatialblend,.inform=T,.parallel=F,.paropts=list(.export=ls(), .packages=.packages(all=T)))

#quit(save='no')
