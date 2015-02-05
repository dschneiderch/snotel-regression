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

 # setup parallel proecessing -----------
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
allocated.cores=2

# get modelling data --------------------------------------------------------
mdldata=swedata[swedata$mth<6,]#if changing mth<6, need to change dates2model too
mdldata$recondate=mdldata$date

# select dates to model ----------
# Option A will only model dates for dates selected from modscag images.
# Option B will model selected dates from modscag for months prior to March and then 1st and 15th of March, April, May
dateselect=dates2model(opt='A')
ind=mdldata$date %in% dateselect
doidata=mdldata[ind,]

# select recon dates to evaluate model with
recondata=subset(mdldata, (dy==1 | dy==15) & yr!=2012 ) #recondata is used to iterate through recondates. these dates can be different than olddata/mdldata if wanted but must include at least 1 date before the mdl date if used in real-time mode

# run model -------------------
cost='cor'
style='real-time'#'real-time'
spatialblend='blend'#flag to do geostatistical blending or not (prediction stage only; always does in the CV stage). blending takes long time for the entire domain..
#find which recondate gives best estimate
which_recon_date=ddply(doidata,.(yrdoy),doDOYfit,cost,style,.parallel=F, .drop=F,.inform=T)
# ,.paropts=list(.export=c('rundata','doXVAL'),
#                .packages=c('ipred','MASS')
# )
# )

write.table(which_recon_date,paste0('diagnostics/rswe_',recon.version,'/',style,'_recondate_selection_',cost,'.txt'),sep='\t',row.names=F,quote=F)     

##### To continue from previous run
which_recon_date=read.table(paste0('diagnostics/rswe_',recon.version,'/',style,'_recondate_selection_',cost,'.txt'),sep='\t',header=T,stringsAsFactors=F)
which_recon_date$date=as.POSIXct(strptime(which_recon_date$date,'%Y-%m-%d',tz='MST'))
which_recon_date$phvrcn_recondate=as.POSIXct(strptime(which_recon_date$phvrcn_recondate,'%Y-%m-%d',tz='MST'))
which_recon_date$reconcordate=as.POSIXct(strptime(which_recon_date$reconcordate,'%Y-%m-%d',tz='MST'))
str(which_recon_date)
#####


#define domain for prediction -------
newdata=as.data.frame(scale(ucophv))#ucophv is automatically loaded with project and  contains newdata for prediction to domain
newdatalocs=SpatialPoints(ucophv[,c('Long','Lat')])
proj4string(newdatalocs)='+proj=longlat +datum=NAD83'
newdatalocs.usgs=spTransform(newdatalocs,CRS('+init=epsg:5070'))
#
#geope=projectExtent(ucophv.stack[[1]],crs=projection(ucophv.stack))
#
#rskel=ucophv.stack[[1]]
#values(rskel)=F
#names(rskel)=NULL
#rskel=projectRaster(rskel,crs=CRS('+init=epsg:5070'))

#output dF of GLobal Moran I and objective functions for best model estimates for each yrdoy (by year). save predicted surfaces to netcdf by year.
d_ply(which_recon_date,.(yr),predict_wrapper,style,newdata,newdatalocs.usgs,spatialblend,.inform=F,.parallel=F,.paropts=list(.export=ls(), .packages=.packages(all=T)))

quit(save='no')
