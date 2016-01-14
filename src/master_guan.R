library('ProjectTemplate')
#setwd('~/GoogleDrive/snotel-regression_project')
#setwd('~/Documents/snotel-regression_project')
library(doMC)
load.project()

# select recon version ----------------------------------------------------
recon.version='v3.1'
ryrs=seq(2010,2012)
snotelrecon=get_snotelrecondata(recon.version,ryrs)

# generate swedata dataframe
swedata=generate_swedata(snotelrecon)
#subset snotel locations to those we have. snotelrec is loaded from cached files and should have only the stations we are using
swedata=swedata[swedata$Station_ID %in% unique(snotelrec$Station_ID),]
snotellocs=snotellocs[snotellocs$Station_ID %in% unique(swedata$Station_ID),]
snotelrecon=snotelrecon[snotelrecon$Station_ID %in% unique(swedata$Station_ID),]

# # setup parallel processing -----------
 getnodes <- function() {
        f <- Sys.getenv('PBS_NODEFILE')
       x <- if (nzchar(f)) readLines(f) else rep('localhost', 2)
        as.data.frame(table(x), stringsAsFactors=FALSE)
       }
 nodes <- getnodes()
 cl <- makeSOCKcluster(nodes$x)

# print(cl)
# print(nodes)

 registerDoSNOW(cl)

 setcores <- function(cl, nodes) {
      cores=nodes$Freq
    f <- function(cores) assign('allocated.cores', cores, envir=.GlobalEnv)
    clusterApply(cl, nodes$Freq, f)
  }
 setcores(cl, nodes)


## Select dates to model ----------
# Option A will only model dates for dates selected from modscag images.
# Option B will model selected dates from modscag for months prior to March and then 1st, 8th, 15th, 22nd of March, April, May
#st=as.POSIXct(paste0(ryrs,'-03-01'),tz='MST')
#aaply(as.integer(st),function(x){
#	yr=strftime(as.POSIXct(as.integer(x),origin='1970-01-01',tz='MST'),'%Y')
#	ed=as.POSIXct(paste0(yr,'-05-31'),tz='MST')
#	as.integer(seq(x,ed,3600*24))
#})

#dateselect=as.POSIXct(expand.grid(ryrs,dates2model(opt='B')
#ind=mdldata$date %in% dateselect
#doidata=mdldata[ind,]

## Do crossvalidation

#define domain for prediction -------
snotellocs.usgs=spTransform(snotellocs,CRS('+init=epsg:5070'))
newdata=as.data.frame(scale(ucophv))#ucophv is automatically loaded with project and  contains newdata for prediction to domain
newdatalocs=SpatialPoints(ucophv[,c('Long','Lat')])
proj4string(newdatalocs)='+proj=longlat +datum=WGS84'
newdatalocs.usgs=spTransform(newdatalocs,CRS('+init=epsg:5070'))
newdatalocs.agg=SpatialPoints(raster('data/gmted_combined_uco6km.tif'))
proj4string(newdatalocs.agg)='+proj=longlat +datum=WGS84'
newdatalocs.agg.usgs=spTransform(newdatalocs.agg,CRS('+init=epsg:5070'))

covrange='idp0.5'
print(ryrs)

l_ply(as.list(ryrs),.inform=F,.parallel=F,.paropts=list(.export=ls(), .packages=.packages(all=T)),
	.fun=function(yr){
print(yr)
recondates=data.frame(
		date=seq(as.POSIXct(paste0(yr,'-03-01'),tz='MST'),as.POSIXct(paste0(yr,'-08-31'),tz='MST'),by=3600*24),
		yr=yr)
	recondates$yrdoy=strftime(recondates$date,'%j')
	blend_recon(recondates,newdatalocs.usgs,newdatalocs.agg.usgs,recon.version,covrange)
})

stopCluster(cl)
quit('no')
