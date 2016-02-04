# setwd('/Users/dosc3612/Documents/snotel-regression_project')
library(ProjectTemplate)
reload.project(list(cache_loading=F))

library(raster)
library(dplyr)
library(tidyr)
library(ggplot2)
library(doMC)
library(ncdf4)
library(RColorBrewer)


NCARv='1'
recon.version='v3.1'
covrange='idp1'
cost='r2'
residblend=''
scalesnotel='scale'
fscaMatch='wofsca'
dateflag='B'
style='real-time'#nalysis'
gridding='_modisgrid'#'_snodasgrid' or '_modisgrid'
postscaled=''#'' or '_postscaled'
predictor='rcn' #'fsca' or 'rcn'

nccoordid=nc_open('data/NCAR_WRF/original_runs/geographic_data.nc')
lat=ncvar_get(nccoordid,'XLAT')
long=ncvar_get(nccoordid,'XLONG')
# coords=data_frame(long=as.numeric(long[,,1]),lat=as.numeric(lat[,,1]))

wrfcorner=c(min(long),max(long),min(lat),max(lat))
wrfpoly=as(extent(wrfcorner),'SpatialPolygons')
proj4string(wrfpoly)='+proj=longlat +ellps=sphere +a=6370000 +b=6370000'
wrfpolylcc=spTransform(wrfpoly,CRS('+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39.0000038146973 +lon_0=-107 +ellps=sphere +a=6370000 +b=6370000 +units=m +no_defs'))
plot(wrfpoly)
plot(wrfpolylcc)
(xmin(extent(wrfpolylcc))-xmax(extent(wrfpolylcc)))/318
xmin(extent(wrfpolylcc)) + 4000*317
xmax(extent(wrfpolylcc))

extent(wrflcc)
wrfcds=data.frame(long=as.numeric(long[,,1]),lat=as.numeric(lat[,,1]))
wrfgrid=SpatialPoints(wrfcds)
proj4string(wrfgrid)='+proj=longlat +ellps=sphere +a=6370000 +b=6370000'
wrflcc=spTransform(wrfgrid,CRS('+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39.0000038146973 +lon_0=-107 +ellps=sphere +a=6370000 +b=6370000 +units=m +no_defs'))
wrflcc$swe=as.numeric(dt2[,,417])
SpatialGridDataFrame(wrflcc,dt2[,,417])

tmp=data.frame(x=c(min(long),min(long),max(long),max(long)),y=c(min(lat),max(lat),max(lat),min(lat)))
p4s='+proj=longlat +ellps=sphere +a=6370000 +b=6370000'
p4crs=CRS('+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39.0000038146973 +lon_0=-107 +ellps=sphere +a=6370000 +b=6370000 +units=m +no_defs')
tmppts=SpatialPoints(tmp)
proj4string(tmppts)=p4s
tmppts.lcc=spTransform(tmppts,p4crs)
coords=coordinates(tmppts.lcc)

## small example
tmp=data.frame(x=c(min(long),min(long),min(long)+5,max(long),max(long)),y=c(min(lat),max(lat),min(lat)+5,max(lat),min(lat)))
extent(tmp)
p4s='+proj=longlat +ellps=sphere +a=6370000 +b=6370000'
p4crs=CRS('+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39.0000038146973 +lon_0=-107 +ellps=sphere +a=6370000 +b=6370000 +units=m +no_defs')
tmppoly=as(extent(tmp),'SpatialPolygons')
tmppts=SpatialPoints(tmp)
proj4string(tmppoly)=p4s
proj4string(tmppts)=p4s
plot(tmppoly)
plot(tmppts,add=T,col='red')
tmppts.lcc=spTransform(tmppts,p4crs)
tmppoly.lcc=spTransform(tmppoly,p4crs)
plot(tmppoly.lcc)
plot(tmppts.lcc,add=T,col='red')
###


r=raster(matrix(as.numeric(lat[,,1]),byrow=T,nrow=263))
extent(r)=c(min(long),max(long),min(lat),max(lat))
r
plot(r)

ncswe=nc_open('data/NCAR_WRF/original_runs/SWE_daily.nc')
ncswe
data_temp=ncvar_get(ncswe,'SNOW')/917#kg/m^3 #convert to meters from kg/m^2
dt2=aperm(data_temp,c(2,1,3))
dt2=dt2[nrow(dt2):1,,]
rb=brick(dt2)
extent(rb)=c(min(long),max(long),min(lat),max(lat))
crs(rb)='+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39.0000038146973 +lon_0=-107 +a=6370000 +b=6370000 +units=m +no_defs'


#
startdate=as.POSIXct('2000-10-01',tz='MST')
timevec=as.numeric(strftime(seq(startdate,startdate+(3014-1)*3600*24,by=3600*24),format = '%Y%m%d'))

names(rb)=timevec

setZ(rb,timevec,'time')

#
fn='~/Downloads/swetest21.tif'
writeRaster(rb,fn,overwrite=TRUE,NAflag=-99,options=c("BLOCKXSIZE=317","BLOCKYSIZE=263","compress=LZW"))
,varname='swe',varunit='meters',zname='time',zunit='%Y%m%d')

brick(fn)
rb


ncid=nc_open(fn,write=TRUE)
ncatt_put(ncid,0,attname='mask_flagging',attval='values >200')
ncatt_put(ncid,0,attname='recon.version',attval = 'v3.1')
ncatt_put(ncid,0,attname='covrange',attval = 'idp1')
ncatt_put(ncid,0,attname='cost',attval = 'r2')
ncatt_put(ncid,0,attname='residblend',attval = 'none')
ncatt_put(ncid,0,attname='scalesnotel',attval = 'scale')
ncatt_put(ncid,0,attname='fscaMatch',attval = 'wofsca')
ncatt_put(ncid,0,attname='dateflag',attval = 'B')
ncatt_put(ncid,0,attname='style',attval = 'real-time')


startdate=as.POSIXct('2000-10-01',tz='MST')
timevec=as.numeric(strftime(seq(startdate,startdate+(3014-1)*3600*24,by=3600*24),tz='MST'))
wantdate=as.POSIXct('2002-01-01',tz='MST')
tind=as.numeric(wantdate-startdate)+1
# image(data_temp[,,tind])
datamat=as.matrix(data_temp[,,tind], by.row = F)
datamat=aperm(datamat,c(2,1))
datamat=datamat[nrow(datamat):1,]
rr=raster(datamat)
extent(rr)=c(min(long),max(long),min(lat),max(lat))
crs(rr)='+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39.0000038146973 +lon_0=-107 +a=6370000 +b=6370000 +units=m +no_defs'
rr

latdim=ncdim_def('lat',units='degrees_north',vals=unique(coordinates(rb)[,2]),longname='latitude')
longdim=ncdim_def('long',units='degrees_east',vals=rev(unique(coordinates(rb)[,1])),longname='longitude')
timedim=ncdim_def('time',units='%Y%m%d',vals=timevec)
ncswevar=ncvar_def(name='swe',dim=list(longdim,latdim,timedim),chunksizes = c(317,263,1),units='meters',missval =-99,longname='snow water equivalent')
ncprojvar=ncvar_def(name='lambert_conformal_conic',units='',dim=list(),prec='char')
ncnew=nc_create('~/Downloads/swetest-ncdf4.nc',var=list(ncprojvar))
ncvar_put(ncnew,ncswevar,getValues(rb))
# ncatt_put(ncnew,'lambert_conformal_conic','grid_mapping=lambert_conformal_conic')
ncatt_put(ncnew,varid='lambert_conformal_conic',attname='latitude_of_projection_origin',attval=39.0000038146973)


lambert_conformal_conic#false_easting=0
lambert_conformal_conic#false_northing=0
lambert_conformal_conic#GeoTransform=-114.8643951416016 0.04961758682780461 0 43.73397064208984 0 -0.0368182831390729
lambert_conformal_conic#grid_mapping_name=lambert_conformal_conic
lambert_conformal_conic#inverse_flattening=0
lambert_conformal_conic#latitude_of_projection_origin=39.0000038146973
lambert_conformal_conic#longitude_of_central_meridian=-107
lambert_conformal_conic#longitude_of_prime_meridian=0
lambert_conformal_conic#semi_major_axis=6370000

ncatt_put(ncnew,0,'proj4string','+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39.0000038146973 +lon_0=-107 +a=6370000 +b=6370000 +units=m +no_defs')
nc_close(ncnew)

left=-105.75
right=-105.25
top=40.125
bottom=39.875
bb=c(left,right,top,bottom)
box=as(extent(bb), "SpatialPolygons")
proj4string(box)='+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39.0000038146973 +lon_0=-107 +a=6370000 +b=6370000 +units=m +no_defs'

nc1=nc_open('~/Downloads/swetest-ncdf4.nc',write=TRUE)
ncatt_put(nc1,0,'spatial_ref','PROJCS["unnamed",GEOGCS["unnamed ellipse",DATUM["unknown",SPHEROID["unnamed",6370000,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Lambert_Conformal_Conic_2SP"],PARAMETER["standard_parallel_1",33],PARAMETER["standard_parallel_2",45],PARAMETER["latitude_of_origin",39.0000038146973],PARAMETER["central_meridian",-107],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]]]')
ncatt_put(nc1,0,'proj4string','+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39.0000038146973 +lon_0=-107 +a=6370000 +b=6370000 +units=m +no_defs')
ncatt_put(nc1,0,'grid_mapping','lambert_conic_conformal')
nc_close(nc1)
          


brick('~/Downloads/swetest-ncdf4.nc')      
new=nc_open('~/Downloads/swetest-ncdf4.nc')
new
plot(new[[tind]])

for(tind in 1:dim(data_temp)[3]){
  datamat=as.matrix(data_temp[,,tind],by.row=F)
  datamat=aperm(datamat,c(2,1))
  datamat=datamat[nrow(datamat):1,]
  rr=raster(datamat)
  ncvar_put(ncnew,ncswevar,getValues(rr),start=c(1,1,tind),count=c(263,317,1))

}

# test=data_frame(coords,swe=datavec)
#
# for(postscaled in c('')){
#      for(style in c('real-time')){
#           for(dateflag in c('B')){
#                for(residblend in c('')){
#                     # scalesnotel='noscale'
#                     for(scalesnotel in c('scale')){
#                          # fscaMatch='wofsca'
#                          if(scalesnotel=='scale' && postscaled=='_postscaled') next
#
#                          for(fscaMatch in c('wofsca')){
#                               # if(cost=='rmse' && fscaMatch=='fsca'){
#                               #   next#rmse/fsca and r2/fsca will have the same predictions
#                               # }
#
#                               config=''#'bic-v2.5-removed_NWbdiff'
#                               basepath=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',scalesnotel,'/',dateflag,'/fullpreds/',cost,'/netcdf/',style,'/',config,'/')
#
#
#
#                          }}}}}}
