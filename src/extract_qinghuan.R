library(raster)
library(ncdf4)

left=-105.75
right=-105.25
top=40.125
bottom=39.875
bb=c(left,right,top,bottom)
box=as(extent(bb), "SpatialPolygons")
proj4string(box)='+proj=longlat +datum=WGS84'

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

config=''#
basepath=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',scalesnotel,'/',dateflag,'/fullpreds/',cost,'/netcdf/',style,'/',config,'/')

for(yr in 2001:2012){
     phv=brick(paste0(basepath,residblend,'preds-phvrcn_',yr,'_blend-',fscaMatch,'.nc'))
     rr=crop(phv,box)
     fn=paste0('~/Downloads/phvrcn-',yr,'.nc')
     writeRaster(rr,fn,overwrite=TRUE)
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
     nc_close(ncid)
}


