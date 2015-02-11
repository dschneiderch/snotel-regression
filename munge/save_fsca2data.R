setwd('~/GoogleDrive/snotel-regression_project')
library(raster)
library(ncdf4)


# ## attempt 1 ----  uses writeRaster. could work well for multiband geotiff or uncompressed raster format. however, with netcdf, it doesn't allow you to choose t dimension values which are assigned the %Y%j (yeardoy) value from the filenames.
# dn0=dir('/Volumes/hydroData/WestUS_Data/UCO_FSCA/',pattern='20',full.names=T)
# dn0=dn0[seq(1,length(dn0),2)]
# dn0=dn0[strtrim(basename(dn0),4) %in% c(2007,2008,2012)]#use smb
# #dn0=dn0[strtrim(basename(dn0),4) %in% c(seq(2000,2006),2009,2010,2011)]#use afp
# 
# lapply(dn0,function(d){
# 	yr=strtrim(basename(d),4)
# 	print(yr)
# 	fn=file.path('data','selectdates','modscag',paste0('fsca',yr,'.nc'))
# 	fns=list.files(d,full.names=T,pattern='*.tif$')
# 	s=stack(fns,quick=T)
# #writeRaster(s,file.path('data','selectdates','modscag',paste0('fsca',yr,'.tif')),options=c("COMPRESS=LZW", "TILED=YES"))
# #writeRaster(s,file.path('data','selectdates','modscag',paste0('fsca',yr,'.grd')))
# 
# writeRaster(s,fn, varname='fsca', varunit='percent', longname='fractional snow covered area', xname='Long', yname='Lat', zname='day', zunit='yeardoy')
# n=nc_open(fn,write=T)
# ncatt_put(n,varid=0,'dates',names(s))
# nc_close(n)
# 	})



## attempt 2  --- this seems to work better. but still can't give multiple attribute values
dn0=dir('/Volumes/hydroData/WestUS_Data/UCO_FSCA/',pattern='20',full.names=T)
dn0=dn0[seq(1,length(dn0),2)]
dn0=dn0[strtrim(basename(dn0),4) %in% c(2007,2008,2012)]#use smb to mount hydroData
#dn0=dn0[strtrim(basename(dn0),4) %in% c(seq(2000,2006),2009,2010,2011)]#use afp

lapply(dn0,function(d){
	yr=strtrim(basename(d),4)
	print(yr)
     fn=file.path('data','selectdates','modscag',paste0('fsca',yr,'.nc'))
	fns=list.files(d,full.names=T,pattern='*.tif$')
	s=stack(fns,quick=T)
	dim1=ncdim_def('Long','degree',seq(-112.25,-104.125,0.00416666667))
	dim2=ncdim_def('Lat','degree',seq(33,43.75,0.00416666667))
	dim3=ncdim_def('time','yrdoy',unlim=T,vals=as.numeric(gsub('X','',names(s))))
	var=ncvar_def('fsca','percent',dim=list(dim1,dim2,dim3),missval=-99,longname='fractional snow covered area',prec='integer',compression=6)
		ncnew=nc_create(fn,var)
		ncvar_put(ncnew, var, getValues(s))
		ncatt_put(ncnew,0,'proj4string',projection(s))#'+proj=longlat +datum=NAD83')
		ncatt_put(ncnew,0,'dates',names(s))
		nc_close(ncnew)
})