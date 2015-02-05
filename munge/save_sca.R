setwd('~/GoogleDrive/snotel-regression_project')

library(raster)
# dn1=dir('/Volumes/hydroData-1/WestUS_Data/UCO_FSCA/',pattern='20',full.names=T)
# dn1=dn1[seq(1,length(dn1),2)]
# dn1=dn1[strtrim(basename(dn1),4) %in% 2006]

dn0=dir('/Volumes/hydroData/WestUS_Data/UCO_FSCA/',pattern='20',full.names=T)
dn0=dn0[seq(1,length(dn0),2)]
dn0=dn0[strtrim(basename(dn0),4) %in% seq(2007,2012)]


# lapply(dn1,function(d){
# 	yr=strtrim(basename(d),4)
# print(yr)
# 	fn=list.files(d,full.names=T,pattern='*.tif$')
# #	print(fn)
# s=stack(fn,quick=T)
# writeRaster(s,file.path('data','selectdates','modscag',paste0('fsca',yr)))
# 	})
# 	
	
lapply(dn0,function(d){
	yr=strtrim(basename(d),4)
print(yr)
	fn=list.files(d,full.names=T,pattern='*.tif$')
#	print(fn)
s=stack(fn,quick=T)
writeRaster(s,file.path('data','selectdates','modscag',paste0('fsca',yr)))
	})
	
	
	