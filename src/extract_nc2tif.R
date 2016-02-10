## regression predictions

recon.version='v3.1'
covrange='idp1'
cost='r2'
product='phvrcn'
fscaMatch='wofsca'
dateflag='B'
resid=''
style='real-time'
config=''#bic-v2.5-removed_NWbdiff'#'bic-v2.2-nosnoteltransform'
# yr=2010
# mth=4
dte=as.POSIXct('1900-03-01',tz='MST')
for(scalesnotel in c('scale')){
for(product in c('phvrcn')){#
for(fscaMatch in c('wofsca')){
for(resid in c('')){
for(yr in seq(2001,2012)){
	# for(mth in strftime(dte,'%m')){
		# month=strftime(dte,'%B')
		# for(dy in strftime(dte,'%d')){
			predfile=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',scalesnotel,'/',dateflag,'/fullpreds/',cost,'/netcdf/',style,'/',config,'/',resid,'preds-',product,'_',yr,'_blend-',fscaMatch,'.nc')
			s=brick(predfile)
			#band=grep(paste0(yr,mth,dy),names(s))
			# writeRaster(x=s,filename='/Volumes/fromDominik/forLeanne/20160203/preds_phvrcn.tif',bylayer=TRUE,NAflag=-99,suffix='names',overwrite=TRUE)
			# file2eval=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',snotelscale,'/',dateflag,'/fullpreds/',cost,'/netcdf/',style,'/',config,'/',resid,'preds-',product,'_',dy,month,yr,'-',fscaMatch,'.tif')	
			# cmd1=paste0('gdal_translate -b ',band, ' -a_srs +init=epsg:4326 -a_nodata -99 ',predfile,' ',file2eval)
			# system(cmd1)
			s[s>200]=NA
			smax=max(s,na.rm=TRUE)
			writeRaster(x=smax,filename=paste0('/Volumes/fromDominik/forLeanne/20160203/preds_phvrcn_max-',yr,'.tif'),NAflag=-99,overwrite=TRUE)
		}
	}
}
}
}
# }}


# ## recon fsca in data/recon_v3.1
# dte=as.POSIXct('1900-06-01',tz='MST')
# for(yr in seq(2000,2012)){
# 	for(mth in strftime(dte,'%m')){
# 		month=strftime(dte,'%B')
# 		for(dy in strftime(dte,'%d')){
# 			predfile=paste0('data/recon_',recon.version,'/reconfscadata_',yr,'_',recon.version,'.nc')
# 			s=stack(predfile)
# 			band=grep(paste0(yr,mth,dy),names(s))
# 			file2eval=paste0('data/recon_',recon.version,'/reconfscadata_',dy,month,yr,'.tif')
# 			cmd1=paste0('gdal_translate -b ',band, ' -a_srs +init=epsg:4326 -a_nodata -99 ',predfile,' ',file2eval)
# 			system(cmd1)
# 		}
# 	}
# }