#allows iteration of predict_surfaces by year to create yearly stacks and ncdf output
predict_wrapper <- function(rdtsdF,style,newdata,newdatalocs,spatialblend){
     
     registerDoMC(allocated.cores)
     print(paste('--',unique(rdtsdF$yr),'--'))
     locs=spTransform(snotellocs,CRS('+init=epsg:5070'))
     bpth=paste0('diagnostics/rswe_',recon.version,'/fullpreds')
     
     fullpreds=ddply(rdtsdF,.(yrdoy),predict_surfaces,locs,newdata,newdatalocs,spatialblend,bpth,.inform=F,.parallel=F)#,.paropts=list(.export=ls(), .packages=.packages(all=T)))

#     #output rasters of each model type.
#      rasterlist=dlply(fullpreds,.(yrdoy),function(dF){
#           staticr=raster(dF$phv.fullpred,CRS('+proj=longlat +datum=WGS'))
#           dynr=raster(dF$phvrcn.fullpred,CRS('+proj=longlat +datum=WGS'))
#           return(list(staticr,dynr))
#           })
     phvsurf=lapply(rasterlist,'[[',1)
     phvrcnsurf=lapply(rasterlist,'[[',2)
     
     phvstack=stack(phvsurf)
     names(phvstack)=strftime(rdtsdF$date,'%Y%m%d')
     phvrcnstack=stack(phvrcnsurf)
     names(phvrcnstack)=strftime(rdtsdF$date,'%Y%m%d')
#     
     yr=unique(rdtsdF$yr)
     pfn=file.path(bpth,'netcdf',paste0('fullpreds-phv_',yr,'_',spatialblend,'.ncdf'))
     prfn=file.path(bpth,'netcdf',paste0('fullpreds-phvrcn_',yr,'_',spatialblend,'.ncdf'))
#     
writeRaster(phvstack,filename=pfn,format='CDF',varname='swe',varunit='meters',xname='long',xunit='deg',yname='lat',yunit='deg',zname='time',zunit='day',overwrite=T)
writeRaster(phvrcnstack,filename=prfn,format='CDF',varname='swe',varunit='meters',xname='long',xunit='deg',yname='lat',yunit='deg',zname='time',zunit='day',overwrite=T)
#      GIdf=ddply(fullpreds,.(yrdoy),plot_localmoran,locs,recon.version)#
}