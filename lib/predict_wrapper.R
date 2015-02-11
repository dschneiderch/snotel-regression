#allows iteration of predict_surfaces by year to create yearly stacks and ncdf output
predict_wrapper <- function(rdtsdF,style,newdata,newdatalocs,spatialblend){
     
     registerDoMC(allocated.cores)
     print(paste('--',unique(rdtsdF$yr),'--'))
     bpth=paste0('diagnostics/rswe_',recon.version,'/fullpreds')
     
     fullpreds=ddply(rdtsdF,.(yrdoy),predict_surfaces,snotellocs.usgs,newdata,newdatalocs,spatialblend,bpth,.inform=F,.parallel=F)#,.paropts=list(.export=ls(), .packages=.packages(all=T)))     
     
     #output rasters of each model type.
     dim1=ncdim_def('Long','degree',seq(-112.25,-104.125,0.00416666667))
     dim2=ncdim_def('Lat','degree',seq(33,43.75,0.00416666667))
     dim3=ncdim_def('time','yrdoy',unlim=T,vals=unique(fullpreds$yrdoy))
     var=ncvar_def('swe','meters',dim=list(dim1,dim2,dim3),missval=-99,longname='snow water equivalent',compression=6)
     yr=unique(rdtsdF$yr)
     phv.fn=file.path(bpth,'netcdf',paste0('fullpreds-phv_',yr,'_',spatialblend,'.nc'))
     phvrcn.fn=file.path(bpth,'netcdf',paste0('fullpreds-phvrcn_',yr,'_',spatialblend,'.nc'))
     phv.nc=nc_create(phv.fn,var)
     phvrcn.nc=nc_create(phvrcn.fn,var)
     tind=unique(fullpreds$yrdoy)
     
     d_ply(fullpreds,.(yrdoy),function(dF){
          ncvar_put(phv.nc, var, vals=dF$phv.fullpred,start=c(1,1,which(unique(dF$yrdoy)==tind)),count=c(-1,-1,1))
          ncvar_put(phvrcn.nc,var,vals=dF$phvrcn.fullpred,start=c(1,1,which(tind==unique(dF$yrdoy))),count=c(-1,-1,1))
          ncatt_put(phv.nc,0,'proj4string','+proj=longlat +datum=NAD83')
          ncatt_put(phvrcn.nc,0,'proj4string','+proj=longlat +datum=NAD83')
     })
     
     nc_close(phvnc)
     nc_close(phvrcnnc)
     
     #     
     #      
}