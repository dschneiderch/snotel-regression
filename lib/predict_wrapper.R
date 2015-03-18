#allows iteration of predict_surfaces by year to create yearly stacks and ncdf output
predict_wrapper <- function(rdtsdF,style,newdata,newdatalocs,newdatalocs.agg.usgs,spatialblend){

     registerDoMC(allocated.cores)
     print(paste('**',unique(rdtsdF$yr),'**'))
     bpth=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/fullpreds')

     newdatapredsfull=ddply(rdtsdF,.(yrdoy),predict_surfaces,snotellocs.usgs,newdata,newdatalocs,newdatalocs.agg.usgs,spatialblend,bpth,.inform=T,.parallel=F)#,.paropts=list(.export=ls(), .packages=.packages(all=T)))     

if(output!='points'){
     watermask=getValues(raster('data/cs_NHD_MOD44_water_mask.tif'))
     # wm=raster('data/cs_NHD_MOD44_water_mask.tif')
     # plot(rr,zlim=c(1,1),colNA='black')            
     
     #output rasters of each model type.
     dim1=ncdim_def('Long','degree',seq(-112.25,-104.125,0.00416666667))
     dim2=ncdim_def('Lat','degree',seq(43.75,33,-0.00416666667))
     dim3=ncdim_def('time','yrdoy',unlim=T,vals=wrd$yrdoy)#as.numeric(unique(newdatapredsfull$yrdoy)))
     var=ncvar_def('swe','meters',dim=list(dim1,dim2,dim3),missval=-99,longname='snow water equivalent',compression=6)
     yr=unique(rdtsdF$yr)
     phv.fn=file.path(bpth,'netcdf',paste0('fullpreds-phv_',yr,'_',spatialblend,'.nc'))
     phvrcn.fn=file.path(bpth,'netcdf',paste0('fullpreds-phvrcn_',yr,'_',spatialblend,'.nc'))
     phv.nc=nc_create(phv.fn,var)
     phvrcn.nc=nc_create(phvrcn.fn,var)
     tind=unique(newdatapredsfull$yrdoy)

     d_ply(newdatapredsfull,.(yrdoy),function(dF){
               #do some data gymnastics to get it in the right order. for dim1 and dim2 min to max, data should be left to right, bottom to top
          val=dF$phv.fullpred*watermask
          ncvar_put(phv.nc, var, vals=val,start=c(1,1,which(unique(dF$yrdoy)==tind)),count=c(-1,-1,1))
               #
          val=dF$phvrcn.fullpred*watermask
          ncvar_put(phvrcn.nc,var,vals=val,start=c(1,1,which(tind==unique(dF$yrdoy))),count=c(-1,-1,1))
               #
          ncatt_put(phv.nc,0,'proj4string','+proj=longlat +datum=WGS84')
          ncatt_put(phvrcn.nc,0,'proj4string','+proj=longlat +datum=WGS84')
     })
     nc_close(phv.nc)
     nc_close(phvrcn.nc)
} else {
     if(nrow(newdatapredsfull)%%237 != 0) stop()
     # xy=grep('^x$|^y$',names(newdatapredsfull)) 
     # newdatapredsfull=newdatapredsfull[,-xy]
     fn2w=paste0(bpth,'/snotellocs/snotel_predictions_',unique(rdtsdF$yr),'.txt')
     write.table(x=newdatapredsfull,file=fn2w,sep='\t',row.names=F,quote=F)

}
}