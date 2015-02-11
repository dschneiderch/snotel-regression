#use static model and best dynamic model to produce surface estimates for a doy. iterate by year. returns rasters  for stacking       
predict_surfaces=function(rdatedF,snotellocs.usgs,newdata,newdatalocs,spatialblend,bpth){
     #    rdatedF=which_recon_date[2,]#for testing only!
     # get reconyear and load the stack with that recon data.
     rdate=rdatedF$phvrcn_recondate
     mdate=rdatedF$date
if(!is.na(mdate)){
     get_sca=function(dte){
          yr=strftime(dte,'%Y')
          if(dte<max(recondata$recondate) & as.numeric(strftime(dte,'%m')) >= 3){
               rfn=paste0('data/recon_',recon.version,'/recondata_',yr,'_',recon.version,'.RData')
               if(!exists(paste0('recon',yr),envir=.GlobalEnv))  load(rfn,envir=.GlobalEnv)
               tmp=get(paste0('recon',yr),envir=.GlobalEnv)
               stckdate=strftime(dte,'X%Y%m%d')
          } else {     
               print('using modscag fsca for masking. may not be completely cloudfree')
               tmp=stack(file.path('data','selectdates','modscag',paste0('fsca',yr,'.nc')))    
               stckdate=strftime(dte,'X%Y%j')
          }
          rind=grep(stckdate,names(tmp))
          vals=tryCatch({getValues(tmp[[rind]])},error=function(x) {
          print(paste0('the date ',dte,' was not available in the selected sca image'))})
     }
     
     recon=get_sca(rdate)
     sca=get_sca(mdate)
     
     #save indices of sca to insert predictions
     scaind=which(sca!=0 & sca<=100)
     cloudind=which(sca>100)
     
     #combine ucophv dataframe and recon data
     newdata2=cbind(newdata,recon)
     newdata2$recon=scale(newdata2$recon)#convoluted but necessary order of operations
     newdata2=newdata2[scaind,]

     # subset mdldata for date
     doydF=mdldata[mdldata$date==rdatedF$date,]
                     
     # fit models for phv and phvrcn for recondate in rdatedF
     static_mdl=doPHVfit(doydF) #I guess this is happening twice...
     dyn_mdl=doPHVRCNfit(rdate,doydF,static_mdl)[[1]]
                     
     #predict best models for whole domain
     static.pred=predict(static_mdl,newdata2)# 
     dyn.pred=predict(dyn_mdl,newdata2)        
                     
     #setup dataframes for spatial analysis
     newdatapreds=data.frame(
          as.data.frame(newdatalocs.usgs)[scaind,],
          phv.predictions=static.pred,#newdata was already subset for sca above
          phvrcn.predictions=dyn.pred)
     locpreds=data.frame(
          date=doydF$date,
          yrdoy=doydF$yrdoy,
          phv.predictions=predict(static_mdl),
          phvrcn.predictions=predict(dyn_mdl),
          snotel=doydF$snotel,
          recon=doydF$recon)
                     
     
     
     #geostatistical blending of residuals
     if(spatialblend=='blend'){
          fullpreds=doSPATIAL(locpreds,snotellocs.usgs, 'predict', newdatapreds, bpth, spatialblend)# newdatapreds has newdata locations and regression predictions. loc
          fullpreds[fullpreds<0]=NA
          newdatapredsfull=data.frame(phv.predictions=rep(NA,length(sca)),phvrcn.predictions=NA,phv.fullpred=NA,phvrcn.fullpred=NA)
          newdatapredsfull$phv.predictions[scaind]=fullpreds$phv.predictions
          newdatapredsfull$phvrcn.predictions[scaind]=fullpreds$phvrcn.predictions
          newdatapredsfull$phv.fullpred[scaind]=fullpreds$phv.fullpred
          newdatapredsfull$phvrcn.fullpred[scaind]=fullpreds$phvrcn.fullpred
          if(length(cloudind)>0) {
               newdatapredsfull=lapply(newdatapredsfull,function(col) {
                    col[cloudind]=sca[cloudind]
                    return(col)
               })
          }
          newdatapredsfull=cbind(newdatapredsfull,as.data.frame(newdatalocs.usgs))
          
#           pfn=paste0('fullpred_phv-',unique(rdatedF$date),'-',spatialblend,'.grd')
#           prfn=paste0('fullpred_phvrcn-',unique(rdatedF$date),'-',spatialblend,'.grd')
          
          #           writeRaster(num2ucoRaster(newdatapredsfull$phv.fullpred),
          #                             filename=file.path(bpth,pfn),overwrite=T)
          #                           #          phv.fullpreds=raster(file.path(bpth,pfn))
          #           writeRaster(num2ucoRaster(newdatapredsfull$phvrcn.fullpred),
          #                             filename=file.path(bpth,pfn),overwrite=T)
          #                           #          phvrcn.fullpreds=raster(file.path(bpth,prfn))                
     
       #          phvrcn.fullpreds=raster(file.path(bpth,prfn))                
     } else {
          newdatapreds[newdatapreds<0]=NA
          newdatapredsfull=data.frame(as.data.frame(newdatalocs),phv.predictions=NA,phvrcn.predictions=NA)
          newdatapredsfull$phv.predictions[scaind]=newdatapreds$phv.predictions
          newdatapredsfull$phvrcn.predictions[scaind]=newdatapreds$phvrcn.predictions
     }
                     
#      pfn=paste0('glmpred_phv-',unique(rdatedF$date),'.grd')
#      prfn=paste0('glmpred_phvrcn-',unique(rdatedF$date),'.grd')
#                     #
#      writeRaster(num2ucoRaster(newdatapredsfull$phv.predictions),
#                        filename=file.path(bpth,pfn),overwrite=T)
#      #phv.fullpreds=raster(file.path(bpth,pfn))
#      writeRaster(num2ucoRaster(newdatapredsfull$phvrcn.predictions),
#                        filename=file.path(bpth,pfn),overwrite=T)
     #phvrcn.fullpreds=raster(file.path(bpth,prfn))
                                      
     return(newdatapredsfull)
                     #           pfn=paste0('fullpred_phv-',unique(rdatedF$date),'-',spatialblend,'.grd')
                     #           prfn=paste0('fullpred_phvrcn-',unique(rdatedF$date),'-',spatialblend,'.grd')
                     #           #
                     #           writeRaster(num2ucoRaster(static.pred),
                     #                       filename=file.path(bpth,pfn),overwrite=T)
                     #           phv.fullpreds=raster(file.path(bpth,pfn))
                     #           writeRaster(num2ucoRaster(dyn.pred),
                     #                       filename=file.path(bpth,pfn),overwrite=T)
                     #           phvrcn.fullpreds=raster(file.path(bpth,prfn))
                     #           
                     #           return(list(phv.fullpreds,phvrcn.fullpreds))
    
} else {
     return(data.frame())
}
}