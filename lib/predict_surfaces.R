#use static model and best dynamic model to produce surface estimates for a doy. iterate by year. returns rasters  for stacking       
predict_surfaces=function(rdatedF,locs,newdata,newdatalocs,spatialblend,bpth){
#    rdatedF=which_recon_date[2,]#for testing only!
     # get reconyear and load the stack with that recon data.
     rdate=rdatedF$phvrcn_recondate
     mdate=rdatedF$date
     get_sca=function(date){
          yr=strftime(date,'%Y')
          if(as.numeric(strftime(date,'%m')) >= 3){
               rfn=paste0('data/recon_',recon.version,'/recondata_',yr,'_',recon.version,'.RData')
               if(!exists(paste0('recon',yr),envir=.GlobalEnv))  load(rfn,envir=.GlobalEnv)
               tmp=get(paste0('recon',yr),envir=.GlobalEnv)
          } else {     
               print('using modscag fsca for masking. may not be completely cloudfree')
               tmp=stack(file.path('data','selectdates','modscag',paste0('fsca',yr)))    
          }
          stckdate=strftime(date,'X%Y%m%d')
          rind=grep(stckdate,names(tmp))
          return(getValues(tmp[[rind]]))
     }
     
     recon=get_sca(rdate)
     sca=get_sca(mdate)
     
     #combine ucophv dataframe and recon data
     newdata=cbind(newdata,recon)
     newdata$recon=scale(newdata$recon)#convoluted but necessary order of operations
     newdata=newdata[sca!=0,]
     scaind=which(sca!=0)
     
     # subset mdldata for date
     doydF=mdldata[mdldata$date==rdatedF$date,]
     
     # fit models for phv and phvrcn for recondate in rdatedF
     static_mdl=doPHVfit(doydF) #I guess this is happening twice...
     dyn_mdl=doPHVRCNfit(rdate,doydF,static_mdl)[[1]]
          
     #predict best models for whole domain
     static.pred=predict(static_mdl,newdata)# 
     dyn.pred=predict(dyn_mdl,newdata)         
     
     #setup dataframes for spatial analysis
     newdatapreds=data.frame(
          as.data.frame(newdatalocs)[sca!=0,],
          phv.predictions=static.pred,
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
     rasterlist=doSPATIAL(locpreds,locs, 'predict', newdatapreds,bpth,spatialblend)# newdatapreds has newdata locations and regression predictions. loc
          return(rasterlist)
     } else {
          pfn=paste0('fullpred_phv-',unique(rdatedF$date),'-',spatialblend,'.grd')
          prfn=paste0('fullpred_phvrcn-',unique(rdatedF$date),'-',spatialblend,'.grd')
          #
          writeRaster(num2ucoRaster(static.pred),
                      filename=file.path(bpn,pfn),overwrite=T)
          phv.fullpreds=raster(file.path(bpn,pfn))
          writeRaster(num2ucoRaster(dyn.pred),
                      filename=file.path(bpn,pfn),overwrite=T)
          phvrcn.fullpreds=raster(file.path(bpn,prfn))
          
          return(list(phv.fullpreds,phvrcn.fullpreds))
     }
}
