#use static model and best dynamic model to produce surface estimates for a doy. iterate by year. returns rasters  for stacking       
predict_surfaces=function(rdatedF,snotellocs.usgs,newdata,newdatalocs,newdatalocs.agg.usgs,spatialblend,bpth){
     #    rdatedF=which_recon_date[2,]#for testing only!

          # get reconyear and load the stack with that recon data.
 mdate=rdatedF$date
 rdate=rdatedF$phvrcn_recondate
 print(paste0('-- model date: ',mdate))
 print(paste0('---- recon date: ',rdate))

 if(!is.na(rdate)){

  get_sca=function(dte){
    yr=strftime(dte,'%Y')
    if(dte<max(snotelrecon$date) & as.numeric(strftime(dte,'%m')) >= 3){
     rfn=paste0('data/recon_',recon.version,'/recondata_',yr,'_',recon.version,'.nc')
     ncstack=stack(rfn)
     stckdate=strftime(dte,'X%Y%m%d')
    } else {     
     print('using modscag fsca for masking. may not be completely cloudfree')
     sfn=file.path('data','selectdates','modscag',paste0('fsca',yr,'.nc'))
     ncstack=stack(sfn)
     stckdate=strftime(dte,'X%Y%j')
    }
    rind=grep(stckdate,names(ncstack))
    vals=tryCatch({getValues(ncstack[[rind]])},error=function(x) {#this error should never happen
    print(paste0('the date ',dte,' was not available in the selected sca image'))})
  }
 sca=get_sca(mdate)

 if(output!='points'){ 
   recon=get_sca(rdate)

   #save indices of sca to insert predictions
   scaind=which(!is.na(sca) & sca!=0 & sca<=100)#modscag is 0-100, recon is 0-peakswe (~8m)
   cloudind=which(sca>100)

   #combine ucophv dataframe and recon data
   newdata2=cbind(newdata,recon)
   newdata2$recon=scale(newdata2$recon)#convoluted but necessary order of operations
   newdata2=newdata2[scaind,]

 } else {
  newdata2=mdldata[mdldata$date==rdatedF$date,]
  newdatalocs.usgs=data.frame(as.data.frame(snotellocs.usgs)[,c('Station_ID','x','y')],recon=newdata2$recon)
  newdata2$recon=scale(newdata2$recon)
  sca=raster(matrix(sca,nrow=2580,byrow=T),xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
  projection(sca)='+proj=longlat +datum=WGS84'
  sca=extract(sca, snotellocs)#get values at snotel pixel locations. named sca so we can use its length later
  scaind=which(sca>0)#convert to binary. if cloudy then assume sca>0. this should only be jan-feb if it happens
  cloudind=NULL
  newdata2=newdata2[scaind,]
}
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
  reconrt=recondata[recondata$recondate==rdatedF$recon_costdate,'recon'])


newdatapredsfull=data.frame(phv.predictions=rep(0,length(sca)),phvrcn.predictions=0,phv.fullpred=0,phvrcn.fullpred=0)
#geostatistical blending of residuals
if(spatialblend=='blend'){
  fullpreds=doSPATIAL(locpreds,snotellocs.usgs, 'predict', newdatapreds, newdatalocs.agg.usgs, bpth, spatialblend)# newdatapreds has newdata locations and regression predictions. loc
  fullpreds[fullpreds<0]=0
  newdatapredsfull$phv.predictions[scaind]=fullpreds$phv.predictions
  newdatapredsfull$phvrcn.predictions[scaind]=fullpreds$phvrcn.predictions
  newdatapredsfull$phv.fullpred[scaind]=fullpreds$phv.fullpred
  newdatapredsfull$phvrcn.fullpred[scaind]=fullpreds$phvrcn.fullpred
  # newdatapredsfull=cbind(newdatapredsfull,as.data.frame(newdatalocs.usgs))            
} else {
  newdatapreds$phv.predictions[newdatapreds$phv.predictions<0]=0#individual masking because newdatapreds needs x y
  newdatapreds$phvrcn.predictions[newdatapreds$phvrcn.predictions<0]=0
#  
  # newdatapredsfull=data.frame(phv.fullpred=0,phvrcn.fullpred=0,as.data.frame(newdatalocs.usgs))
  newdatapredsfull$phv.fullpred[scaind]=newdatapreds$phv.predictions
  newdatapredsfull$phvrcn.fullpred[scaind]=newdatapreds$phvrcn.predictions
}
newdatapredsfull$phv.predictions[cloudind]=sca[cloudind]
newdatapredsfull$phvrcn.predictions[cloudind]=sca[cloudind]
newdatapredsfull$phv.fullpred[cloudind]=sca[cloudind]
newdatapredsfull$phvrcn.fullpred[cloudind]=sca[cloudind]
return(newdatapredsfull)    

} else {
  print(paste(mdate,'empty model'))
  return(data.frame())
}
}