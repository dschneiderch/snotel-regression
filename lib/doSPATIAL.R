#takes estimates residuals at snotel locations and fits krige model. substracts krige estimates from predictions (either at snotel snotellocs.usgs or for entire domain)
doSPATIAL=function(dF, snotellocs.usgs,covrange, app, newdatapreds=NULL,newdatalocs.agg.usgs,bpth,spatialblend){
     # print(str(dF))
     dF=mutate(dF,
          phvresid=phv.predictions-snotel,#
          phvrcnresid=phvrcn.predictions-snotel,
          # reconresid=recon-snotel,
          reconrtresid=reconrt-snotel)
    locsmat=as.matrix(as.data.frame(snotellocs.usgs)[,c('x','y')])
    covrange=as.numeric(gsub('km','',covrange))*1000
    rlambda=summary(Krig(locsmat,dF$reconrtresid,theta=covrange,cov.function="wendland.cov"))$lambda
    rresidkr=mKrig(locsmat,dF$reconrtresid,theta=covrange,lambda=rlambda,cov.function="wendland.cov")
    slambda=summary(Krig(locsmat,dF$phvresid,theta=covrange,cov.function="wendland.cov"))$lambda     
    sresidkr=mKrig(locsmat,dF$phvresid,theta=covrange,lambda=slambda,cov.function="wendland.cov")
    dlambda=summary(Krig(locsmat,dF$phvrcnresid,theta=covrange,cov.function="wendland.cov"))$lambda     
    dresidkr=mKrig(locsmat,dF$phvrcnresid,theta=covrange,lambda=dlambda,cov.function="wendland.cov")
    
    if(app=='xval'){
        phv.fullpred=dF$phv.predictions-predict(sresidkr,x=locsmat)
        phvrcn.fullpred=dF$phvrcn.predictions-predict(dresidkr,x=locsmat)
        reconrt.fullpred=dF$reconrt-predict(rresidkr,x=locsmat)

        snotelpreds= data.frame(
                        date=dF$date,
                        Station_ID=dF$Station_ID,
                        yrdoy=dF$yrdoy,
                        snotel=dF$snotel,
                        recondate=dF$recondate,
                        reconrt=dF$reconrt,
                        phvresid=dF$phvresid,
                        phvrcnresid=dF$phvrcnresid,
                        reconrtresid=dF$reconrtresid,
                        phv.pred=dF$phv.predictions,
                        phvrcn.pred=dF$phvrcn.predictions, 
                        phv.fullpred,
                        phvrcn.fullpred,
                        reconrt.fullpred)
          return(snotelpreds)
     }
     if(app=='predict'){
          #bpn=file.path(bpn,'raster')
#           pfn=paste0('fullpred_phv-',unique(dF$date),'-',spatialblend,'.grd')
#           prfn=paste0('fullpred_phvrcn-',unique(dF$date),'-',spatialblend,'.grd')
          xpred=as.data.frame(newdatalocs.agg.usgs)[,c('x','y')]
          grid.list=list(x=newdatapreds[,'x'],y=newdatapreds[,'y'])
          
          predagg=predict(sresidkr,x=xpred)
          obj=list(x=xpred[,'x'],y=xpred[,'y'],z=predagg)
          phv.fullpred=newdatapreds$phv.predictions-interp.surface.grid(obj,grid.list)
          
          predagg=predict(dresidkr,x=xpred)
          obj=list(x=xpred[,'x'],y=xpred[,'y'],z=predagg)
          phvrcn.fullpred=newdatapreds$phvrcn.predictions-interp.surface.grid(obj,grid.list)
          return(
               data.frame(
                    phv.predictions=newdatapreds$phv.predictions,
                    phvrcn.predictions=newdatapreds$phvrcn.predictions,
                    phv.fullpred,
                    phvrcn.fullpred)
               )
#           rs=interpolate(rskel,sresidkr,xyOnly=T)
#           rs=projectRaster(rs,geope,method='bilinear',crs=projection(geope))
#           rpred=num2ucoRaster(newdatapreds$phv.predictions)
#           #plot(rpred-rs)
#           writeRaster(rpred-rs,filename=file.path(bpn,pfn),overwrite=T)
#           phv.fullpreds=raster(file.path(bpn,pfn))
# #              
#           dresidkr=Krig(locsmat,dF$phvrcnresid,theta=200000)
#           rs=interpolate(rskel,dresidkr,xyOnly=T)
#           rs=projectRaster(rs,geope,method='bilinear',crs=projection(geope))
#           rpred=num2ucoRaster(newdatapreds$phvrcn.predictions)
#           #plot(rpred-rs)
#           writeRaster(rpred-rs,filename=file.path(bpn,prfn),overwrite=T)
#           phvrcn.fullpreds=raster(file.path(bpn,prfn))
#           #
#           return(list(phv.fullpreds, phvrcn.fullpreds))
#          
     }
}