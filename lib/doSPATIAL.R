#takes estimates residuals at snotel locations and fits krige model. substracts krige estimates from predictions (either at snotel snotellocs.usgs or for entire domain)
doSPATIAL=function(dF, snotellocs.usgs, app, newdata=NULL,bpth,spatialblend){
     dF=mutate(dF,
          phvresid=phv.predictions-snotel,#
          phvrcnresid=phvrcn.predictions-snotel,
          reconresid=recon-snotel)
     locsmat=as.matrix(as.data.frame(snotellocs.usgs)[,c('x','y')])
     #print(nrow(locsmat))
     #print(str(dF))
     #kr=Krig(locsmat,dF$phvresid,theta=200000)
     rresidkr=Krig(locsmat,dF$reconresid,theta=200000)
     sresidkr=Krig(locsmat,dF$phvresid,theta=200000)
     dresidkr=Krig(locsmat,dF$phvrcnresid,theta=200000)
     #residtps=Tps(locsmat,dF$phvresid)
     #print(residkr)
     if(app=='xval'){
          phv.fullpred=dF$phv.predictions-predict(sresidkr,x=locsmat)
          phvrcn.fullpred=dF$phvrcn.predictions-predict(dresidkr,x=locsmat)
          recon.fullpred=dF$recon-predict(rresidkr,x=locsmat)
          
          snotelpreds= data.frame(
                         date=dF$date,
                         yrdoy=dF$yrdoy,
                         snotel=dF$snotel,
                         recondate=dF$recondate,
                         recon=dF$recon,
                         phvresid=dF$phvresid,
                         phvrcnresid=dF$phvrcnresid,
                         reconresid=dF$reconresid,
                         phv.pred=dF$phv.predictions,
                         phvrcn.pred=dF$phvrcn.predictions, 
                         phv.fullpred,
                         phvrcn.fullpred,
                         recon.fullpred)
          return(snotelpreds)
     }
     if(app=='predict'){
          #bpn=file.path(bpn,'raster')
#           pfn=paste0('fullpred_phv-',unique(dF$date),'-',spatialblend,'.grd')
#           prfn=paste0('fullpred_phvrcn-',unique(dF$date),'-',spatialblend,'.grd')
          #
          phv.fullpred=newdata$phv.predictions-predict(sresidkr,x=newdata[,c('x','y')])
          phvrcn.fullpred=newdata$phvrcn.predictions-predict(dresidkr,x=newdata[,c('x','y')])        
          return(
               data.frame(
                    phv.predictions=newdata$phv.predictions,
                    phvrcn.predictions=newdata$phvrcn.predictions,
                    phv.fullpred,
                    phvrcn.fullpred)
               )
#           rs=interpolate(rskel,sresidkr,xyOnly=T)
#           rs=projectRaster(rs,geope,method='bilinear',crs=projection(geope))
#           rpred=num2ucoRaster(newdata$phv.predictions)
#           #plot(rpred-rs)
#           writeRaster(rpred-rs,filename=file.path(bpn,pfn),overwrite=T)
#           phv.fullpreds=raster(file.path(bpn,pfn))
# #              
#           dresidkr=Krig(locsmat,dF$phvrcnresid,theta=200000)
#           rs=interpolate(rskel,dresidkr,xyOnly=T)
#           rs=projectRaster(rs,geope,method='bilinear',crs=projection(geope))
#           rpred=num2ucoRaster(newdata$phvrcn.predictions)
#           #plot(rpred-rs)
#           writeRaster(rpred-rs,filename=file.path(bpn,prfn),overwrite=T)
#           phvrcn.fullpreds=raster(file.path(bpn,prfn))
#           #
#           return(list(phv.fullpreds, phvrcn.fullpreds))
#          
     }
}