#takes estimates residuals at snotel locations and fits krige model. substracts krige estimates from predictions (either at snotel locs or for entire domain)
doSPATIAL=function(dF, locs, app, newdata=NULL,bpth,spatialblend){
     dF=mutate(dF,
          phvresid=phv.predictions-snotel,#
          phvrcnresid=phvrcn.predictions-snotel,
          reconresid=recon-snotel)
     locsmat=as.matrix(as.data.frame(locs)[,c('x','y')])
     #print(nrow(locsmat))
     #print(str(dF))
     #kr=Krig(locsmat,dF$phvresid,theta=200000)
     sresidkr=Krig(locsmat,dF$phvresid,theta=200000)

     #dresidkr=Krig(locsmat,dF$phvrcnresid,theta=200000)
     #residtps=Tps(locsmat,dF$phvresid)
     #print(residkr)
     if(app=='xval'){
          x=locsmat
          phv.fullpred=dF$phv.predictions-predict(sresidkr,x=x)
          phvrcn.fullpred=dF$phvrcn.predictions-predict(sresidkr,x=x,y=dF$phvrcnresid)
          return(
               data.frame(
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
                    phvrcn.fullpred)
          )
     }
     if(app=='predict'){
          bpn=file.path(bpth,'raster')
          pfn=paste0('fullpred_phv-',unique(dF$date),'-',spatialblend,'.grd')
          prfn=paste0('fullpred_phvrcn-',unique(dF$date),'-',spatialblend,'.grd')
          #
          rs=interpolate(rskel,sresidkr,xyOnly=T)
          rs=projectRaster(rs,geope,method='bilinear',crs=projection(geope))
          rpred=num2ucoRaster(newdata$phv.predictions)
          #plot(rpred-rs)
          writeRaster(rpred-rs,filename=file.path(bpn,pfn),overwrite=T)
          phv.fullpreds=raster(file.path(bpn,pfn))
#              
          dresidkr=Krig(locsmat,dF$phvrcnresid,theta=200000)
          rs=interpolate(rskel,dresidkr,xyOnly=T)
          rs=projectRaster(rs,geope,method='bilinear',crs=projection(geope))
          rpred=num2ucoRaster(newdata$phvrcn.predictions)
          #plot(rpred-rs)
          writeRaster(rpred-rs,filename=file.path(bpn,prfn),overwrite=T)
          phvrcn.fullpreds=raster(file.path(bpn,prfn))
          #
          return(list(phv.fullpreds, phvrcn.fullpreds))
#           x=as.matrix(as.data.frame(newdata)[,c('x','y')])       		  									    
#           phv.fullpred=newdata$phv.predictions-predict(residkr,x=x)
#           phvrcn.fullpred=newdata$phvrcn.predictions-predict(residkr,x=x,y=dF$phvrcnresid)        
#           return(
#                data.frame(
#                     phv.fullpred,
#                     phvrcn.fullpred
#                     )
#                )
     }
}