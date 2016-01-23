# setwd('/Users/dosc3612/Documents/snotel-regression_project')
library(ProjectTemplate)
reload.project(list(cache_loading=F))

library(raster)
library(reshape2)
library(plyr)
library(ggplot2)
library(doMC)
library(RColorBrewer)

get_fsca=function(dte){
     yr=strftime(dte,'%Y')
     if( as.numeric(strftime(dte,'%m')) >= 3){
          rfn=paste0('data/recon_',recon.version,'/reconfscadata_',yr,'_',recon.version,'.nc')
          ncstack=stack(rfn)
          stckdate=strftime(dte,'X%Y%m%d')
          rind=grep(stckdate,names(ncstack))
          r=ncstack[[rind]]
     } else {
          print('using modscag fsca for snotel pixel sca. may not be completely cloudfree')
          sfn=file.path('data','selectdates','modscag',paste0('fsca',yr,'.nc'))
          ncstack=stack(sfn)
          stckdate=strftime(dte,'X%Y%j')
          rind=grep(stckdate,names(ncstack))
          r=tryCatch({ncstack[[rind]]},error=function(x) {#this error should never happen
               print(paste0('the date ',dte,' was not available in the selected sca image'))})
          r=shift(r,0.0041666667/2,-0.004166667/2)
          valind <- (r<=100 & r!=0)
          # plot(valind)
          r[valind]=r[valind]/100
     }
     return(r)
}

get_modelswe=function(metadata,surf,gridding){
     ## create polygon of extent of survey area
     e <- extent(surf)
     box <- as(e, "SpatialPolygons")
     proj4string(box)=CRS(projection(surf))

     # surfpoly=rasterToPolygons(r2,dissolve=T)
     surfpoly=spTransform(box,CRS('+init=epsg:4326'))

     dte=as.POSIXct(strptime(names(surf),'X%Y%m%d',tz='MST'))
     # if(dte==as.POSIXct('2001-04-23',tz='MST')) dte=as.POSIXct('2001-04-22',tz='MST') #upper san juan
     # if(dte==as.POSIXct('2001-04-27',tz='MST')) dte=as.POSIXct('2001-04-26',tz='MST') #
     # if(dte==as.POSIXct('2001-04-27',tz='MST')) dte=as.POSIXct('2001-04-26',tz='MST') #
     # if(dte==as.POSIXct('2008-04-05',tz='MST')) dte=as.POSIXct('2008-04-04',tz='MST')
     if(dte==as.POSIXct('2009-02-27',tz='MST')) dte=as.POSIXct('2009-02-28',tz='MST') #joe wright survey
     #if(dte==as.POSIXct('2002-04-03',tz='MST')) dte=as.POSIXct('2002-04-02',tz='MST') #lilly pond
     #if(dte==as.POSIXct('2002-04-05',tz='MST')) dte=as.POSIXct('2002-04-04',tz='MST') #upper san juan
     # if(dte==as.POSIXct('2002-04-06',tz='MST')) dte=as.POSIXct('2002-04-04',tz='MST') #
     print(paste(dte,metadata[1]))
     yr=strftime(dte,'%Y')

     ## snodas
     layername=paste0('X',strftime(dte,'%Y%m%d'))
     newlayername=strftime(dte,'X%Y%j')
     if(yr>2003){
          snodas=stack(paste0('data/spatialvalidation/snodas/snodas.survey.nc'))
          layerind=which(names(snodas) %in% layername)
          snodassub=snodas[[layerind]]
          snodas.crop=crop(snodassub,extent(surfpoly)*1.2)
          snodaspoly=spTransform(rasterToPolygons(snodas.crop),CRS(proj4string(surf)))

          # # print(metadata[[2]])
          # png(filename=paste0('data/spatialvalidation/',metadata[[2]],'.png'))
          # plot(surf,main=metadata,colNA='grey80')
          # phvpoly=spTransform(rasterToPolygons(phv.crop),proj4string(surf))
          # plot(phvpoly,add=T)
          # plot(snodaspoly,add=T,lty=3)
          # dev.off()

          vm=extract(snodas.crop,snodaspoly)
          vo=extract(surf,snodaspoly,cellnumbers=F)
          snodasswe=ldply(as.list(1:length(vm)), function(i){
               x=vo[[i]]
               if (!is.null(x)) {
                    data.frame(cell=i,swe.obs=x,swe.model=vm[[i]]/1000,model.type='SNODAS')
               } else { NULL}
          })
     }


     if(scalesnotel=='noscale' && postscaled=='_postscaled'){
          fsca=get_fsca(dte)
          forden=get_forden('umd_forden',yr)
          forden[forden>.99] = .99
          fscaind <- (fsca<=1 & fsca!=0)
          badind <- (fsca>200)
          fsca_vgfcorrected <- fsca
          fsca_vgfcorrected[fscaind] <- fsca[fscaind] / (1 - forden[fscaind])
          fsca_vgfcorrected[fsca_vgfcorrected>1]=1
          # print(fsca_vgfcorrected)
          fsca_vgfcorrected[badind]=NA#fsca[badind]
          # print(fsca_vgfcorrected)
     }

     ## PHV
     layername=strftime(dte,'X%Y%m%d')
     ## phv
     #      print('--phv')
     phv=stack(paste0(basepath,residblend,'preds-phv_',yr,'_blend-',fscaMatch,'.nc'))
     layerind=which(names(phv) %in% layername)
     phvsub=phv[[layerind]]
     if(scalesnotel=='noscale' & postscaled=='_postscaled') phvsub=phvsub*fsca_vgfcorrected
     phv.crop=crop(phvsub,extent(surfpoly)*1.2)
     if(yr>2003 & grepl('snodasgrid',gridding)) phv.crop=projectRaster(phv.crop,snodas.crop)
     phvpoly=spTransform(rasterToPolygons(phv.crop),CRS(proj4string(surf)))
     vm=extract(phv.crop,phvpoly)
     vo=extract(surf,phvpoly,cellnumbers=F,weights=F)
     phvswe=ldply(as.list(1:length(vm)), function(i){
          x=vo[[i]]
          if (!is.null(x)) {
               data.frame(cell=i,swe.obs=x,swe.model=vm[[i]],model.type='PHV')
          } else { NULL}
     })

     ## phvrcn
     #    print('--phvrcn')
     phvrcn=stack(paste0(basepath,residblend,'preds-phvrcn_',yr,'_blend-',fscaMatch,'.nc'))
     layerind=which(names(phvrcn) %in% layername)
     phvrcnsub=phvrcn[[layerind]]
     if(scalesnotel=='noscale' & postscaled=='_postscaled') phvrcnsub=phvrcnsub*fsca_vgfcorrected
     phvrcn.crop=crop(phvrcnsub,extent(surfpoly)*1.2)
     if(yr>2003 & grepl('snodasgrid',gridding)) phvrcn.crop=projectRaster(phvrcn.crop,snodas.crop)
     phvrcnpoly=spTransform(rasterToPolygons(phvrcn.crop),CRS(proj4string(surf)))
     vm=extract(phvrcn.crop,phvrcnpoly)
     vo=extract(surf,phvrcnpoly,cellnumbers=F)
     phvrcnswe=ldply(as.list(1:length(vm)), function(i){
          x=vo[[i]]
          if (!is.null(x)) {
               data.frame(cell=i,swe.obs=x,swe.model=vm[[i]],model.type='PHVRCN')
          } else { NULL}
     })

     if(predictor=='fsca' && fscaMatch=='fsca'){

          ## phvfsca
          #    print('--phvrcn')
          phvfsca=stack(paste0(basepath,residblend,'preds-phvfsca_',yr,'_blend-',fscaMatch,'.nc'))
          layerind=which(names(phvrcn) %in% layername)
          phvfscasub=phvfsca[[layerind]]
          if(scalesnotel=='noscale' && postscaled=='_postscaled') phvrcnsub=phvfscasub*fsca_vgfcorrected
          phvfsca.crop=crop(phvrcnsub,extent(surfpoly)*1.2)
          if(yr>2003 && grepl('snodasgrid',gridding)) phvfsca.crop=projectRaster(phvfsca.crop,snodas.crop)
          phvfscapoly=spTransform(rasterToPolygons(phvfsca.crop),CRS(proj4string(surf)))
          vm=extract(phvfsca.crop,phvfscapoly)
          vo=extract(surf,phvfscapoly,cellnumbers=F)
          phvfscaswe=ldply(as.list(1:length(vm)), function(i){
               x=vo[[i]]
               if (!is.null(x)) {
                    data.frame(cell=i,swe.obs=x,swe.model=vm[[i]],model.type='PHVFSCA')
               } else { NULL}
          })
     } else {
          phvfscaswe=phvswe
          phvfscaswe$model.type='PHVFSCA'
          phvfscaswe$swe.model=NA
     }
     ## recon
     #  print('--rcn')
     layername=paste0('X',strftime(dte,'%Y%m%d'))
     newlayername=strftime(dte,'X%Y%j')
     if(as.numeric(strftime(dte,'%m'))>=3){
          recon=stack(paste0('data/recon_',recon.version,'/recondata_',yr,'_',recon.version,'.nc'))
          layerind=which(names(recon) %in% layername)
          reconsub=recon[[layerind]]
          if(scalesnotel=='noscale' & postscaled=='_postscaled') reconsub=reconsub*fsca_vgfcorrected
          recon.crop=crop(reconsub,extent(surfpoly)*1.2)
          if(yr>2003 & grepl('snodasgrid',gridding)) recon.crop=projectRaster(recon.crop,snodas.crop)
          reconpoly=spTransform(rasterToPolygons(recon.crop),CRS(proj4string(surf)))
          vm=extract(recon.crop,reconpoly)
          vo=extract(surf,reconpoly,cellnumbers=F)
          reconswe=ldply(as.list(1:length(vm)), function(i){
               x=vo[[i]]
               if (!is.null(x)) {
                    data.frame(cell=i,swe.obs=x,swe.model=vm[[i]],model.type='RCN')
               } else { NULL}
          })
     } else {
          reconswe=phvswe
          reconswe$model.type='RCN'
          reconswe$swe.model=NA
     }

     if(yr<=2003){
          snodasswe=phvswe
          snodasswe$model.type='SNODAS'
          snodasswe$swe.model=NA
     }

     mdlswe=rbind(phvswe,phvrcnswe,phvfscaswe,reconswe,snodasswe)
     mdlswe$date=dte#as.POSIXct(strptime(metadata[3],'%Y%m%d',tz='MST'))
     mdlswe$site=metadata[1]
     mdlswe$yr=strftime(mdlswe$date,'%Y')
     mdlswe$mth=strftime(mdlswe$date,'%m')
     mdlswe$swe.model[mdlswe$swe.model>10]=NA
     return(mdlswe)
}


recon.version='v3.1'
covrange='idp1'
cost='r2'
residblend=''
scalesnotel='scale'
fscaMatch='fsca'
dateflag='surveyvalidation'
style='real-time'#nalysis'
gridding='_modisgrid'#'_snodasgrid' or '_modisgrid'
postscaled=''#'' or '_postscaled'
predictor='fsca'

for(postscaled in c('')){
     for(style in c('real-time')){
          for(dateflag in c('surveyvalidation')){
               for(residblend in c('')){
                    # scalesnotel='noscale'
                    for(scalesnotel in c('scale','noscale')){
                         # fscaMatch='wofsca'
                         if(scalesnotel=='scale' && postscaled=='_postscaled') next

                         for(fscaMatch in c('wofsca','fsca')){
                              # if(cost=='rmse' && fscaMatch=='fsca'){
                              #   next#rmse/fsca and r2/fsca will have the same predictions
                              # }

                              config=''#'bic-v2.5-removed_NWbdiff'
                              basepath=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',scalesnotel,'/',dateflag,'/fullpreds/',cost,'/netcdf/',style,'/',config,'/')

                              ## CSU
                              ## ----
                              basein='data/spatialvalidation/csu_COWY_snowsurvey0809/Meromyetal_2'
                              surflist=list(c('Niwot','nis08apr_swe','20080407'),
                                            c('Niwot','nis08may_swe','20080505'),
                                            c('Niwot','nis09mar_swe','20090306'),
                                            c('Lizard Head','lh08mar_swe','20080316'),
                                            c('South Brush Creek','bcs08apr_swe','20080405'),
                                            c('South Brush Creek','bcs09mar_swe','20090329'),
                                            c('Dry Lake','dl08apr_swe','20080404'),
                                            c('Dry Lake','dl08may_swe','20080502'),
                                            c('Dry Lake','dl09feb_swe','20090228'),
                                            c('Dry Lake','dl09mar_swe','20090328'),
                                            c('Joe Wright','jw08apr_swe','20080403'),
                                            c('Joe Wright','jw08may_swe','20080501'),
                                            c('Joe Wright','jw09feb_swe','20090227'),
                                            c('Joe Wright','jw09may_swe','20090502'))

                              # surflist=surflist[c(12,13,14)]

                              surfrasters=lapply(surflist,function(x){
                                   r=raster(file.path(basein,paste0(x[2],'.tif')))
                                   names(r)=x[3]
                                   return(r)
                                   #  assign(paste(x[1],x[3],'meromysurf',sep='-'),r,envir=.GlobalEnv)
                              })
                              # metadata=surflist[[2]]
                              # surf=surfrasters[[2]]

                              if(grepl('colorado.edu',Sys.getenv("HOSTNAME"))){
                                   registerDoMC(8)
                              } else {
                                   registerDoMC(3)
                              }
                              
                              csuswe=ldply(seq(1:length(surflist)),function(x) get_modelswe(surflist[[x]], surfrasters[[x]], gridding), .parallel=T,.inform=F)
                              csuswe$swe.obs=csuswe$swe.obs/100

                              if( !grepl('snodas',gridding)) {
                                   ########  urg surveys
                                   basein='data/spatialvalidation/urg.validation/Rio_Grande_Headwaters/interpolated_depths_swe'
                                   surflist=list(c('Lilly Pond','sw_lilly0401','20010427'),
                                                 c('Lilly Pond','sw_lilly0402','20020403'),
                                                 c('Slumgullion','sw_slum0401','20010422'),
                                                 c('Upper San Juan','sw_usj0401','20010423'),
                                                 c('Wolf Creek Summit','sw_wc0401','20010424'),
                                                 c('Slumgullion','sw_slum0402','20020406'),
                                                 c('Upper San Juan','sw_usj0402','20020405'),
                                                 c('Wolf Creek Summit','sw_wc0402','20020404'))

                                   surfrasters=lapply(surflist,function(x){
                                        r=raster(file.path(basein,x[2]))
                                        names(r)=x[3]
                                        return(r)
                                        #  assign(paste(x[1],x[3],'meromysurf',sep='-'),r,envir=.GlobalEnv)
                                   })
                                   
                                   urgswe=ldply(seq(1:length(surflist)),function(x) get_modelswe(surflist[[x]], surfrasters[[x]], gridding), .parallel=T,.inform=F)
                                   urgswe$swe.obs=urgswe$swe.obs/100
                              }

                              ##### GLV
                              glvdates=read.table('data/spatialvalidation/glv.validation/snow_survey_dates.txt',sep='\n',stringsAsFactor=F)
                              colnames(glvdates)='date'
                              glvdates$date=as.POSIXct(glvdates$date,tz='MST')
                              glvdates$yr=strftime(glvdates$date,'%Y')
                              glvdates$mth=strftime(glvdates$date,'%m')
                              if(grepl('snodas',gridding)) glvdates=glvdates[glvdates$yr>2003,]
                              #--
                              glvdensity=read.table('data/spatialvalidation/glv.validation/snow_survey_densities.txt',sep=',',header=F,stringsAsFactor=F)
                              colnames(glvdensity)=c('date','density')
                              glvdensity$date=as.POSIXct(glvdensity$date,tz='MST')
                              glvdensity$yr=strftime(glvdensity$date,'%Y')
                              #
                              glvsurvey=merge(glvdates,glvdensity)
                              glvsurvey=subset(glvsurvey,yr>2000 & yr <= 2012)
                              glvsurvey=glvsurvey[!is.na(glvsurvey$density),]
                              # glvsurvey=subset(glvsurvey,yr<2003)

                              glvpath='data/spatialvalidation/glv.validation/swe_surfaces'
                              glvswe=ddply(glvsurvey,.(yr),.parallel=F,function(dF){
                                   yr=dF$yr
                                   dte=dF$date
                                   # print(yr)
                                   fn=file.path(glvpath,paste0('glv_30m_erickson_swe_wmsk_',yr,'.tif'))
                                   glv=raster(fn)
                                   glv=crop(glv,extent(c(444345.1,448000,4432000,4436000)))
                                   names(glv)=strftime(dte,'X%Y%m%d',tz='MST')
                                   glvswe=get_modelswe(c('Green Lakes Valley',paste('glv',yr,sep='-'),strftime(dte,'%Y%m%d')),glv,gridding)
                                   return(glvswe)
                              })


                              if( !grepl('snodas',gridding)) {
                                   swe=rbind(csuswe,urgswe,glvswe)
                              } else {
                                   swe=rbind(csuswe,glvswe)
                              }
                              if(residblend=='full') saveRDS(swe,file=paste0('data/spatialvalidation/blended/',dateflag,'/',style,'/surveygriddedswe',gridding,'_rswe',recon.version,'_',covrange,'_snotel',scalesnotel,postscaled,'_',cost,'_',fscaMatch,'-withphvfsca.rds'))
                              if(residblend=='') saveRDS(swe,file=paste0('data/spatialvalidation/unblended/',dateflag,'/',style,'/surveygriddedswe',gridding,'_rswe',recon.version,'_',covrange,'_snotel',scalesnotel,postscaled,'_',cost,'_',fscaMatch,'-withphvfsca.rds'))

                         }}}}}}
