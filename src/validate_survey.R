setwd('~/Documents/snotel-regression_project')
library(raster)
library(reshape2)
library(plyr)
library(ggplot2)

rswe='v3.1'
snotelscale='scale'
get_modelswe=function(dF,swesp){
      yrdoy=unique(dF$yrdoy)
      print(yrdoy)
      swesp=swesp[swesp$yrdoy==yrdoy,]##needs to be first
      
      if(yrdoy==2001113) yrdoy=2001112
      if(yrdoy==2001117) yrdoy=2001116
      if(yrdoy==2002093) yrdoy=2002092
      if(yrdoy==2002095) yrdoy=2002094
      if(yrdoy==2002096) yrdoy=2002097
      if(yrdoy==2008096) yrdoy=2008095
      if(yrdoy==2009058) yrdoy=2009059
                  
      yr=as.numeric(unique(dF$yr))
      ## phv
      phv=stack(paste0(basepath,'fullpreds-phv_',yr,'_blend.nc'))
      layername=paste0('X',yrdoy)
      layerind=which(names(phv) %in% layername)
      phvsub=phv[[layerind]]
      swesp$model.type='PHV'
      phvswe=extract(phvsub,swesp,sp=T,cellnumbers=T)
      names(phvswe)=c(names(swesp),'cell',layername)

      ## phvrcn
      phvrcn=stack(paste0(basepath,'fullpreds-phvrcn_',yr,'_blend.nc'))
      layername=paste0('X',yrdoy)
      layerind=which(names(phvrcn) %in% layername)
      phvrcnsub=phvrcn[[layerind]]
      swesp$model.type='PHVRCN'
      phvrcnswe=extract(phvrcnsub,swesp,sp=T,cellnumbers=T)
      names(phvrcnswe)=c(names(swesp),'cell',layername)

      ## recon
dte=as.POSIXct(strptime(yrdoy,'%Y%j',tz='MST'))
layername=paste0('X',strftime(dte,'%Y%m%d'))
newlayername=strftime(strptime(layername,'X%Y%m%d',tz='MST'),'X%Y%j')
if(as.numeric(strftime(dte,'%m'))>=3){
    recon=stack(paste0('data/recon_',rswe,'/recondata_',yr,'_',rswe,'.nc'))
    layerind=which(names(recon) %in% layername)
    reconsub=recon[[layerind]]
    names(reconsub)=newlayername
    swesp$model.type='RCN'
    reconswe=extract(reconsub,swesp,sp=T,cellnumbers=T)
    names(reconswe)=c(names(swesp),'cell',newlayername)
    } else {
      reconswe=phvswe
      reconswe$model.type='RCN'
      reconswe[[newlayername]]=NA
    }

    if(yr>2003){
    snodas=stack(paste0('data/spatialvalidation/snodas/snodas.survey.nc'))
    layerind=which(names(snodas) %in% layername)
    snodassub=snodas[[layerind]]
    names(snodassub)=newlayername
    swesp$model.type='SNODAS'
    snodasswe=extract(snodassub,swesp,sp=T,cellnumbers=T)
    names(snodasswe)=c(names(swesp),'cell',newlayername)
    snodasswe[[newlayername]]=snodasswe[[newlayername]]/1000
} else {
      snodasswe=phvswe
      snodasswe$model.type='SNODAS'
      snodasswe[[newlayername]]=NA
    }

    mdlswe=as.data.frame(rbind(phvswe,phvrcnswe,reconswe,snodasswe))
    mdlswe.m=melt(mdlswe,c(names(swesp),'cell','x','y'))
    mdlswe.m=mdlswe.m[,-(ncol(mdlswe.m)-1)]
    colnames(mdlswe.m)[ncol(mdlswe.m)]='swe.model'
    mdlswe.m$swe.model[mdlswe.m$swe.model>10]=NA
    return(mdlswe.m)
}


## CSU
## ----
source('data/spatialvalidation/csu.csv2sp.R')
basepath=paste0('diagnostics/rswe_',rswe,'/covrangeidp0.5/fullpreds/',snotelscale,'/netcdf/')
swe$yrdoy=strftime(swe$date,'%Y%j')
swe$yr=strftime(swe$date,'%Y')
swe=spTransform(swe,CRS('+proj=longlat +datum=WGS84'))
coordnames(swe)=c('x','y')
swe=swe[,c('site','swe.obs','date','yrdoy','yr')]
swe.df=as.data.frame(swe)

csuswe=ddply(swe.df,.(yrdoy),get_modelswe,swe)
str(csuswe)


########  urg surveys
urg=read.table('data/spatialvalidation/urg.validation/rghsurveys.2001.2002.NAD83.csv',sep=',',header=T)
coordinates(urg)=~long+lat
proj4string(urg)=CRS('+proj=longlat +datum=NAD83')
urg=spTransform(urg,CRS('+proj=longlat +datum=WGS84'))
coordnames(urg)=c('x','y')
urg$date=as.POSIXct(strptime(as.character(urg$date),'%B-%d-%Y',tz='MST'))
urg$yr=as.numeric(strftime(urg$date,'%Y'))
urg$yrdoy=as.numeric(strftime(urg$date,'%Y%j'))
urg=urg[,c('site','swe.corr.std.cm','date','yrdoy','yr')]#same var and order as csuswe
names(urg)=c('site','swe.obs','date','yrdoy','yr')
urg$swe.obs=urg$swe.obs/100##convert to meters
urg.df=as.data.frame(urg)
urgswe=ddply(urg.df,.(yrdoy),get_modelswe,urg)
str(urgswe)

## GLV
glvdates=read.table('data/spatialvalidation/glv.validation/snow_survey_dates.txt',sep='\n',stringsAsFactor=F)
colnames(glvdates)='date'
glvdates$date=as.POSIXct(glvdates$date,tz='MST')
glvdates$yr=strftime(glvdates$date,'%Y')
glvdates$mth=strftime(glvdates$date,'%m')
#--
glvdensity=read.table('data/spatialvalidation/glv.validation/snow_survey_densities.txt',sep=',',header=F,stringsAsFactor=F)
colnames(glvdensity)=c('date','density')
glvdensity$date=as.POSIXct(glvdensity$date,tz='MST')
glvdensity$yr=strftime(glvdensity$date,'%Y')
#
glvsurvey=merge(glvdates,glvdensity)
glvsurvey=subset(glvsurvey,yr>2000 & yr <= 2012)
glvsurvey=glvsurvey[!is.na(glvsurvey$density),]

glv.df=ddply(glvsurvey,.(yr),function(dF){
  yr=dF$yr
  print(yr)
  fn=paste0('ss',yr)
  glv=readOGR('data/spatialvalidation/glv.validation/survey_points_qc',fn)
  coordnames(glv)=c('x','y')
  glv=spTransform(glv,CRS('+init=epsg:4326'))
  glv$swe.obs=dF[,'density']*as.numeric(glv$snowdepth)/100
  glv=glv[,'swe.obs']
  glv$date=dF$date
  glv$yrdoy=strftime(dF$date,'%Y%j')
  glv$yr=yr
  glv$site='Green Lakes Valley'
  assign(paste0('glv',yr),envir = .GlobalEnv,glv)
  glvsd=as.data.frame(glv)
  return(glvsd)
})
glvlist=llply(ls(pattern='glv2'),function(x) get(x))
glv=do.call(rbind,glvlist)
ind=which(glv$swe.obs<10)
glv=glv[ind,]
glvswe=ddply(glv.df,.(yrdoy),get_modelswe,glv)
str(glvswe)
str(csuswe)
#
swe=rbind(csuswe,urgswe,glvswe)
saveRDS(swe,file=paste0('data/spatialvalidation/surveypointswe_rswe',rswe,'_',snotelscale,'.rds'))


# summ.df=ddply(swe,.(model.type,swe.model,site,date,yr),function(x)   summarise(x,obs.min=min(swe.obs,na.rm=T), obs.max=max(swe.obs,na.rm=T), obs.mean=mean(swe.obs,na.rm=T),swe.model=unique(swe.model)))
# dev.new()
# ggplot(summ.df)+
#     geom_point(aes(x=obs.mean,y=swe.model,colour=model.type),size=4)+
#     geom_errorbarh(aes(x=obs.mean,y=swe.model,xmin=obs.min,xmax=obs.max,colour=model.type),size=1)+
# #    stat_smooth(aes(x=obs.mean,y=swe.model,color=model.type),method='lm', se=F,size=2)+
#     geom_abline()+
#     theme_bw()+
#     xlab('Observed SWE (m)')+
#     ylab('Modeled SWE (m)')+
#     facet_grid(~yr,scales='free')+
#     guides(color=guide_legend(title='Model'),shape=guide_legend(title='Model'))+
#     ggtitle('Snow Surveys')

# ggsave('graphs/rswe_v3.2/scatter_wrangebars_obs_vs_mdlswe.pdf')

#----
str(swe)
swe=subset(swe,as.character(site)!='Togwotee Pass')#no model data
swe$mth=strftime(swe$date,'%m')
docost=function(y,yhat,cost){
if(cost=='mae') { #statistics to minimize
          costfun<-function(y, yhat) mean(abs(yhat-y),na.rm=T)
} else if(cost=='rmse') {
          costfun<-function(y,yhat) sqrt(mean((yhat-y)^2,na.rm=T))
     }
     costfun(y,yhat)
   }
swecells=ddply(swe,.(cell,model.type,site,date,mth,yr),function(dF) return(data.frame(swe.model=unique(dF$swe.model),swe.obs.avg=mean(dF$swe.obs,na.rm=T)) )  )
sweerror=ddply(swecells,.(model.type,site,date,mth,yr),function(dF) with(dF,docost(swe.obs.avg,swe.model,'rmse')))
names(sweerror)=c('source','site','date','mth','yr','objcost')
obs=swe[,c('site','date','yrdoy','yr','mth','swe.obs')]
obs$source='OBS'
mdl=swe[,c('site','date','yrdoy','yr','mth','swe.model','model.type')]
names(obs) <- names(mdl) <- c('site','date','yrdoy','yr','mth','swe','source')
sweboxplot=rbind(obs,mdl)
sweavg=ddply(sweboxplot,.(site,date,mth,yr,source),function(dF) summarise(dF, avg=mean(swe,na.rm=T)))

ggplot()+
  geom_boxplot(data=sweboxplot,aes(x=source,y=swe),outlier.shape=NA)+
  geom_point(data=sweavg,aes(x=source,y=avg,colour=site),size=3)+
  geom_point(data=sweerror,aes(x=source,y=objcost,colour=site),shape=4,size=3)+
  scale_colour_brewer(palette='Set1')+
  guides(colour=guide_legend(override.aes=list(size=5)))+
  theme_bw()+
  facet_grid(yr~mth,scales='free')+
    theme(legend.position=c(.4,1),
            legend.justification=c(0,1),
             legend.key.size=grid::unit(1.5,'lines'),
             legend.background=element_rect(colour='black',size=1))

ggsave(paste0('graphs/rswe_',rswe,'/boxplot_models_vs_obs_facet-yr&mth_snotel',snotelscale,'.pdf'))

ggplot()+
  geom_boxplot(data=sweboxplot,aes(x=source,y=swe),outlier.shape=NA)+
  geom_point(data=sweavg,aes(x=source,y=avg,colour=site),size=3)+
  geom_point(data=sweerror,aes(x=source,y=objcost,colour=site),shape=4,size=3)+
  scale_colour_brewer(palette='Paired')+
  labs(x='',y='SWE (m)',title='')+
  #scale_y_continuous(lim=c(0,1.5))+
  guides(colour=guide_legend(override.aes=list(size=5),nrow=1))+
  theme_minimal()+
  facet_wrap(~date+site,scales='free_y')+
  theme(legend.position='bottom',
            # legend.justification=c(1,0),
            strip.text=element_text(face='bold'),
            axis.line=element_line(size=.5),
            panel.border=element_rect(fill=NA,colour='grey50'))

ggsave(paste0('graphs/rswe_',rswe,'/boxplot_models_vs_pointobs_facet-datesite_snotel',snotelscale,'.pdf'))