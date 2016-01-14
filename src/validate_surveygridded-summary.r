# setwd('/Volumes/Dominik/Documents/snotel-regression_project')
setwd('~/Documents/snotel-regression_project')
library(raster)
library(reshape2)
library(plyr)
library(ggplot2)
library(doMC)
library(RColorBrewer)
library(dplyr)
library(tidyr)

## define cost function
docost=function(y,yhat,cost){
if(cost=='mae') { #statistics to minimize
          costfun<-function(y, yhat) mean(abs(yhat-y),na.rm=T)
} else if(cost=='rmse') {
          costfun<-function(y,yhat) sqrt(mean((yhat-y)^2,na.rm=T))
} else if (cost=='relrmse'){
          costfun<-function(y,yhat) sqrt(mean((yhat-y)^2,na.rm=T))/mean(y,na.rm=T)
} else if (cost=='bias'){
          costfun<-function(y,yhat) mean(yhat-y,na.rm=T)
} else if (cost=='pctbias'){
              costfun<-function(y,yhat) mean(yhat-y,na.rm=T)/mean(y,na.rm=T)*100
          }
              
     costfun(y,yhat)
}

allcost=function(surveytype,cost){
  if(surveytype=='alpine') swe=filter(cellavg,site=='Green Lakes Valley')
  if(surveytype=='forested') swe=filter(cellavg,site!='Green Lakes Valley')
  if(postscaled==''){
    postscaled='notpostscaled'
  } else {
    postscaled='postscaled'
  }
  if(gridding==''){
    gridding='modis'
  } else {
    gridding='snodas'
  }

  swe %>%
    group_by(model.type,date,site,yr,mth) %>%
    do(
     data.frame(cost=with(.,docost(swe.obs.avg,swe.model,cost)))
     ) %>%
    mutate(
      'surveytype'=surveytype,
      'scalesnotel'=scalesnotel,
      'postscaled'=postscaled,
      'fscaMatch'=fscaMatch,
      'blending'=blending,
      'style'=style,
      'gridding'=gridding
      )
}

recon.version='v3.1'
covrange='idp1'
cost='r2'
dateflag='surveyvalidation'

fscaMatch='wofsca'
scalesnotel='scale'
blending='unblended'
style='real-time'
gridding='_snodasgrid'#_snodasgrid'
postscaled=''#_postscaled'

objcost='rmse'
for(objcost in c('rmse','pctbias')){
  erralldF=data_frame()
  errdF=data_frame()
  errtypedF=data_frame()
  print(objcost)

for(cost in 'r2'){
  for(blending in c('unblended')){
    for(style in 'real-time'){
      for(gridding in c('','_snodasgrid')){
        for(fscaMatch in c('wofsca','fsca')){
          for(scalesnotel in c('noscale','scale')){
            for(postscaled in c('','_postscaled')){
if(scalesnotel=='scale' & postscaled=='_postscaled') next
if(dateflag=='surveyvalidation' & gridding=='_snodasgrid') next
basepath=paste0('data/spatialvalidation/',blending,'/',dateflag,'/',style,'/')
graphpath=paste0('graphs/rswe_',recon.version,'/covrange',covrange,'/snotel',scalesnotel,'/',dateflag,'/',fscaMatch,'/',style,'/')

swe=readRDS(file=paste0(basepath,'surveygriddedswe',gridding,'_rswe',recon.version,'_',covrange,'_snotel',scalesnotel,postscaled,'_',cost,'_',fscaMatch,'.rds'))
swe_all=swe
if(gridding!='_snodasgrid'){
cellcount=ddply(swe_all,.(model.type),function(dF) plyr::count(dF$cell))
freqdf=dcast(cellcount,model.type~x,value.var='freq')
ind=which(freqdf[1,2:ncol(freqdf)]<140)#about 278 30m cells in an 500m pixel #rows 1,2,3 will be the same
cell2rm=as.numeric(colnames(freqdf)[-1][ind])
ind2=which(freqdf[nrow(freqdf),2:ncol(freqdf)]<560)#about 1111 30m cells in 1km pixel #last row is snodas
cell2rm2=as.numeric(colnames(freqdf)[-1][ind2])
swe=subset(swe_all,(model.type!='SNODAS' & !(cell %in% cell2rm)) | (model.type=='SNODAS' & !(cell %in% cell2rm2)))
}
swe$swe.model[swe$swe.model<0.001]=0

cellavg=ddply(swe,.(model.type,date,site,yr,mth,cell),function(dF){ 
    summarise(dF,
            swe.obs.avg=mean(swe.obs,na.rm=T),#the cell number in swe represents the model cell number so there can be muliple obs cell values (30m)
                                        # na.rm=F will average only for model pixels that are completely covered by observations, otherwise there can be overlap
              swe.model=unique(swe.model))# there should only be one value for each cell
    })

if(dateflag=='B2'){
  cellavg=filter(cellavg,model.type=='PHV' | model.type=='PHVRCN' | model.type=='SNODAS')  
} else if(dateflag=='surveyvalidation') {
  cellavg=filter(cellavg,model.type=='PHV' | model.type=='PHVRCN')  
}

## Green Lakes
err_glv=allcost('alpine',objcost)

## Forested
err_for=allcost('forested',objcost)

## All Errors
err_all=rbind(err_glv,err_for) 
err_all=arrange(err_all,model.type,surveytype,date)
fn=paste0('boxplot_avg',objcost,'_',blending,'models',postscaled,gridding,'_vs_griddedobs_facet-datesite_',cost,'-rcnselect.txt')
write.table(err_all,paste0(graphpath,fn),sep='\t',row.names=F,quote=F)

erralldF=rbind(erralldF,err_all)

## Mean error all
avgerr_all=err_all %>%
  group_by(model.type,scalesnotel,postscaled,fscaMatch,blending,style,gridding) %>%
    summarise(
      err=mean(cost,na.rm=T)
        ) 

errdF=rbind(errdF,avgerr_all)

## Mean Error by survey type
avgerr_type=err_all %>%
  group_by(model.type,surveytype,scalesnotel,postscaled,fscaMatch,blending,style,gridding) %>%
    summarise(
      err=mean(cost,na.rm=T)
        )
errtypedF=rbind(errtypedF,avgerr_type)

}}}}}}}

if(dateflag=='B2'){
## count best model
if(objcost=='rmse'){
diffdF=erralldF %>% 
  spread(model.type,cost) %>%
  filter(yr>2003) %>%
  mutate(
    phvrcndiff=PHVRCN-SNODAS,
    phvdiff=PHV-SNODAS,
    phvrcnbest=ifelse(phvrcndiff<=0,'PHVRCN','SNODAS'),
    phvbest=ifelse(phvdiff<=0,'PHV','SNODAS')
    )   
} else if(objcost=='pctbias'){
diffdF=erralldF %>% 
  spread(model.type,cost) %>%
  filter(yr>2003) %>%
  mutate(
    phvrcndiff=abs(PHVRCN)-abs(SNODAS),
    phvdiff=abs(PHV)-abs(SNODAS),
    phvrcnbest=ifelse(phvrcndiff<=0,'PHVRCN','SNODAS'),
    phvbest=ifelse(phvdiff<=0,'PHV','SNODAS')
    )     
}
diffdF  %>%
  group_by(scalesnotel,postscaled,fscaMatch,blending,style,gridding,phvrcnbest) %>% 
  tally %>%
  write.table(file=paste0(basepath,'phvrcn_bestfreq_',objcost,gridding,'_',cost,'-rcnselect.txt'),row.names=F,quote=F,sep='\t')

diffdF  %>%
  group_by(scalesnotel,postscaled,fscaMatch,blending,style,gridding,phvbest) %>% 
  tally %>%
  write.table(file=paste0(basepath,'phv_bestfreq_',objcost,gridding,'_',cost,'-rcnselect.txt'),row.names=F,quote=F,sep='\t') 
}

avgobjerr_all = errdF %>% spread(model.type,err)
fn=paste0('avg',objcost,'_',blending,'models',gridding,'_vs_griddedobs_facet-datesite_',cost,'-rcnselect.txt')
write.table(format(avgobjerr_all,digits=2),file.path(basepath,fn),sep='\t',row.names=F,quote=F)

avgobjerr_type = errtypedF %>% spread(model.type,err)
fn=paste0('avg',objcost,'_surveytype_',blending,'models',gridding,'_vs_griddedobs_facet-datesite_',cost,'-rcnselect.txt')
write.table(format(avgobjerr_type,digits=2),file.path(basepath,fn),sep='\t',row.names=F,quote=F)
}

