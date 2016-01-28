## uco_basin_vol.R
library('ProjectTemplate')
# setwd('/Volumes/Dominik/Documents/snotel-regression_project')
load.project()
library(doMC)
library(raster) 
# library(rgdal)
# library(plyr)
# library(reshape2)
# library(gtable)
# library(gridExtra)

dem.m=raster('data/gis/gmted_uco_recongrid.tif')
demm=dem.m
demm[demm<2000]=1000
demm[demm>=2000 & demm <2500]=2000
demm[demm>=2500 & demm<3000]=2500
demm[demm>=3000 & demm<3500]=3000
demm[demm>=3500 & demm<4000]=3500
demm[demm>=4000 & demm<4500]=4000

writeRaster(demm,'data/gis/gmted_uco_recongrid_classified_metric.tif',overwrite=T)

elev_area=freq(demm)
elev_area=mutate(as.data.frame(elev_area),
                 area=count*500*500,
                 areaper=count/sum(count))#

recon.version='v3.1'
covrange='idp1'
cost='r2'
product='phvrcn'
fscaMatch='wofsca'
snotelscale='scale'
residblending='unblended'
if(residblending=='unblended') {
     resid='' 
} else {
     resid='full'
}
config=''#bic-v2.5-removed_NWbdiff'#'bic-v2.2-nosnoteltransform'
dateflag='B'
style='real-time'
dte=as.POSIXct('1900-04-01',tz='MST')

graphbase=paste0('graphs/rswe_',recon.version,'/covrange',covrange,'/snotel',snotelscale,'/',dateflag,'/',fscaMatch,'/',style)

units='metric'
yr=2012
mth='04'
dy='01'

# huc4=huc4[!(huc4$HU_4_Name %in% c('Rio Grande Closed Basins','Upper Canadian','Upper Gila','Upper Pecos')),]

## create dataframe of all data
dte=as.POSIXct('1900-04-01',tz='MST')
fordenstack=stack('data/forestdensity/umd_forden.nc')
huc4raster=raster('data/gis/UpperCRB_rasterized.tif')
huc4=readOGR('data/gis','UpperCRB')
#
basinsize=as.data.frame(freq(huc4raster))
#
basinelev=as.data.frame(zonal(dem.m,huc4raster,'mean'))
colnames(basinelev)=c('zone','elevavg')
basinelev$elevavg=floor(basinelev$elevavg)
#
yr=2001
registerDoMC(3)
swestats=ldply(seq(2001,2012),.parallel=T,function(yr){
     print(yr)
     mth='04'
     dy='01'
     
     layerid=paste0('X',yr-2011+12)
     layernum=grep(layerid,names(fordenstack))
     forden=fordenstack[[layernum]]
     basinforden=as.data.frame(zonal(forden,huc4raster,'mean'))
     basinforden=merge(as.data.frame(huc4[,c('HUC_4','HU_4_Name')]),basinforden,by.x='HUC_4',by.y='zone')
     colnames(basinforden)=c('zone','basin','fordenavg')
     basinforden$fordenavg=floor(basinforden$fordenavg*100)/100
     
     basin_static=merge(basinelev,basinforden,by='zone')
     if(snotelscale=='scale'){
          fscastack=stack(paste0('data/recon_',recon.version,'/reconfscadata_',yr,'_',recon.version,'.nc'))
          band=grep(paste0(yr,mth,dy),names(fscastack))
          fsca=fscastack[[band]]
     }
     
     swetmp=data.frame()
     for(product in c('phv','phvrcn')){#			
          predfile=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',snotelscale,'/',dateflag,'/fullpreds/',cost,'/netcdf/',style,'/',config,'/',resid,'preds-',product,'_',yr,'_blend-',fscaMatch,'.nc')
          s=stack(predfile)
          band=grep(paste0(yr,mth,dy),names(s))
          r=s[[band]]
          r[r==253]=NA
          if(snotelscale=='scale'){
               r=r*fsca
          }
          #
          r.ext=as.data.frame(zonal(r,huc4raster,'mean'))#swe avg is meters
          colnames(r.ext)=c('zone','sweavg')
          r.ext$vol=r.ext$sweavg*basinsize$count*500*500/10^9#cu km
          if(units=='metric'){
               basin_stats=merge(basin_static,r.ext,by='zone')
               swedF=data.frame(date=as.POSIXct(paste0(yr,'-',mth,'-',dy),tz='MST'),yr=yr,model=product,basin_stats)
          }
          swetmp=rbind(swetmp,swedF)
     }
     return(swetmp)
})

saveRDS(swestats,file=paste0(graphbase,'/huc4_',residblending,'swestats_2001-2012.rds'))

swestats=readRDS(paste0(graphbase,'/huc4_',residblending,'swestats_2001-2012.rds'))

##Section 5.5- Basin-wide volume differences
### paragraph 1
yrlyvol=dcast(swestats,yr~model,value.var='vol',fun.agg=sum)
yrlyvol=mutate(yrlyvol,
               diff=phvrcn-phv,
               diffpct=diff/phv*100)
arrange(yrlyvol,diffpct)

cor(yrlyvol$diffpct,yrlyvol$phv)

### paragraph 2
dcast(swestats,basin+yr~model,value.var='sweavg')
basinvol=dcast(swestats,basin+yr~model,value.var='vol')
basinvol=mutate(basinvol,
                diff=phvrcn-phv,
                diffpct=diff/phv*100)
basindiff=dcast(basinvol,basin~.,value.var='diffpct',fun.agg=mean)
colnames(basindiff)=c('basin','avgdiffpct')
basindiff$gtmodel=ifelse(basindiff[,2]>0,'PHVRCN','PHV')
ddply(basindiff,.(gtmodel),summarise,
      avg=mean(avgdiffpct))


# 6.1 comparision with 


##-----------------

swerange=ddply(swestats,.(model,basin),function(dF){
     summarise(dF,
               avg=mean(vol,na.rm=F),
               min=min(vol,na.rm=F),
               max=max(vol,na.rm=F),
               forden=mean(fordenavg,na.rm=F),
               elev=mean(elevavg,na.rm=F))
})
avgdiff=dcast(swerange,basin~model,value.var='avg')


avgdiff
swerange2=ddply(swerange,.(basin,model),function(dF){
     summarise(dF,
               range=max-min,
               avg=avg,
               forden=forden,
               elev=elev)
})
swerange2
ggplot(swerange2)+geom_point(aes(x=forden,y=avg,size=elev))+facet_wrap(~model)
dcast(swerange,basin~model,value.var='avg')


dcast(swerange,model~.,value.var='avg',sum)
basinrangediff=dcast(swerange2,basin~model,value.var='range')
basinrangediff=mutate(basinrangediff,
                      diff=phvrcn-phv,
                      diffpct=diff/phv*100)
mean(basinrangediff$phv)
mean(basinrangediff$phvrcn)
mean(basinrangediff$diffpct)

sum(swerange$avg)
head(swerange)
swerange2=swerange#[!(swerange$basin %in% c('Rio Grande Closed Basins','Upper Canadian','Upper Gila','Upper Pecos')),]

ggplot(subset(swerange,model=='phvrcn'))+geom_bar(aes(x=basin,y=forden),stat='identity')
ggplot(subset(swerange,model=='phvrcn'))+geom_bar(aes(x=basin,y=elev),stat='identity')


dcast(swerange2,model~.,fun.aggregate=mean,na.rm=T,value.var='avg')

swerange$basin=paste0(seq(1,8),'. ',swerange$basin)
ggplot(swerange)+
     geom_bar(aes(x=basin,y=avg,fill=model),stat='identity',position=position_dodge(width=0.9))+
     geom_errorbar(aes(x=basin,ymin=min,ymax=max,color=model),position=position_dodge(width=0.9),width=.5)+
     scale_color_manual(values=c('black','black'),guide=F)+
     scale_fill_manual(values=c('red','blue'),labels=c('PHV','PHV+RCN'))+
     guides(fill=guide_legend('Model'))+
     labs(x='Basin',y='SWE Volume [cu. km]')+
     theme_minimal()+
     theme(
          axis.line=element_line(colour='grey10'),
          axis.text.x=element_text(size=6,angle=45,hjust=1,vjust=1),
          axis.text.y=element_text(size=6),
          axis.title.x=element_text(size=8),
          axis.title.y=element_text(size=8),
          legend.title=element_text(size=8),
          legend.text=element_text(size=6),
          legend.key=element_rect(colour='white'),
          legend.key.size=unit(.5,'lines'),
          legend.position = c(.8,.85) )

ggsave(filename=paste0(graphbase,'/huc4_',residblending,'swe.pdf'),width=6.6,height=3)

write.table(format(swestats,digits=2,nsmall=0,scientific=F),paste0('reports/',snotelscale,'/huc4_',resid,'swetable_',units,'_',dy,mnth,yr,'-',product,'-',covrange,'-',cost,'-',fscaMatch,'.txt'),sep='\t',quote=F,row.names=F)
swestats_all=rbind(swestats_all,swestats)




