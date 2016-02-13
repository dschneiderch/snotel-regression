# setwd('/Users/dosc3612/Documents/snotel-regression_project')
library(ProjectTemplate)
reload.project(list(cache_loading=F))

library(raster)
library(dplyr)
library(tidyr)
library(ggplot2)
library(doMC)
library(ncdf4)
library(RColorBrewer)
library(broom)
library(rasterVis)

NCARv='1'
product='phvrcn'
recon.version='v3.1'
covrange='idp1'
cost='r2'
residblend=''
scalesnotel='scale'
fscaMatch='wofsca'
dateflag='B'
style='real-time'#nalysis'
gridding='_modisgrid'#'_snodasgrid' or '_modisgrid'
postscaled=''#'' or '_postscaled'
predictor='rcn' #'fsca' or 'rcn'
config=''

rgeo=raster('data/NCAR_WRF/headwaters_grid/Geo_Info_WRF_Headwaters/hgt_m1.tif')

mae=function(yhat, yobs){
     abs(yhat-yobs)
}

## original coord info. dont use ----
nccoordid=nc_open('data/NCAR_WRF/original_runs/geographic_data.nc')
lat=ncvar_get(nccoordid,'XLAT')
long=ncvar_get(nccoordid,'XLONG')
# coords=data_frame(long=as.numeric(long[,,1]),lat=as.numeric(lat[,,1]))

wrfcorner=c(min(long),max(long),min(lat),max(lat))
wrfpoly=as(extent(wrfcorner),'SpatialPolygons')
proj4string(wrfpoly)='+proj=longlat +ellps=sphere +a=6370000 +b=6370000'
wrfpolylcc=spTransform(wrfpoly,CRS('+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39.0000038146973 +lon_0=-107 +ellps=sphere +a=6370000 +b=6370000 +units=m +no_defs'))
plot(wrfpoly)
plot(wrfpolylcc)
(xmin(extent(wrfpolylcc))-xmax(extent(wrfpolylcc)))/318
xmin(extent(wrfpolylcc)) + 4000*317
xmax(extent(wrfpolylcc))

extent(wrflcc)
wrfcds=data.frame(long=as.numeric(long[,,1]),lat=as.numeric(lat[,,1]))
wrfgrid=SpatialPoints(wrfcds)
proj4string(wrfgrid)='+proj=longlat +ellps=sphere +a=6370000 +b=6370000'
wrflcc=spTransform(wrfgrid,CRS('+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39.0000038146973 +lon_0=-107 +ellps=sphere +a=6370000 +b=6370000 +units=m +no_defs'))
wrflcc$swe=as.numeric(dt2[,,417])
SpatialGridDataFrame(wrflcc,dt2[,,417])

tmp=data.frame(x=c(min(long),min(long),max(long),max(long)),y=c(min(lat),max(lat),max(lat),min(lat)))
p4s='+proj=longlat +ellps=sphere +a=6370000 +b=6370000'
p4crs=CRS('+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39.0000038146973 +lon_0=-107 +ellps=sphere +a=6370000 +b=6370000 +units=m +no_defs')
tmppts=SpatialPoints(tmp)
proj4string(tmppts)=p4s
tmppts.lcc=spTransform(tmppts,p4crs)
coords=coordinates(tmppts.lcc)

## small example ----
tmp=data.frame(x=c(min(long),min(long),min(long)+5,max(long),max(long)),y=c(min(lat),max(lat),min(lat)+5,max(lat),min(lat)))
extent(tmp)
p4s='+proj=longlat +ellps=sphere +a=6370000 +b=6370000'
p4crs=CRS('+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39.0000038146973 +lon_0=-107 +ellps=sphere +a=6370000 +b=6370000 +units=m +no_defs')
tmppoly=as(extent(tmp),'SpatialPolygons')
tmppts=SpatialPoints(tmp)
proj4string(tmppoly)=p4s
proj4string(tmppts)=p4s
plot(tmppoly)
plot(tmppts,add=T,col='red')
tmppts.lcc=spTransform(tmppts,p4crs)
tmppoly.lcc=spTransform(tmppoly,p4crs)
plot(tmppoly.lcc)
plot(tmppts.lcc,add=T,col='red')

## get wrf swe data ----
ncswe=nc_open('data/NCAR_WRF/headwaters_grid/SWE_daily.nc')
ncswe
data_temp=ncvar_get(ncswe,'SNOW')/917#kg/m^3 #convert to meters from kg/m^2
dt2=aperm(data_temp,c(2,1,3))
dt2=dt2[nrow(dt2):1,,]
rb=brick(dt2)
extent(rb)=extent(rgeo)
crs(rb)=crs(rgeo)
startdate=as.Date('2000-10-01',tz='MST')
ncdates=startdate+0:(nlayers(rb)-1)
names(rb)=ncdates
setZ(rb,ncdates,'time')


## reproject phvrcn and write to disk ----
yr=2001
for(yr in 2001:2008){
     print(yr)
     predfile=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',scalesnotel,'/',dateflag,'/fullpreds/',cost,'/netcdf/',style,'/',config,'/',residblend,'preds-',product,'_',yr,'_blend-',fscaMatch,'.nc')
     phvrcn=brick(predfile)
     layernames=names(phvrcn)
     phvrcn[phvrcn>200]=NA
     phvrcn=aggregate(phvrcn,c(8,8),fun=mean,na.rm=TRUE)
     names(phvrcn)=layernames
     phvrcn=setZ(phvrcn,layernames,name = 'time')
     phvrcn.lcc=projectRaster(phvrcn,rgeo,method='bilinear')
     # phvrcn.lcc=brick(paste0('data/NCAR_WRF/headwaters_grid/phvrcn-',yr,'-wrfgrid.tif'))
     names(phvrcn.lcc)=layernames
     writeRaster(phvrcn.lcc,paste0('data/NCAR_WRF/headwaters_grid/phvrcn-',yr,'-wrfgrid.tif'),overwrite=TRUE)
     #
}

## difference phvrcn and wrf swe as rasters and save spatial maps. calc bias too. ----
biasdf=data_frame()
yr=2001
for(yr in 2001:2008){
     predfile=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',scalesnotel,'/',dateflag,'/fullpreds/',cost,'/netcdf/',style,'/',config,'/',residblend,'preds-',product,'_',yr,'_blend-',fscaMatch,'.nc')
     phvrcn=brick(predfile)
     layernames=names(phvrcn)
     phvrcn.lcc=brick(paste0('data/NCAR_WRF/headwaters_grid/phvrcn-',yr,'-wrfgrid.tif'))     
     names(phvrcn.lcc)=layernames
     #
     wrfformat=strftime(as.Date(layernames,'X%Y%m%d',tz='MST'),'X%Y.%m.%d')
     wrflayers=which(names(rb) %in% wrfformat)
     wrfraster=rb[[wrflayers]]
     #
     if(nlayers(wrfraster)==nlayers(phvrcn.lcc)){
          swediffraster=wrfraster-phvrcn.lcc
     } else stop()
     
     biasdf=bind_rows(biasdf,
          data.frame(avgdiff=cellStats(swediffraster,mean)) %>%
               mutate(dte=rownames(.))
     )
     
## plot spatial map of differences
     # i=5
     # for(i in 1:nlayers(swediffraster)){
     #      zrange=seq(-.525,.525,by=0.05)#should be an even number of colors
     #      numcols=length(zrange)+1
     #      mypal=colorRampPalette(c('red','white','blue'))(numcols)
     #      mypala=rev(brewer.pal(length(zrange)/2,'Blues'))
     #      mypalb=brewer.pal(length(zrange)/2,'Reds')
     #      # mypal=c(mypala,'grey99',mypalb)
     #      myTheme <- BTCTheme()
     #      myTheme$panel.background$col='grey50'
     #      png(filename=paste0('data/NCAR_WRF/headwaters_grid/graphs/swediff_rasters/swediff_raster_',layernames[i],'.png'),res=600,height=6,width=6,units = 'in')
     #           print(
     #                levelplot(swediffraster[[i]],main=names(swediffraster)[i],par.settings=myTheme,col.regions=mypal, at=zrange)
     #           )
     #                dev.off()
     # }
}

## plot bias for each year
biasdf=biasdf %>%
     mutate(dte=as.Date(strptime(dte,'X%Y.%m.%d')),
            yr=strftime(dte,'%Y'),
            doy=strftime(dte,'%j'),
            dymnth=strftime(dte,'%d%b'))

avgbias=biasdf %>%
     group_by(yr) %>%
     summarise(avg=sprintf(fmt='%.3f',mean(avgdiff))) %>%
     mutate(y=seq(-0.05,-0.12,length.out=length(unique(.$yr))),
            x=as.Date('2012-06-01'))

ggplot(biasdf,aes(x=as.Date(paste(2012,strftime(dte,format="%m-%d"),sep='-')),y=avgdiff,colour=yr))+
     labs(x='date',y='bias [m]')+
     geom_point(size=2)+
     geom_smooth(method='loess',se=FALSE)+
     geom_text(data=avgbias,aes(x,y,label=avg,colour=yr),size=4)+
     scale_x_date(date_breaks='month',date_labels='%d-%B')+
     scale_colour_brewer(palette='Set1')+
     theme_bw(base_size=14)+
     theme(panel.background=element_rect(fill='grey20'),
           panel.grid=element_line(colour='grey30'),
           axis.text.x=element_text(angle=45,hjust=1))

ggsave('data/NCAR_WRF/headwaters_grid/graphs/bias_wrf-phvrcn_coloryears.png',dpi=600,height=6,width=6,units='in')



## combine all reprojected phvrcn in dataframe and save dates----
# writeRaster isn't writing the layernames properly so read in original phvrcn and use those names
phvrcndates=data_frame()
phvrcndata=data_frame()
for(yr in 2001:2008){
     predfile=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',scalesnotel,'/',dateflag,'/fullpreds/',cost,'/netcdf/',style,'/',config,'/',residblend,'preds-',product,'_',yr,'_blend-',fscaMatch,'.nc')
     phvrcn=brick(predfile)
     layernames=names(phvrcn)
     phvrcn.lcc=brick(paste0('data/NCAR_WRF/headwaters_grid/phvrcn-',yr,'-wrfgrid.tif'))     
     names(phvrcn.lcc)=layernames
     #
     phvrcndata=as.data.frame(phvrcn.lcc,xy=TRUE) %>%
          tbl_df %>%
          gather(dte,swe,-x,-y) %>%
          mutate(dte=as.Date(dte,'X%Y%m%d',tz='MST'),
                 model='phvrcn') %>%
          bind_rows(phvrcndata)

     phvrcndates=bind_rows(phvrcndates,data_frame(wrfformat=strftime(as.Date(layernames,'X%Y%m%d',tz='MST'),'X%Y.%m.%d')))
}     

## extract days from wrf that match phvrcn. convert to dataframe
wrflayers=which(names(rb) %in% phvrcndates$wrfformat)
wrfdata=as.data.frame(rb[[wrflayers]],xy=TRUE)
wrfdata=gather(wrfdata,dte,swe,-x,-y) %>%
     mutate(dte=as.Date(dte,'X%Y.%m.%d',tz='MST'),
            model='wrf'
     )

## wrf elevation data
wrfelev=as.data.frame(rgeo,xy=TRUE) 

## join phvrcn, wrf and wrfelev
alldata=bind_rows(wrfdata,phvrcndata) %>% 
     full_join(wrfelev)

## substract rasters


## extract 01Apr from alldata for plotting. can't spread all dates (memory issues) ----
alldata2=alldata %>%
     mutate(dymnth=strftime(dte,'%d%b'),
            yr=strftime(dte,'%Y')) %>%
     filter(dymnth=='01Apr') %>%
     spread(model,swe)

## fit trend line for each year between wrf and phvrcn swe
eqndf=alldata2 %>%
     group_by(yr) %>%
     do(
          data_frame(eqn=lm_eqn(glm(data=.,wrf~phvrcn,family=gaussian(link='identity'))))
     )

## plot phvrcn vs wrf by year
ggplot(alldata2,aes(x=phvrcn,y=wrf))+
     geom_point()+
     geom_abline()+
     geom_smooth(method='lm',colour='red')+
     geom_text(data=eqndf,aes(x=0,y=1.5,label=eqn),hjust=0,parse=TRUE,colour='red')+
     coord_equal()+
     facet_wrap(~yr)+
     ggtitle('April 1')

ggsave('data/NCAR_WRF/headwaters_grid/graphs/phvrcn-wrf-facetyrs-01Apr.png',dpi=600,width=12,height=8,units='in')

## plot phvrcn - wrf difference by year and by elevation
alldata2 %>%
     mutate(swediff=wrf-phvrcn) %>%
     {
          ggplot(.)+
               geom_point(aes(x=hgt_m1,y=swediff))+
               geom_smooth(aes(x=hgt_m1,y=swediff),method='lm')+
               facet_wrap(~yr)+
               ggtitle('swe difference\nwrf-phvrcn\npositive difference=wrf greater')
     }

ggsave('data/NCAR_WRF/headwaters_grid/graphs/elev-phvrcn_wrf_diff-facetyrs-01Apr.png',dpi=600,width=12,height=8,units='in')

## trendlines of swediff ~ elev are significantly different than 0 ----
alldata2 %>%
     mutate(swediff=wrf-phvrcn) %>% 
     group_by(yr) %>%
     do(
          tidy(anova(lm(data=.,swediff~hgt_m1)))
     )

## 

