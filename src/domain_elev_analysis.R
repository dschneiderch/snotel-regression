# setwd('/Users/dosc3612/Documents/snotel-regression_project')
library(ProjectTemplate)
load.project(list(cache_loading=TRUE))

library(raster)
library(dplyr)
library(tidyr)
library(ggplot2)
library(doMC)
library(ncdf4)
library(RColorBrewer)


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

yr=2001
mthdy='0401'
for(mthdy in c('0401')){
for(thresh in c(0.01,0.05,0.1,0.15)){
s1=stack()
for(yr in seq(2001,2012)){
     # for(mth in strftime(dte,'%m')){
     # month=strftime(dte,'%B')
     # for(dy in strftime(dte,'%d')){
     predfile=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',scalesnotel,'/',dateflag,'/fullpreds/',cost,'/netcdf/',style,'/',config,'/',residblend,'preds-',product,'_',yr,'_blend-',fscaMatch,'.nc')
     s=brick(predfile)
     dateind=grep(paste0(yr,mthdy),names(s))
     s1=addLayer(s1,s[[dateind]])
}

s1
s1[s1>200]=NA
s1avg=mean(s1,na.rm=T)
# thresh=0.1
s1avg.masked=s1avg
s1avg.masked[s1avg.masked<thresh]=NA

# plot(s1avg.masked)

dem.m=raster('data/gis/gmted_uco_recongrid.tif')

dem.m.sca=mask(dem.m,s1avg.masked,maskvalue=NA)
# plot(dem.m.sca)



glvsurvey=data_frame(Site_Name='Green Lakes Valley',Latitude=40.0516,Longitude=-105.6344,Elevation_m=as.integer(3636))
#
surveysites=c('SLUMGULLION','UPPER SAN JUAN', 'WOLF CREEK SUMMIT','UPPER RIO GRANDE','NIWOT','JOE WRIGHT','LIZARD HEAD','LILY POND','DRY LAKE','SOUTH BRUSH CREEK')
surveylocs=snotellocs %>% as.data.frame %>% tbl_df %>%
     filter(Site_Name %in% surveysites)  %>%
     select(Site_Name,Latitude,Longitude,Elevation_m) %>%
     full_join(glvsurvey) %>%
     rename(elev=Elevation_m)

surveylocs.df=surveylocs %>% 
     select(elev) %>%
     mutate(source='SURVEY')

snotellocs.df=snotellocs %>% as.data.frame %>% tbl_df %>%
     select(Elevation_m) %>%
     mutate(source='SNOTEL',
            elev=as.numeric(Elevation_m)) %>%
     select(-Elevation_m)

demelev=dem.m.sca %>% as.data.frame %>% tbl_df %>%
     mutate(source='DEM')
colnames(demelev)=c('elev','source')

elevdf=bind_rows(snotellocs.df,demelev,surveylocs.df)

ggplot(elevdf)+
     geom_freqpoly(aes(x=elev,colour=source,y=..density..),binwidth=250,size=2)+
     # guides(colour=guide_legend(''))+
     labs(x='Elevation [m]')+
     theme_bw(base_size = 16)+
     theme(legend.title=element_blank(),
           legend.key=element_rect(colour=NA),
           legend.key.size = unit(.05, "npc"),
           legend.position=c(1,1),
           legend.justification=c(1,1))

ggsave(paste0('graphs/elevation_distribution-',mthdy,'-',thresh,'.png'),dpi=300)



}

}