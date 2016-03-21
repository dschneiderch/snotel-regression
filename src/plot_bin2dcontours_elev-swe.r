library('ProjectTemplate')
# setwd('/Volumes/Dominik/Documents/snotel-regression_project')
reload.project(override.config = list(cache_loading=F))
# library(doMC)
# library(scales)
library(dplyr)
library(tidyr)

dem.m=raster('data/gis/gmted_uco_recongrid.tif')

# ----Read in some gis files
ucomask=raster('data/gis/UpperCRB_mask.tif')
dem.m.uco=mask(dem.m,ucomask)

dem.m.uco.df=as.data.frame(dem.m.uco,xy=FALSE) %>% 
     tbl_df %>%
     rename(elev=gmted_uco_recongrid) %>%
     mutate(id=1:nrow(.)) %>%
     filter(!is.na(elev))

units='metric'
#--- Regression Analysis ---
recon.version='v3.1'
covrange='idp1'
scalesnotel='scale'
product='phvrcn'
cost='r2'
fscaMatch='wofsca'
style='real-time'
residblending='unblended'
if(residblending=='blended'){
     resid='full'
} else {
     resid=''
}
dateflag='B'
config=''
graphbase=paste0('graphs/rswe_',recon.version,'/covrange',covrange,'/snotel',scalesnotel,'/',dateflag,'/',fscaMatch,'/',style,'/')


yr=2010
dte=as.POSIXct('1900-03-01',tz='MST')
dy=strftime(dte,'%d')
month=strftime(dte,'%B')
mnth=strftime(dte,'%b')
mth=strftime(dte,'%m')
dtestack=stack()
for(product in c('phv','phvrcn')){
     for (yr in seq(2011,2012) ) {#yr=2010#2011#2012
          predfile=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',scalesnotel,'/',dateflag,'/fullpreds/',cost,'/netcdf/',style,'/',config,'/',resid,'preds-',product,'_',yr,'_blend-',fscaMatch,'.nc')
          s=stack(predfile)
          band=grep(paste0(yr,mth,dy),names(s))
          r=s[[band]]
          r[r==253]=NA
          r=mask(r,ucomask)
          names(r)=paste(product,names(r),sep='.')
          # plot(r)
          # projection(r)='+init=epsg:4326'
          if(scalesnotel=='noscale'){
               fscastack=stack(paste0('data/recon_',recon.version,'/reconfscadata_',yr,'_',recon.version,'.nc'))
               band=grep(paste0(yr,mth,dy),names(fscastack))
               fsca=fscastack[[band]]
               r=r*fsca
          }
          dtestack=addLayer(dtestack,r)
     }
}


swedat=as.data.frame(dtestack,xy=TRUE) %>% 
     tbl_df %>%
     mutate(id=1:nrow(.)) %>%
     inner_join(.,dem.m.uco.df,by=c('id')) %>%
     # rename(x=x.x,y=y.x) %>%
     # select(-x.y,-y.y) %>%
     gather(model,swe,-elev,-id) %>%
     separate(col=model,into=c('model','dte'),sep='.X')
  
saveRDS(swedat,paste0(graphbase,'March01_SWEdepth_elev.rds'))

swedat %>%
     filter(dte=='20120301' | dte=='20110301') %>% 
     {
          ggplot(.)+
               geom_density2d(aes(x=elev,y=swe,colour=model),h=c(100,0.1),contour=TRUE)+
               facet_wrap(~dte)
     }

ggsave(paste0(graphbase,'densitycontour_elev_swedepth.png'),dpi = 300)
