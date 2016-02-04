library('ProjectTemplate')
# setwd('/Volumes/Dominik/Documents/snotel-regression_project')
load.project()
library(dplyr)
library(tidyr)

dem.m=raster('data/gis/gmted_uco_recongrid.tif')

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

demelev=dem.m %>% as.data.frame %>% tbl_df %>%
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

ggsave('graphs/elevation_distribution.png',dpi=300)


ggplot()+
     geom_freqpoly(data=filter(elevdf,source=='DEM'),aes(x=elev,colour=source),binwidth=500,size=2)+
     geom_dotplot(data=filter(elevdf,source=='SNOTEL' | source=='SURVEY'),aes(x=elev,fill=source),binaxis='x',binwidth=250,method="dotdensity")
     # guides(colour=guide_legend(''))+
     labs(x='Elevation [m]')+
     theme_bw(base_size = 16)+
     theme(legend.title=element_blank(),
           legend.key=element_rect(colour=NA),
           legend.key.size = unit(.05, "npc"),
           legend.position=c(1,1),
           legend.justification=c(1,1))



