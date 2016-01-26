library('ProjectTemplate')
# setwd('/Volumes/Dominik/Documents/snotel-regression_project')
load.project()
library(dplyr)
library(tidyr)



snotellocs.df=snotellocs %>% as.data.frame %>% tbl_df %>%
     select(Elevation_m) %>%
     mutate(source='SNOTEL',
            elev=as.numeric(Elevation_m)) %>%
     select(-Elevation_m)
demelev=dem.m %>% as.data.frame %>% tbl_df %>%
     mutate(source='DEM')
colnames(demelev)=c('elev','source')
elevdf=rbind(snotellocs.df,demelev)

ggplot(elevdf)+
     geom_freqpoly(aes(x=elev,colour=source,y=..density..),binwidth=200)+
     # guides(colour=guide_legend(''))+
     labs(x='Elevation [m]')+
     theme_bw(base_size = 18)+
     theme(legend.title=element_blank(),
           legend.key=element_rect(colour=NA),
           legend.key.size = unit(.1, "cm"),
           legend.position=c(1,1),
           legend.justification=c(1,1))

ggsave('graphs/elevation_distribution.png',dpi=600,width=3.3,height=2,units = 'in')
