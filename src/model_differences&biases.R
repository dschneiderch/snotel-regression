library('ProjectTemplate')
# setwd('/Volumes/Dominik/Documents/snotel-regression_project')
load.project()
library(dplyr)

recon.version='v3.1'
cost='r2'
style='real-time'
covrange='idp1'#'300km-20150313Fri'
scalesnotel='scale'
fscaMatch='wofsca'
dateflag='B'
residblending='unblended'

swediffmap=point_plot_setup(recon.version,style,cost,covrange,scalesnotel,fscaMatch,dateflag)

swediffmap2=swediffmap %>%
     select(Station_ID,date,mnth,yrdoy,phv.pred,phvrcn.pred,swediff.phv,swediff.phvrcn,fsca,yr,snotel,swe,elev) %>%     mutate(
          phv.phvrcn.diff=phv.pred-phvrcn.pred)
     
ggplot(swediffmap2)+geom_point(aes(x=elev,y=phv.phvrcn.diff))

ggplot(swediffmap2)+
     geom_point(aes(x=swediff.phv,y=swediff.phvrcn,colour=mnth),alpha=0.5)+
     geom_abline()+
     scale_colour_brewer(palette='Set1')+
     coord_equal()

ggplot(swediffmap2)+
     geom_bin2d(aes(x=swe,y=phv.pred),binwidth=0.1)+
     geom_abline()

ggplot(swediffmap2)+
     geom_bin2d(aes(x=snotel,y=phvrcn.pred),binwidth=0.1)+
     geom_abline()
