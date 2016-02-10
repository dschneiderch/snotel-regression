library('ProjectTemplate')
# setwd('/Volumes/Dominik/Documents/snotel-regression_project')
load.project()
library(dplyr)
library(tidyr)

recon.version='v3.1'
cost='r2'
style='real-time'
covrange='idp1'#'300km-20150313Fri'
scalesnotel='scale'
fscaMatch='wofsca'
dateflag='B'
residblending='unblended'

graphbase=paste0('graphs/rswe_',recon.version,'/covrange',covrange,'/snotel',scalesnotel,'/',dateflag,'/',fscaMatch,'/',style,'/')
swediffmap=point_plot_setup(recon.version,style,cost,covrange,scalesnotel,fscaMatch,dateflag)

swediffmap2=swediffmap %>%
     select(Station_ID,date,mnth,yrdoy,phv.pred,phvrcn.pred,swediff.phv,swediff.phvrcn,fsca,yr,snotel,swe,elev) %>%     
     mutate(
          phv.phvrcn.diff=phv.pred-phvrcn.pred,
          phvdiff.phvrcndiff.diff=swediff.phv-swediff.phvrcn,
          doy=as.numeric(strftime(date,'%j')),
          elevcut=cut(elev,breaks=seq(1500,3500,500),dig.lab=5))

indh=which(swediffmap2$phv.phvrcn.diff > 0.35)
unique(swediffmap2$Station_ID[indh])
ind=which(swediffmap2$phv.phvrcn.diff < -0.75)
unique(swediffmap2$Station_ID[ind])    


## swe-phvdiff_phvrcndiff_diff-bin2d ---- 
ggplot(swediffmap2)+
     geom_bin2d(aes(x=swe,y=phvdiff.phvrcndiff.diff,fill=..density..),binwidth=c(.1,.1))+
     ggtitle('phvdiff - phvrcndiff, density plot')
ggsave(paste0(graphbase,'swe-phvdiff_phvrcndiff_diff-bin2d.png'))


ggplot(swediffmap2)+geom_point(aes(x=elev,y=phv.phvrcn.diff,colour=swe))+
     ggtitle('SWE estimate difference PHV-PHVRCN')
ggsave(paste0(graphbase,'phv_phvrcn_diff-elev_colourswe.png'))

model_names <- c(
     "swediff.phv" = "PHV-baseline",
     "swediff.phvrcn" = "PHV-RCN"
)

## elev-swediff_facetmodel_colourswe.png -------
mypal=c(brewer.pal(7,'PuBu')[2:7])
avgbias=swediffmap2 %>% 
     gather(model,swediff,swediff.phv:swediff.phvrcn) %>% 
     group_by(Station_ID,elev,model) %>%
     summarise(
          avgbias=mean(swediff,na.rm=T)
     )

lineeqns=avgbias %>%
     group_by(model) %>%
     do(
          data_frame(eqn=lm_eqn(glm(formula=avgbias~elev,data=.,family=gaussian(link='identity'))))
                 )
swediffmap2 %>%
     gather(model,swediff,swediff.phv:swediff.phvrcn) %>% 
     mutate(swecut=cut(swe,breaks=c(0,.1,.25,.5,1,2,4))) %>% 
     {
          ggplot()+
               geom_point(data=.,aes(x=elev,y=swediff,colour=swecut),alpha=0.9)+
               geom_point(data=avgbias,aes(x=elev,y=avgbias),colour='red')+
               geom_smooth(data=avgbias,aes(x=elev,y=avgbias),colour='red',method='lm',formula='y~x',se=FALSE)+
               geom_text(data=lineeqns,aes(x=2000,y=0.9,label=eqn),parse=TRUE,hjust=0)+
               facet_wrap(~model,labeller=as_labeller(model_names))+
               scale_colour_manual(values=mypal)+
               labs(y='Bias [m]',x='Elevation [m]')+
               guides(colour=guide_legend('SWE Depth [m]',override.aes=list(size=3)))+
               theme_bw(base_size=18)+
               theme(legend.key=element_rect(colour=NA),
                    strip.background=element_blank())
               # ggtitle('model - snotel*fsca')
     }
ggsave(paste0(graphbase,'elev-swediff_facetmodel_colourswe.png'),dpi=600,width=12,height=7)
    



ggplot(swediffmap2)+
     geom_bin2d(aes(x=swe,y=phv.pred,fill=..density..),binwidth=0.1)+
     geom_abline()
ggsave(paste0(graphbase,'phv_snotel_bin2d.png'))

ggplot(swediffmap2)+
     geom_bin2d(aes(x=swe,y=phvrcn.pred,fill=..density..),binwidth=0.1)+
     geom_abline()
ggsave(paste0(graphbase,'phvrcn_snotel_bin2d.png'))
