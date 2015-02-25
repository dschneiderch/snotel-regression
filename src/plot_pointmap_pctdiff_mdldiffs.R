library('ProjectTemplate')
setwd('~/GoogleDrive/snotel-regression_project')
load.project()
library(dplyr)

recon.version='v3.2'
swediffmap=point_plot_setup(recon.version)

## Day of Year for plots
monthday='Apr01'
swediffmapplot=subset(swediffmap,mnthdy==monthday)#& (yr==2007 | yr==2011))
#

## Some parameters for plotting the data
brkpoints=11
outlier=250
## Map differences between %differences of the regression models facet by year
swediffmapplot$mdldiffpct.phv.phvrcn=(abs(swediffmapplot$swediffpct.phv)-abs(swediffmapplot$swediffpct.phvrcn))*100
swediffmapplot$mdldiffpct.phv.phvrcn[swediffmapplot$mdldiffpct.phv.phvrcn > outlier]=outlier-1
swediffmapplot$mdldiffpct.phv.phvrcn[swediffmapplot$mdldiffpct.phv.phvrcn < -outlier]=-outlier+1
swediffmapplot$cuts=cut(swediffmapplot$mdldiffpct.phv.phvrcn,breaks=seq(-outlier,outlier,length.out=brkpoints),right=F, include.lowest=T)
swediffmapplot$cuts=factor(swediffmapplot$cuts,levels=rev(levels(swediffmapplot$cuts)))
#
legendlabels=levels(swediffmapplot$cuts)
outlabel=outlier*2/(brkpoints-1)*(brkpoints-1)/2-outlier*2/(brkpoints-1)
legendlabels[1]=paste('>= ',outlabel,sep='')
legendlabels[length(levels(swediffmapplot$cuts))]=paste('< ',-outlabel,sep='')
#
extremesz=4
szmin=2
szstep=1
szmax=ceiling((szmin+szstep*length(levels(swediffmapplot$cuts)))/2) 

myplot=function(){
ggplot()+
    geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
    scale_alpha_continuous(range=c(1,0.5))+
    geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
    geom_point(data=subset(swediffmapplot,mdldiffpct.phv.phvrcn >= 0),aes(x=long,y=lat,size=cuts,color=cuts,shape=cuts),alpha=0.75)+
    geom_point(data=subset(swediffmapplot,mdldiffpct.phv.phvrcn < 0),aes(x=long,y=lat,size=cuts,color=cuts,shape=cuts),alpha=0.75)+
    scale_color_manual(values = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(swediffmapplot$cuts)))),drop=F, labels=legendlabels)+
    scale_shape_manual(values=rev(c(rep(16,length(levels(swediffmapplot$cuts))/2),rep(17,length(levels(swediffmapplot$cuts))/2))),drop=F,labels=legendlabels)+
    scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F, labels=legendlabels)+
    labs(x='Longitude', y='Latitude')+
    guides(alpha=F,
           size=guide_legend('Error\nDifferences (%)'),
           colour=guide_legend('Error\nDifferences (%)',reverse=F),
           shape=guide_legend('Error\nDifferences (%)'),reverse=F)+
    coord_fixed(ratio=1,xlim=c(-112.25,-104.125),ylim=c(33,43.75))+
    theme_bw()+
    theme(legend.key=element_rect(fill='grey80'),
          legend.key.size=unit(1.5,'lines'),
          legend.background=element_rect(fill='grey80'),
          legend.text=element_text(size=14),
          legend.title=element_text(size=16),
          axis.text=element_text(size=14),
          plot.title=element_text(size=18),
          axis.title=element_text(size=16),
          strip.text=element_text(size=14,face='bold'))+
    facet_wrap(~yr)+
 ggtitle(paste('April 1 Differences in % Error between \nRegression w/o RCN and Regression w/ RCN\nRelative to Observed SNOTEL SWE\nBlue triangles indicate a smaller error with Regression w/ RCN',sep=''))
}
myplot()   #

ggsave(file=paste0('graphs/rswe_',recon.version,'/mdldiff.phv.phvrcn.allyears.png'),dpi=300,width=12, height=12)
ggsave(file=paste0('graphs/rswe_',recon.version,'/mdldiff.phv.phvrcn.2002.2011.png'),dpi=300,width=14, height=7)