library('ProjectTemplate')
setwd('~/GoogleDrive/snotel-regression_project')
load.project()
library(dplyr)

recon.version='v3.2'
cost='r2'
style='real-time'
swediffmap=point_plot_setup(recon.version,style,cost)

#### Combined mapping of pct diffs- rows = modeltype, columns 4 years.
facet1_names=list(
    'swephvpct.phvfull'='Regression w/o RCN',
    'swephvpct.phvrcnfull'='Regression w/ RCN')
plot_labeller <- function(variable,value){
    if(variable=='yr'){
        return(value)
    } else {
        return(facet1_names[value])
    }
 }

monthday='Apr01'
brkpoints=9#odd number
outlier=100#interval = outlier/floor(brkpoints/2)#outlier=interval*floor(brkpoints/2)
swediffmapplot=subset(swediffmap,mnthdy==monthday)# & yr==2002)# | yr==2011))
swediffmapplot=melt(swediffmapplot[,c('yr','lat','long', 'swediffpct.phvfull','swediffpct.phvrcnfull')],id=c('yr','lat','long'))
#
swediffmapplot$value=swediffmapplot$value*100
swediffmapplot$value[swediffmapplot$value > outlier]=outlier-1
swediffmapplot$value[swediffmapplot$value < -outlier | is.na(swediffmapplot$value)]=-outlier+1
swediffmapplot$cuts=cut(swediffmapplot$value,breaks=seq(-outlier,outlier,length.out=brkpoints),right=F,include.lowest=T)
swediffmapplot$cuts=factor(swediffmapplot$cuts,levels=rev(levels(swediffmapplot$cuts)))
#
legendlabels=levels(swediffmapplot$cuts)
outlabel=outlier*2/(brkpoints-1)*((brkpoints-1)/2-1)
legendlabels[1]=paste('>= ',outlabel,sep='')
legendlabels[length(levels(swediffmapplot$cuts))]=paste('< ',-outlabel,sep='')
#
szmin=2
szstep=1
szmax=ceiling((szmin+szstep*length(levels(swediffmapplot$cuts)))/2)
#
myplot=function(){
	ggplot()+
    geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
    scale_alpha_continuous(range=c(1,0.5))+
    geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
    geom_point(data=subset(swediffmapplot,value < 0),aes(x=long,y=lat,color=cuts,shape=cuts),alpha=0.65,size=3)+
    geom_point(data=subset(swediffmapplot,value >= 0),aes(x=long,y=lat,color=cuts,shape=cuts),alpha=0.65,size=3)+
    scale_color_manual(values = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(swediffmapplot$cuts)))),drop=F,labels=legendlabels)+
    scale_shape_manual(values=c(rep(16,length(levels(swediffmapplot$cuts))/2),rep(17,length(levels(swediffmapplot$cuts))/2)),drop=F,labels=legendlabels)+
    # scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F,labels=legendlabels)+
    labs(x='Longitude',y='Latitude')+
    guides(alpha=FALSE,
           size=guide_legend('SWE\nDifferences (%)',reverse=F),
           colour=guide_legend('SWE\nDifferences (%)',reverse=F),
           shape=guide_legend('SWE\nDifferences (%)'),reverse=F,override.aes=list(alpha=1))+
    coord_fixed(ratio=1,xlim=c(-112.25,-104.125),ylim=c(33,43.75))+
    theme_minimal()+
    theme(legend.key=element_rect(fill='grey40',colour='grey40'),
          legend.key.size=unit(1.5,'lines'),
          legend.background=element_rect(fill='grey40',colour='grey40'),
          legend.text=element_text(size=14,colour='grey90'),
          legend.title=element_text(size=16,colour='grey90'),
          axis.text=element_text(size=14),
          plot.title=element_text(size=18),
          axis.title=element_text(size=16),
          strip.text=element_text(size=14,face='bold'))+
    ggtitle('April 1 % Differences\nwith Observed SNOTEL SWE\n')+
    facet_grid(variable~yr,labeller=plot_labeller)
}
myplot()

# ggsave(paste0('graphs/rswe_',recon.version,'/modelcomparison.pctdiff.2002.2011.png'),dpi=300,width=11,height=10)
ggsave(paste0('graphs/rswe_',recon.version,'/modelcomparison.pctdiff.allyears.png'),dpi=300,width=24,height=9.5)
