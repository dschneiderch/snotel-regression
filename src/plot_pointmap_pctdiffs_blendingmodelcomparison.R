library('ProjectTemplate')
setwd('~/GoogleDrive/snotel-regression_project')
load.project()
library(dplyr)

recon.version='v3.2'
cost='r2'
style='real-time'
covrange='300km'
swediffmap=point_plot_setup(recon.version,style,cost,covrange)

#### Combined mapping of pct diffs- rows = modeltype, columns 4 years.
facet1_names=list(
    'swediffpct.phv'='Regression w/o RCN',
    'swediffpct.phvrcn'='Regression w/ RCN',
    'regressiondiff'='Difference')
facet2_names=list(
  'unblended'='Unblended',
  'blended'='Blended',
  'blenddiff'='Difference')
# plot_labeller <- function(variable,value){
#     if(variable=='yr'){
#         return(value)
#     } else {
#         return(facet1_names[value])
#     }
#  }
plot_labeller <- function(variable,value){
    if(variable=='blended'){
        return(facet2_names[value])
    } else if (variable=='variable'){
        return(facet1_names[value])
    } else {
      return(as.character(value))
    }
 }

yer=2001
for(yer in 2001:2012){
  print(yer)
monthday= 'Apr01'
brkpoints=9#odd number
outlier=20#interval = outlier/floor(brkpoints/2)
swediffmap$phvblendndiff=(abs(swediffmap$swediffpct.phv)-abs(swediffmap$swediffpct.phvfull))
swediffmap$phvrcnblendndiff=(abs(swediffmap$swediffpct.phvrcn)-abs(swediffmap$swediffpct.phvrcnfull))
swediffmap$noblendrcndiff=(abs(swediffmap$swediffpct.phv)-abs(swediffmap$swediffpct.phvrcn))
swediffmap$blendrcndiff=(abs(swediffmap$swediffpct.phvfull)-abs(swediffmap$swediffpct.phvrcnfull))

swediffmapplot=subset(swediffmap,mnthdy==monthday & yr==yer)# | yr==2011))
swediffmapplot=melt(swediffmapplot[,c('yr','lat','long','phvblendndiff','phvrcnblendndiff','noblendrcndiff','blendrcndiff','swediffpct.phv','swediffpct.phvrcn','swediffpct.phvfull','swediffpct.phvrcnfull')],id=c('yr','lat','long'))

swediffmapplot$blended='unblended'
fullind=c(grep('full',as.character(swediffmapplot$variable)),grep('blendrcndiff',as.character(swediffmapplot$variable)))
diffind=grep('ndiff',as.character(swediffmapplot$variable))
swediffmapplot$blended[fullind]='blended'
swediffmapplot$blended[diffind]='diff'
swediffmapplot$blended=factor(swediffmapplot$blended,levels=c('unblended','blended','diff'))

swediffmapplot$variable=factor(swediffmapplot$variable,levels=c('swediffpct.phv','swediffpct.phvrcn','regressiondiff'))
swediffmapplot$variable[swediffmapplot$variable=='swediffpct.phvfull']='swediffpct.phv'
swediffmapplot$variable[swediffmapplot$variable=='swediffpct.phvrcnfull']='swediffpct.phvrcn'
swediffmapplot$variable[swediffmapplot$variable=='phvblendndiff']='swediffpct.phv'
swediffmapplot$variable[swediffmapplot$variable=='phvrcnblendndiff']='swediffpct.phvrcn'
swediffmapplot$variable[swediffmapplot$variable=='noblendrcndiff']='regressiondiff'
swediffmapplot$variable[swediffmapplot$variable=='blendrcndiff']='regressiondiff'
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
myplot=function(yr){
	ggplot()+
    geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
    scale_alpha_continuous(range=c(1,0.5))+
    geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
    geom_point(data=subset(swediffmapplot,value < 0),aes(x=long,y=lat,color=cuts,shape=cuts),alpha=0.75,size=2)+
    geom_point(data=subset(swediffmapplot,value >= 0),aes(x=long,y=lat,color=cuts,shape=cuts),alpha=0.75,size=2)+
    scale_color_manual(values = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(swediffmapplot$cuts)))),drop=F,labels=legendlabels)+
    scale_shape_manual(values=rev(c(rep(16,length(levels(swediffmapplot$cuts))/2),rep(17,length(levels(swediffmapplot$cuts))/2))),drop=F,labels=legendlabels)+
    scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F,labels=legendlabels)+
    labs(x='Longitude',y='Latitude')+
    guides(alpha=FALSE,
           size=guide_legend('SWE\nDifferences (%)',reverse=F),
           colour=guide_legend('SWE\nDifferences (%)',reverse=F),
           shape=guide_legend('SWE\nDifferences (%)'),reverse=F)+
    coord_fixed(ratio=1,xlim=c(-112.25,-104.125),ylim=c(33,43.75))+
    theme_minimal()+
    theme(legend.key=element_rect(fill='grey40',colour='grey40'),
          legend.key.size=unit(1.5,'lines'),
          legend.background=element_rect(fill='grey40',colour='grey40'),
          legend.text=element_text(size=14,colour='grey95'),
          legend.title=element_text(size=16,colour='grey95'),
          axis.text=element_text(size=14),
          plot.title=element_text(size=18),
          axis.title=element_text(size=16),
          strip.text=element_text(size=14,face='bold'))+
    ggtitle(paste0(yr,' April 1 % Differences with Observed SNOTEL SWE\nBlue triangles in Difference map indicate a smaller error with Regression w/ RCN\n'))+
    facet_grid(variable~blended,labeller=plot_labeller)
}
mp=myplot(yer)
show(mp)
ggsave(plot=mp,filename=paste0('graphs/rswe_',recon.version,'/blendingmodelcomparison.5percent.pctdiff.',yer,'.',covrange,'.',cost,'.png'),dpi=300,width=11,height=10)
}

