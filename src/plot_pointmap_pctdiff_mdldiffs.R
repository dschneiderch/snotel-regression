library('ProjectTemplate')
setwd('~/GoogleDrive/snotel-regression_project')
reload.project()
library(dplyr)

recon.version='v3.2'
cost='rmse'
style='real-time'
covrange='300km-20150313Fri'
swediffmap=point_plot_setup(recon.version,style,cost,covrange)


## Day of Year for plots
monthday='Apr15'
swediffmapplot=subset(swediffmap,mnthdy==monthday)#& (yr==2007 | yr==2011))
#
## Some parameters for plotting the data
brkpoints=9#odd number
outlier=20#interval = outlier/floor(brkpoints/2)
## Map differences between %differences of the regression models facet by year
swediffmapplot$mdldiffpct.phv.phvrcn=(abs(swediffmapplot$swediffpct.phvfull)-abs(swediffmapplot$swediffpct.phvrcnfull))*100
swediffmapplot$mdldiffpct.phv.phvrcn[swediffmapplot$mdldiffpct.phv.phvrcn > outlier]=outlier-1
swediffmapplot$mdldiffpct.phv.phvrcn[swediffmapplot$mdldiffpct.phv.phvrcn < -outlier]=-outlier+1
swediffmapplot$cuts=cut(swediffmapplot$mdldiffpct.phv.phvrcn,breaks=seq(-outlier,outlier,length.out=brkpoints),right=F, include.lowest=T)
# swediffmapplot$cuts=factor(swediffmapplot$cuts,levels=levels(swediffmapplot$cuts))
#
legendlabels=levels(swediffmapplot$cuts)
outlabel=outlier*2/(brkpoints-1)*(brkpoints-1)/2-outlier*2/(brkpoints-1)
legendlabels[1]=paste('< ',-outlabel,sep='')
legendlabels[length(levels(swediffmapplot$cuts))]=paste('>= ',outlabel,sep='')
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
    geom_point(data=subset(swediffmapplot,mdldiffpct.phv.phvrcn < 0),aes(x=long,y=lat,color=cuts,shape=cuts),alpha=0.75,size=3.5)+
    geom_point(data=subset(swediffmapplot,mdldiffpct.phv.phvrcn >= 0),aes(x=long,y=lat,color=cuts,shape=cuts),alpha=0.75,size=3.5)+
    scale_color_manual(values = colorRampPalette(brewer.pal(11,"PuOr"))(length(levels(swediffmapplot$cuts))),drop=F, labels=legendlabels)+
    scale_shape_manual(values=c(rep(16,length(levels(swediffmapplot$cuts))/2),rep(17,length(levels(swediffmapplot$cuts))/2)),drop=F,labels=legendlabels)+
    # scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F, labels=legendlabels)+
    labs(x='Longitude', y='Latitude')+
    guides(alpha=F,
           # size=guide_legend('Error\nDifferences (%)'),
           colour=guide_legend('Relative %Bias\nDifferences',reverse=T),
           shape=guide_legend('Relative %Bias\nDifferences',reverse=T),override.aes=list(alpha=1,size=5))+
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
    facet_wrap(~yr)+
 ggtitle(paste('April 1 Differences in Relative %Bias between \nRegression w/o RCN and Regression w/ RCN\nRelative to Observed SNOTEL SWE\nPurple triangles indicate a smaller error with Regression w/ RCN\n',sep=''))
}
myplot()   #

ggsave(file=paste0('graphs/rswe_',recon.version,'/mdldiff.phv.phvrcn.allyears.png'),dpi=300,width=12, height=12)
#ggsave(file=paste0('graphs/rswe_',recon.version,'/mdldiff.phv.phvrcn.2002.2011.png'),dpi=300,width=14, height=7)

swediffmapplot$improvcount=ifelse(swediffmapplot$mdldiffpct.phv.phvrcn>0,1,0)
mdlimprov=ddply(swediffmapplot,.(Station_ID),summarise,
  count=sum(improvcount),
  lat=mean(lat),
  long=mean(long))
mdlimprov=arrange(mdlimprov,count)
mdlimprov$count=factor(mdlimprov$count,levels=rev(seq(0,max(mdlimprov$count))))

ggplot()+
    geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
    scale_alpha_continuous(range=c(1,0.5))+
    geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
    geom_point(data=mdlimprov,aes(x=long,y=lat,color=as.factor(count)),alpha=.95,size=6,shape=16)+
    scale_color_manual(values = rev(c('#b15928','#ffff99','#ff7f00','#fdbf6f','#e31a1c','#fb9a99','#6a3d9a','#cab2d6','#1f78b4','#a6cee3','#33a02c','#b2df8a','black')),drop=F)+
    labs(x='Longitude', y='Latitude')+
    scale_x_continuous(limits=c(-112.25,-104.125))+
    scale_y_continuous(limits=c(33,43.75))+
    coord_fixed(ratio=1)+
    guides(alpha=F,
           colour=guide_legend('# Yrs',reverse=F,override.aes=list(alpha=1,size=5)))+
    theme_minimal()+
    theme(legend.key=element_rect(fill='grey40',colour='grey40'),
          legend.key.size=unit(1.5,'lines'),
          legend.background=element_rect(fill='grey40',colour='grey40'),
          legend.text=element_text(size=14,colour='grey90'),
          legend.title=element_text(size=16,colour='grey90'),
          axis.text=element_text(size=14),
          plot.title=element_text(size=18),
          axis.title=element_text(size=16))+
 ggtitle(paste('# Years at each SNOTEL station where\nRegression w/ RCN decreases % Bias\ncompared to Regression w/o RCN\n',sep=''))

ggsave(file=paste0('graphs/rswe_',recon.version,'/covrange','300km','/yrs_improvement_r2_realtime_',monthday,'.png'),dpi=600,width=12, height=12)

swediffstats=ddply(swediffmapplot,.(yr),summarise,
  avg=mean(mdldiffpct.phv.phvrcn,na.rm=T),
  med=median(mdldiffpct.phv.phvrcn,na.rm=T))

#alldF=data.frame(yr="all",avg=mean(swediffmapplot$mdldiffpct.phv.phvrcn,na.rm=T),med=median(swediffmapplot$mdldiffpct.phv.phvrcn,na.rm=T))
#swediffstats=rbind(swediffstats,alldF)
#swediffstats$yr=factor(swediffstats$yr,levels=unique(swediffstats$yr))
#swediffstats
#str(swediffstats)
datamin=floor(min(swediffmapplot$mdldiffpct.phv.phvrcn,na.rm=T))
datamax=ceiling(max(swediffmapplot$mdldiffpct.phv.phvrcn,na.rm=T))
breakint=5
brks=seq(datamin,datamax,breakint)
ggplot()+
  geom_histogram(data=swediffmapplot,aes(x=mdldiffpct.phv.phvrcn),breaks=brks,right=F,fill='grey80',colour='grey50')+
  geom_linerange(data=swediffstats,aes(x=avg,ymin=0,ymax=65),show_guide=T,colour='red',size=1)+
  geom_linerange(data=swediffstats,aes(x=med,ymin=0,ymax=65),show_guide=T,colour='blue',size=1)+
  geom_text(data=swediffstats,aes(x=datamin,y=50,label=paste('avg=',round(avg,2))),colour='red',hjust=0,vjust=1)+
  geom_text(data=swediffstats,aes(x=datamin,y=40,label=paste('med=',round(med,2))),colour='blue',hjust=0,vjust=1)+
  scale_x_continuous(breaks=brks[seq(1,length(brks),2)])+
  theme_minimal()+
  theme(strip.text=element_text(size=14,face='bold'),
        panel.grid=element_blank())+
  facet_wrap(~yr,nrow=4)#,space='free_y',scales='free_y')

ggsave(filename=paste0('graphs/rswe_',recon.version,'/hist_mdldiffs_yearfacets.pdf'),width=13,height=9)