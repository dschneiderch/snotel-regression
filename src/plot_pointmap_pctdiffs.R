library('ProjectTemplate')
setwd('~/GoogleDrive/snotel-regression_project')
load.project()
library(dplyr)

recon.version='v3.2'
cost='rmse'
style='real-time'
swediffmap=point_plot_setup(recon.version,style,cost)

bias.m=melt(swediffmap[,c('yr','swediff.phv','swediff.phvrcn','swediff.phvfull','swediff.phvrcnfull')],.(yr))
dcast(bias.m,variable~yr,value.var='value',mean,na.rm=T)

##-------- PERCENT DIFFERENCES IN SWE, facet by year.
pctdiff=function(dF,colnm,monthday){
	dFplot=subset(dF,mnthdy==monthday)
	dFplot[,colnm]=dFplot[,colnm]*100
	dFplot$cuts=cut(dFplot[,colnm],breaks=seq(-200,200,length.out=9))
#
	dFplot$extreme=NA
	dFplot$extreme[dFplot[,colnm] > 200] ='> 200%'
	dFplot$extreme[dFplot[,colnm] < -200 | is.na(dFplot[,colnm])] = '< -200%'
#
	extremesz=2
	szmin=2
	szstep=1
	szmax=szmin+szstep*(length(levels(dFplot$cuts))-2)/2
#
	g=ggplot()+
	geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
	scale_alpha_continuous(range=c(1,0.5))+
	geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
	geom_point(data=subset(dFplot,colnm > 200 | colnm < -200 | is.na(colnm)),aes(x=long,y=lat,shape=extreme),size=3,fill='black',alpha=0.5)+
	geom_point(data=subset(dFplot,colnm <= 200 | colnm >= -200),aes(x=long,y=lat,color=cuts,size=cuts),alpha=0.75)+
	scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F)+
	scale_color_manual(values = colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(dFplot$cuts))),drop=F)+
	scale_shape_manual(values=c(25,24),drop=F)+
	labs(x='Longitude',y='Latitude')+
	guides(shape=guide_legend('Extremes',reverse=T),
		alpha=FALSE,
		size=guide_legend('SWE\nDifferences (%)',reverse=T),
		colour=guide_legend('SWE\nDifferences (%)',reverse=T))+
	    coord_fixed(ratio=1,xlim=c(-112.25,-104.125),ylim=c(33,43.75))+	
	facet_wrap(~yr)+
	theme_bw()+
	  theme(legend.key=element_rect(fill='grey20',colour='grey20'),
          legend.key.size=unit(1.5,'lines'),
          legend.background=element_rect(fill='grey20'),
          legend.text=element_text(size=14,colour='grey90'),
          legend.title=element_text(size=16,colour='grey90'),
          axis.text=element_text(size=14),
          plot.title=element_text(size=18),
          axis.title=element_text(size=16),
          strip.text=element_text(size=14,face='bold'))
	# theme(legend.key=element_rect(fill='grey20'),
	# 	legend.background=element_rect(fill='grey20')
	 
	return(g)
}

gphv=pctdiff(swediffmap,'swediffpct.phv','Apr01')
gphv+ggtitle('April 1 % Differences with Observed SNOTEL SWE\nRegression with only Physiographics (unblended)')

ggsave(file=paste0('graphs/rswe_',recon.version,'/swediffpct.phv.map.pdf'),height=12,width=12)

gphvrcn=pctdiff(swediffmap,'swediffpct.phvrcn','Apr01')
gphvrcn+ggtitle('April 1 % Differences with Observed SNOTEL SWE\nRegression with Physiographics and Reconstruction (unblended)')

ggsave(file=paste0('graphs/rswe_',recon.version,'/swediffpct.phvrcn.map.pdf'),height=12,width=12)

grecon=pctdiff(swediffmap,'swediffpct.recon','Apr01')
grecon+ggtitle('April 1 % Differences with Observed SNOTEL SWE\nReconstruction (unblended)')

ggsave(file=paste0('graphs/rswe_',recon.version,'/swediffpct.recon.map.pdf'),height=12,width=12)

gphv=pctdiff(swediffmap,'swediffpct.phvfull','Apr01')
gphv+ggtitle('April 1 % Differences with Observed SNOTEL SWE\nRegression with only Physiographics (Blended)')

ggsave(file=paste0('graphs/rswe_',recon.version,'/swediffpct.phvfull.map.pdf'),height=12,width=12)

gphvrcn=pctdiff(swediffmap,'swediffpct.phvrcnfull','Apr01')
gphvrcn+ggtitle('April 1 % Differences with Observed SNOTEL SWE\nRegression with Physiographics and Reconstruction (Blended)')

ggsave(file=paste0('graphs/rswe_',recon.version,'/swediffpct.phvrcnfull.map.pdf'),height=12,width=12)

grecon=pctdiff(swediffmap,'swediffpct.reconfull','Apr01')
grecon+ggtitle('April 1 % Differences with Observed SNOTEL SWE\nReconstruction (Blended)')

ggsave(file=paste0('graphs/rswe_',recon.version,'/swediffpct.reconfull.map.pdf'),height=12,width=12)