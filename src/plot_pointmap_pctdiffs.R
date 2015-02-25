library('ProjectTemplate')
setwd('~/GoogleDrive/snotel-regression_project')
load.project()
library(dplyr)

recon.version='v3.2'
swediffmap=point_plot_setup(recon.version)

##-------- PERCENT DIFFERENCES IN SWE, facet by year.
pctdiff=function(dF,colnm,monthday){
	dFplot=filter(dF,mnthdy==monthday)
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
	szmax=szmin+szstep*(length(levels(dFplot$cuts))-1)/2
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
		size=guide_legend('SWE Differences (%)',reverse=T),
		colour=guide_legend('SWE Differences (%)',reverse=T))+
	    coord_fixed(ratio=10/10,xlim=c(-112.25,-104.125),ylim=c(33,43.75))+	
	facet_wrap(~yr)+
	theme_bw()+
	theme(legend.key=element_rect(fill='grey80'),
		legend.background=element_rect(fill='grey80'))
	return(g)
}


gphv=pctdiff(swediffmap,'swediffpct.phv','Apr01')
gphv+ggtitle('April 1 % Differences with Observed SNOTEL SWE\nRegression with only Physiographics')

ggsave(file=paste0('graphs/rswe_',recon.version,'/swediffpct.phv.map.pdf'),height=12,width=12)

gphvrcn=pctdiff(swediffmap,'swediffpct.phvrcn','Apr01')
gphv+ggtitle('April 1 % Differences with Observed SNOTEL SWE\nRegression with Physiographics and Reconstruction')

ggsave(file=paste0('graphs/rswe_',recon.version,'/swediffpct.phvrcn.map.pdf'),height=12,width=12)

grecon=pctdiff(swediffmap,'swediffpct.recon','Apr01')
grecon+ggtitle('April 1 % Differences with Observed SNOTEL SWE\nReconstruction')

ggsave(file=paste0('graphs/rswe_',recon.version,'/swediffpct.recon.map.pdf'),height=12,width=12)