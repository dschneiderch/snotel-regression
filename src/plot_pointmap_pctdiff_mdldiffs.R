library('ProjectTemplate')
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
graphbase=paste0('graphs/rswe_',recon.version,'/covrange',covrange,'/snotel',scalesnotel,'/',dateflag,'/',fscaMatch,'/',style,'/')

## Day of Year for plots
monthday='Apr01'
swediffmapplot=subset(swediffmap,mnthdy==monthday)#& (yr==2007 | yr==2011))
#
## Some parameters for plotting the data
brkpoints=9#odd number
outlier=50#interval = outlier/floor(brkpoints/2)
## Map differences between %differences of the regression models facet by year
if(residblending=='unblended'){
    swediffmapplot$mdldiffpct.phv.phvrcn=(abs(swediffmapplot$swediffpct.phv)-abs(swediffmapplot$swediffpct.phvrcn))*100  
} else { 
swediffmapplot$mdldiffpct.phv.phvrcn=(abs(swediffmapplot$swediffpct.phvfull)-abs(swediffmapplot$swediffpct.phvrcnfull))*100  
}

swediffmapplot$mdldiffpct.phv.phvrcn[swediffmapplot$mdldiffpct.phv.phvrcn > outlier]=outlier-1
swediffmapplot$mdldiffpct.phv.phvrcn[swediffmapplot$mdldiffpct.phv.phvrcn < -outlier]=-outlier+1
swediffmapplot$cuts=cut(swediffmapplot$mdldiffpct.phv.phvrcn,breaks=seq(-outlier,outlier,length.out=brkpoints),right=F, include.lowest=T)
# swediffmapplot$cuts=factor(swediffmapplot$cuts,levels=levels(swediffmapplot$cuts))
#

### ---- plot model differencesin percent biases for each year (facted) on a map
##
##
#
# legendlabels=levels(swediffmapplot$cuts)
# outlabel=outlier*2/(brkpoints-1)*(brkpoints-1)/2-outlier*2/(brkpoints-1)
# legendlabels[1]=paste('< ',-outlabel,sep='')
# legendlabels[length(levels(swediffmapplot$cuts))]=paste('>= ',outlabel,sep='')
# #
# extremesz=4
# szmin=2
# szstep=1
# szmax=ceiling((szmin+szstep*length(levels(swediffmapplot$cuts)))/2) 

# # myplot=function(){
# ggplot()+
#     geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
#     scale_alpha_continuous(range=c(1,0.5))+
#     geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
#     geom_point(data=subset(swediffmapplot,mdldiffpct.phv.phvrcn < 0),aes(x=long,y=lat,color=cuts,shape=cuts),alpha=0.75,size=3.5)+
#     geom_point(data=subset(swediffmapplot,mdldiffpct.phv.phvrcn >= 0),aes(x=long,y=lat,color=cuts,shape=cuts),alpha=0.75,size=3.5)+
#     scale_color_manual(values = colorRampPalette(brewer.pal(11,"PuOr"))(length(levels(swediffmapplot$cuts))),drop=F, labels=legendlabels)+
#     scale_shape_manual(values=c(rep(16,length(levels(swediffmapplot$cuts))/2),rep(17,length(levels(swediffmapplot$cuts))/2)),drop=F,labels=legendlabels)+
#     # scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F, labels=legendlabels)+
#     labs(x='Longitude', y='Latitude')+
#     guides(alpha=F,
#            # size=guide_legend('Error\nDifferences (%)'),
#            colour=guide_legend('Relative %Bias\nDifferences',reverse=T),
#            shape=guide_legend('Relative %Bias\nDifferences',reverse=T),override.aes=list(alpha=1,size=5))+
#     coord_fixed(ratio=1,xlim=c(-112.25,-104.125),ylim=c(33,43.75))+
#     theme_minimal()+
#     theme(legend.key=element_rect(fill='grey40',colour='grey40'),
#           legend.key.size=unit(1.5,'lines'),
#           legend.background=element_rect(fill='grey40',colour='grey40'),
#           legend.text=element_text(size=14,colour='grey90'),
#           legend.title=element_text(size=16,colour='grey90'),
#           axis.text=element_text(size=14),
#           plot.title=element_text(size=18),
#           axis.title=element_text(size=16),
#           strip.text=element_text(size=14,face='bold'))+
#     facet_wrap(~yr)+
#  ggtitle(paste('April 1 Differences in Relative %Bias between \nRegression w/o RCN and Regression w/ RCN\nRelative to Observed SNOTEL SWE\nPurple triangles indicate a smaller error with Regression w/ RCN\n',sep=''))
# }
# myplot()   #

# ggsave(file=paste0(graphbase,'/mdldiff.phv.phvrcn.allyears_',cost,'.png'),dpi=300,width=12, height=12)
# #ggsave(file=paste0('graphs/rswe_',recon.version,'/mdldiff.phv.phvrcn.2002.2011.png'),dpi=300,width=14, height=7)

swediffmapplot$improvcount=ifelse(swediffmapplot$mdldiffpct.phv.phvrcn>0,1,0)
mdlimprov=ddply(swediffmapplot,.(Station_ID),summarise,
  count=sum(improvcount),
  lat=mean(lat),
  long=mean(long),
  elev=mean(elev))
mdlimprov=arrange(mdlimprov,count)
mdlimprov$countper=ifelse(mdlimprov$count<=2,'0%-25%',
                                            ifelse(mdlimprov$count>2 & mdlimprov$count<5,'25%-50%',
                                              ifelse(mdlimprov$count==5,'50%',
                                            ifelse(mdlimprov$count>5 & mdlimprov$count<=7,'50%-75%',
                                            ifelse(mdlimprov$count>7 & mdlimprov$count<=12,'75%-100%',NA)))))
# mdlimprov$count=factor(mdlimprov$count,levels=rev(seq(0,max(mdlimprov$count))))

mdlimprov$elevcut=cut(mdlimprov$elev,breaks=seq(1500,4000,500),dig.lab=4)
mdlimprov$per=mdlimprov$count/12
head(mdlimprov)
str(mdlimprov)
tmp2=merge(snotellocs,mdlimprov)
mapView(tmp2)

dF=as.data.frame(table(mdlimprov$count))
sum(dF$Freq[8:13])
sum(dF$Freq[8:13])/237
sum(dF$Freq[1:6])
sum(dF$Freq[1:6])/237
ggplot(mdlimprov)+geom_text(aes(x=long,y=lat,label=count))

table(mdlimprov$elevcut)
outer=subset(mdlimprov,elev<2000 | elev>3500)
outer

mid=subset(mdlimprov,elev>2500 & elev<3000)
mean(mid$per)

dcast(mdlimprov,elevcut~.,value.var='per',fun.agg=mean)

ggplot(mdlimprov)+geom_boxplot(aes(x=elevcut,y=count))+scale_x_discrete(drop=F)
ggplot(mdlimprov)+geom_bar(aes(x=elevcut))


## ind[LMU] subsets represent 1 std deviation cutoff of the snotel elevation distribution.
indL=which(mdlimprov$elev<2440)
mdlimprov[indL,]
length(indL)
mean(mdlimprov[indL,'per'])
indM=which(mdlimprov$elev>2440 & mdlimprov$elev<3160)
length(indM)
mean(mdlimprov[indM,'per'])
indU=which(mdlimprov$elev>3160)
length(indU)
mean(mdlimprov[indU,'per'])

avgper=dcast(mdlimprov,elevcut~.,mean,value.var='per')
names(avgper)=c('ElevationBand','AveragePerformance')
mean(avgper$AveragePerformance[1:3])
dev.new()
ggplot(avgper)+geom_bar(aes(ElevationBand,AveragePerformance),stat='identity')

ggplot()+
    geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
      scale_alpha_continuous(range=c(1,0.5))+
    geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
    geom_point(data=mdlimprov,aes(x=long,y=lat,color=as.factor(countper)),alpha=.95,size=6,shape=16)+
    scale_color_manual(values=(brewer.pal(length(unique(mdlimprov$countper)),'Purples')),drop=F)+
    # scale_color_manual(values = rev(c('#b15928','#ffff99','#ff7f00','#fdbf6f','#e31a1c','#fb9a99','#6a3d9a','#cab2d6','#1f78b4','#a6cee3','#33a02c','#b2df8a','black')),drop=F)+##use this for count variable
    labs(x='Longitude', y='Latitude')+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    coord_fixed(ratio=1,ylim=c(33,43.75),xlim=c(-112.25,-104.125))+
    # expand_limits(x=c(-112.25,-104.125))+
    guides(alpha=F,
           colour=guide_legend('% of Yrs',reverse=T,override.aes=list(alpha=1,size=5)))+
    theme_bw()+
    theme(legend.key=element_rect(fill='grey40',colour='grey40'),
          legend.key.size=unit(1.5,'lines'),
          legend.background=element_rect(fill='grey40',colour='grey40'),
          legend.text=element_text(size=14,colour='grey90'),
          legend.title=element_text(size=16,colour='grey90'),
          legend.justification=c(1,0),
          legend.position=c(1,0),
          axis.ticks=element_line(colour='black',size=1),
          axis.text=element_text(size=14),
          plot.title=element_text(size=18),
          panel.grid=element_blank(),
          panel.background=element_blank(),
          plot.background=element_rect(colour=NA),
          axis.title=element_text(size=16))#+
      # ggtitle(paste('% of Years at each SNOTEL station where\nPHV-RCN decreases % Bias\ncompared to PHV-baseline\n',sep=''))

ggsave(file=paste0(graphbase,'percentyrs_improvement_r2_realtime_',monthday,'_',residblending,'.png'),dpi=600,width=12, height=12)

### ----- Check spatial autocorrelation of regional performance
## are the counts spatially correlated?
mdlimprovsp=mdlimprov
coordinates(mdlimprovsp)=~long+lat
proj4string(mdlimprovsp)='+proj=longlat +datum=WGS84'
mdlimprovsp=spTransform(mdlimprovsp,CRS('+init=epsg:5070'))

     kd=dnearneigh(mdlimprovsp,d1=0,d2=300000,row.names=mdlimprovsp$Station_ID)#d2 should equal abut 200km fro all stations to have some neighbors. min 8 with 200km.
     #summary(kd)
    # plot(kd,coords=coordinates(snotellocs.usgs))#
     getneighbors=function(nm){
          ind=grep(nm,attr(kd,'region.id'))
          vec=kd[[ind]]
          data.frame(station=nm,neighbors=attr(kd,'region.id')[vec])
     } 
     snotelneighbors=ldply(as.list(attr(kd,'region.id')),getneighbors)
#
     dist=nbdists(kd,mdlimprovsp)
     p=.5#fassnacht found 0.5 to be optimal until mid-march, then up to 1.3 at end-april then back down to .5.
     idw=lapply(dist,function(x) 1/(x^p))
     kd.idw=nb2listw(kd,glist=idw,style='B')
     moran.test(mdlimprovsp$count,listw=kd.idw)

### ------


### ---- Plot histogram of pct bias differences
mean(swediffmapplot$mdldiffpct.phv.phvrcn)
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
brks=seq(datamin-1,datamax+1,breakint)
brklabs=brks
brklabs[1]=paste('<= ',-outlier,sep='')
brklabs[length(brklabs)]=paste('>= ',outlier,sep='')
ggplot()+
  geom_histogram(data=swediffmapplot,aes(x=mdldiffpct.phv.phvrcn),breaks=brks,right=F,fill='grey80',colour='grey50')+
  geom_linerange(data=swediffstats,aes(x=avg,ymin=0,ymax=65),show.legend=TRUE,colour='red',size=1)+
  geom_linerange(data=swediffstats,aes(x=med,ymin=0,ymax=65),show.legend=TRUE,colour='blue',size=1)+
  geom_text(data=swediffstats,aes(x=-50,y=60,label=paste('avg=',round(avg,2))),colour='red',hjust=0,vjust=0)+
  geom_text(data=swediffstats,aes(x=-50,y=50,label=paste('med=',round(med,2))),colour='blue',hjust=0,vjust=.5)+
  geom_hline(yintercept=0,size=.8,color='black')+
  geom_text(data=data.frame(yr=seq(2001,2012)),aes(x=0,y=75,label=yr),fontface='bold',vjust=0)+
  labs(x='% Bias Difference')+
  scale_x_continuous(breaks=brks[seq(1,length(brks),2)],labels=brklabs[seq(1,length(brks),2)])+
  scale_y_continuous(expand=c(0,0),limits=c(0,90))+
  theme_minimal(base_size = 15)+
  theme(strip.text=element_blank(),#element_text(size=14,face='bold'),
        panel.grid=element_blank(),
        axis.line=element_line(size=.5))+
  facet_wrap(~yr,nrow=4)#,space='free_y',scales='free_y')

ggsave(filename=paste0(graphbase,'hist_',residblending,'mdldiffs_yearfacets_',cost,'.pdf'),width=13,height=9)

##