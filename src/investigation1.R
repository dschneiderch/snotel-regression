# setwd('/Volumes/PIA/Documents/snotel-regression_project')
setwd('~/Documents/snotel-regression_project')
library('ProjectTemplate')
load.project()

recon.version='v3.1'
yr='2012'
fid=file('~/Downloads/recon_v2015/fsca.dat','rb')
t=32
dFall=data.frame()
for(t in 1:184){
r=readgrads(t,fid)
dF=data.frame(var=extract(r,snotellocs),doy=t+(31+29),yr=2012,Station_ID=snotellocs$Station_ID,Site_Name=snotellocs$Site_Name)
dFall=rbind(dFall,dF)
}
close(fid)

dFall$date=strptime(paste0(dFall$yr,dFall$doy),'%Y%j',tz='MST')
sitenames=c('NIWOT','BERTHOUD SUMMIT','RED MOUNTAIN PASS','WOLF CREEK SUMMIT')
ind=which(dFall$Site_Name %in% sitenames)
fsca=dFall[ind,]
fg=ggplot(fsca)+
	geom_path(aes(x=date,y=var,colour=Site_Name),size=1,alpha=0.5)+
	scale_colour_brewer(palette='Set1')+
	ggtitle(paste0('fsca-',recon.version))+
	labs(y='fsca')+
	theme_bw()

fid=file('~/Downloads/recon_v2015/swe.dat','rb')
t=32
dFall=data.frame()
for(t in 1:184){
r=readgrads(t,fid)
dF=data.frame(var=extract(r,snotellocs),doy=t+(31+29),yr=2012,Station_ID=snotellocs$Station_ID,Site_Name=snotellocs$Site_Name)
dFall=rbind(dFall,dF)
}
dFall$date=strptime(paste0(dFall$yr,dFall$doy),'%Y%j',tz='MST')
sitenames=c('NIWOT','BERTHOUD SUMMIT','RED MOUNTAIN PASS','WOLF CREEK SUMMIT')
ind=which(dFall$Site_Name %in% sitenames)
swe=dFall[ind,]
sg=ggplot(swe)+
	geom_path(aes(x=date,y=var,colour=Site_Name),size=1,alpha=0.5)+
	scale_colour_brewer(palette='Set1')+
	ggtitle(paste0('swe-',recon.version))+
	labs(y='recon swe')+
	theme_bw()

g=arrangeGrob(fg,sg)
ggsave(plot=g,file=paste0('data/recon_',recon.version,'/snotelstation_timeseries_',recon.version,'_',yr,'v2015.pdf'),height=7,width=7)
 


recon.version='v3.2'
yr=2010
load(paste0('data/recon_',recon.version,'/snotelfscadata_',yr,'_',recon.version,'.RData'))
sitenames=c('NIWOT','BERTHOUD SUMMIT','RED MOUNTAIN PASS','WOLF CREEK SUMMIT')
ind=which(snotelrecon$Site_Name %in% sitenames)
sr=snotelrecon[ind,]
fg=ggplot(sr)+
	geom_path(aes(x=date,y=recon,colour=Site_Name),size=1,alpha=0.5)+
	scale_colour_brewer(palette='Set1')+
	ggtitle(paste0('fsca-',recon.version))+
	labs(y='fsca')+
	theme_bw()

load(paste0('data/recon_',recon.version,'/snotelrecondata_',yr,'_',recon.version,'.RData'))
sitenames=c('NIWOT','BERTHOUD SUMMIT','RED MOUNTAIN PASS','WOLF CREEK SUMMIT')
ind=which(snotelrecon$Site_Name %in% sitenames)
sr=snotelrecon[ind,]
rg=ggplot(sr)+
	geom_path(aes(x=date,y=recon,colour=Site_Name),size=1,alpha=0.5)+
	scale_colour_brewer(palette='Set1')+
	ggtitle(paste0('swe-',recon.version))+
	labs(y='recon swe')+
	theme_bw()
g=arrangeGrob(fg,rg)
ggsave(plot=g,file=paste0('data/recon_',recon.version,'/snotelstation_timeseries_',recon.version,'_',yr,'.pdf'),height=7,width=7)
 
yr=2012
s=raster::stack(paste0('data/recon_v3.1/reconfscadata_',yr,'_v3.1.nc'))
ind3=grep(paste0('X',yr,'0301'),names(s))
ind4=grep(paste0('X',yr,'0401'),names(s))
ind5=grep(paste0('X',yr,'0501'),names(s))
ind6=grep(paste0('X',yr,'0601'),names(s))
rfirst=s[[c(ind3,ind4,ind5,ind6)]]

snotelrecon=subset(snotelrecon,date==as.POSIXct('2010-04-01',tz='MST'))
coordinates(snotelrecon)=~Longitude+Latitude
f=extract(r,snotelrecon,sp=T)
ii=which(f$recon!=f$X20100401)
str(f[ii,])

dF=data.frame(reconf=snotelrecon2012$recon,f)
qplot(data=f@data,recon,X20100401)
snotelrecon2012

str(snotelrecon2012)

dF=data.frame(date=snotelrecon2012$date,fsca=f,Station_ID=snotelrecon2012$Station_ID)


huc4=readOGR('data/gis/','WBD_HU4_subset')
unique(huc4$HU_4_Name)

# dates=as.POSIXct(paste0(yr,'-',sprintf('%02d',seq(3,6)),'-01'),tz='MST')
dates=seq(from=as.POSIXct(paste0(yr,'-03-01'),tz='MST'),to=as.POSIXct(paste0(yr,'-06-01'),tz='MST'),by=3600*24)
s.sub=s[[1:length(dates)]]

load('data/recon_v3.1/2012basinfsca.rda') 
# --> basinfsca from >> firstbasinfsca=extract(s.sub,huc4)
# basinfsca=ldply(firstbasinfsca,function(x){
# 	x[x>1]=NA
# 	data.frame(fscaavg=apply(x,2,mean,na.rm=T),date=dates)
# })
# # colnames(basinfsca)=c('avg_fsca')
# basinfsca$basin=rep(unique(huc4$HU_4_Name),each=length(dates))

basinfsca$region=ifelse(as.character(basinfsca$basin)=='Upper Colorado-Dirty Devil','west',
			ifelse(as.character(basinfsca$basin)=='Lower Green','west',
			ifelse(as.character(basinfsca$basin)=='Escalante Desert-Sevier Lake','west',
			ifelse(as.character(basinfsca$basin)=='Great Salt Lake','west',
			ifelse(as.character(basinfsca$basin)=='Bear','west',
			ifelse(as.character(basinfsca$basin)=='Great Divide-Upper Green','west',
			ifelse(as.character(basinfsca$basin)=='Big Horn','west',
			ifelse(as.character(basinfsca$basin)=='Upper Snake','west',
			ifelse(as.character(basinfsca$basin)=='North Platte','east',
			ifelse(as.character(basinfsca$basin)=='White-Yampa','east',
			ifelse(as.character(basinfsca$basin)=='Colorado Headwaters','east',
			ifelse(as.character(basinfsca$basin)=='Upper Colorado-Dolores','east',
			ifelse(as.character(basinfsca$basin)=='San Juan','east',
			ifelse(as.character(basinfsca$basin)=='Rio Grande-Elephant Butte','east',
			ifelse(as.character(basinfsca$basin)=='Upper Pecos','east',
			ifelse(as.character(basinfsca$basin)=='Upper Canadian','east',
			ifelse(as.character(basinfsca$basin)=='Rio Grande Headwaters','east',
			ifelse(as.character(basinfsca$basin)=='Gunnison','east',
			ifelse(as.character(basinfsca$basin)=='Upper Arkansas','east',
			ifelse(as.character(basinfsca$basin)=='South Platte','east',
			ifelse(as.character(basinfsca$basin)=='Salt','south',
			ifelse(as.character(basinfsca$basin)=='Upper Gila','south',
			ifelse(as.character(basinfsca$basin)=='Little Colorado','south',
			ifelse(as.character(basinfsca$basin)=='Rio Grande Closed Basins','east',NA))))))))))))))))))))))))

gw=ggplot(subset(basinfsca,region=='west'))+
	geom_path(aes(x=date,y=fscaavg,colour=basin))+
	facet_wrap(~region)+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1))
ge=ggplot(subset(basinfsca,region=='east'))+
	geom_path(aes(x=date,y=fscaavg,colour=basin))+
	scale_colour_manual(values=c(brewer.pal(12,'Paired'),'black'))+
	facet_wrap(~region)+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1))
gs=ggplot(subset(basinfsca,region=='south'))+
	geom_path(aes(x=date,y=fscaavg,colour=basin))+
	facet_wrap(~region)+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1))
grid.arrange(gw,ge,gs)
g=arrangeGrob(gw,ge,gs,nrow=3)
ggsave(plot=g,file='data/recon_v3.1/huc4_east-west-south_basinfsca_timeseries2012.pdf',height=12,width=8)

basinfsca$dy=as.numeric(strftime(basinfsca$date,'%d'))
basinfscafirst=subset(basinfsca,dy==1)
ggplot(basinfscafirst)+
	geom_bar(aes(x=basin,y=fscaavg,fill=as.factor(date)),stat='identity',position='dodge')+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1))
ggsave('data/recon_v3.1/huc4_basinfsca_firstofmonth2012.pdf',height=6,width=12)

## if you ever extract all the years of data, use this to plot.
# ggplot(basinfsca)+
# 	geom_raster(aes(x=as.numeric(doy),y=strftime(date,'%Y'),fill=fscaavg)+
# 	facet_grid(region~basin)


#### raster plots of snotel timeseries
snotellocs.df=as.data.frame(snotellocs)
snotelrec$doy=as.numeric(strftime(snotelrec$date,'%j'))
apr1doy=data.frame(yr=seq(1999,2012))
apr1doy$doy=as.numeric(strftime(as.POSIXct(paste0(1999:2012,'-04-01',tz='MST')),'%j'))
snotelrec$swecut=cut(snotelrec$swe,c(seq(0,1,0.1),1.5,max(snotelrec$swe)),right=T)
l_ply(unique(snotelrec$Station_ID),function(x){
g=ggplot(subset(snotelrec,Station_ID==x))+
	geom_raster(aes(x=doy,y=as.factor(yr),fill=swecut),rev=T)+
	#geom_segment(data=apr1doy,aes(x=doy,xend=doy,y=as.factor(yr),yend=as.factor(yr),group=1))
	# geom_path(data=apr1doy,aes(x=doy,y=as.factor(yr),group=as.factor(yr)))#xend=doy,y=as.factor(yr),yend=as.factor(yr)))
	scale_x_continuous(breaks=c(1,seq(30,360,30)))+
	facet_wrap(~Station_ID)+
	scale_fill_manual(name='SWE [m]',
		values=c('grey98',brewer.pal(name='PuBu',n=9),'grey15','black'),
		drop=F,
		guide = guide_legend(reverse=TRUE))+
	theme(legend.key=element_rect(colour='grey90'))
ggsave(plot=g,file=paste0('data/snotelQAQC/plots/snotelrec_rasterlowvals_',x,'_',snotellocs.df[snotellocs.df$Station_ID==x,'Site_Name'],'_2000-2012.pdf'),width=12,height=6)
})

### raster plots of recon timeseries
recon.version='v3.2'
for(yr in 2000:2012){
load(paste0('data/recon_',recon.version,'/snotelrecondata_',yr,'_',recon.version,'.RData'))
}
snotelrecon=do.call(rbind,mget(ls(pattern='snotelrecon2')))
snotelrecon$doy=as.numeric(strftime(snotelrecon$date,'%j'))
snotelrecon$yr=strftime(snotelrecon$date,'%Y')
snotelrecon$reconcut=cut(snotelrecon$recon,c(seq(0,1,0.1),1.5,max(snotelrecon$recon)))
l_ply(unique(snotelrecon$Station_ID),function(x){
dF=subset(snotelrecon,Station_ID==x)
if(length(unique(dF$reconcut))>1){
gr=ggplot(dF)+
	geom_raster(aes(x=doy,y=as.factor(yr),fill=reconcut))+
	#geom_segment(data=apr1doy,aes(x=doy,xend=doy,y=as.factor(yr),yend=as.factor(yr),group=1))
	# geom_path(data=apr1doy,aes(x=doy,y=as.factor(yr),group=as.factor(yr)))#xend=doy,y=as.factor(yr),yend=as.factor(yr)))
	scale_x_continuous(limits=c(1,365),breaks=c(1,seq(30,360,30)))+
	facet_wrap(~Station_ID)+
	scale_fill_manual(name='SWE [m]',
		values=c('grey98',brewer.pal(name='PuBu',n=9),'grey15','black'),
		drop=FALSE,
		guide = guide_legend(reverse=TRUE))+
	theme(legend.key=element_rect(colour='grey90'))
ggsave(plot=gr,paste0('data/recon_',recon.version,'/snotelrecon_raster_',x,'_',snotellocs.df[snotellocs.df$Station_ID==x,'Site_Name'],'_2000-2012.pdf'),width=12,height=6)
}
})

#### raster plot of recon-snotel difference timeseries
recon.version='v3.2'
for(yr in 2000:2012){
load(paste0('data/recon_',recon.version,'/snotelrecondata_',yr,'_',recon.version,'.RData'))
}
snotelrecon=do.call(rbind,mget(ls(pattern='snotelrecon2')))
snotelrecon$doy=as.numeric(strftime(snotelrecon$date,'%j'))
snotelrecon$yr=strftime(snotelrecon$date,'%Y')
swe=merge(snotelrecon,snotelrec,by=c("Station_ID",'date','doy','yr'))
swe$diff=with(swe,recon-swe)
swe$diff[swe$diff==0]=NA
swe$diffcut=cut(swe$diff,seq(-1.25,1.25,.25))
l_ply(unique(swe$Station_ID),function(x){
dF=subset(swe,Station_ID==x)
if(length(unique(dF$diffcut))>1){
gd=ggplot(dF)+
	geom_raster(aes(x=doy,y=as.factor(yr),fill=diffcut))+
	scale_x_continuous(limits=c(1,365),breaks=c(1,seq(30,360,30)))+
	annotate('text',label='0 not shown\nonly occurs when both are 0\n+ -> recon greater\n- -> snotel greater',x=300,y=as.factor(2011))+
	facet_wrap(~Station_ID)+
	scale_fill_manual(name='SWE [m]',
		values=c(brewer.pal(name='PRGn',n=10)),
		drop=FALSE,
		guide = guide_legend(reverse=TRUE))+
	theme(legend.key=element_rect(colour='grey90'))
ggsave(plot=gd,paste0('data/recon_',recon.version,'/diffsnotelrecon_raster_',x,'_',snotellocs.df[snotellocs.df$Station_ID==x,'Site_Name'],'_2000-2012.pdf'),width=12,height=6)
}
})


### where is 2012 swe higher than 2011
snotelrec.sub=subset(snotelrec,yr==2012 | yr==2011)
maxswe=ddply(snotelrec.sub,.(yr,Station_ID),function(dF){
	summarise(dF,
		peakswe=max(dF$swe,na.rm=T))
})
peakswe=dcast(maxswe,Station_ID~yr)
peakswe$diff=peakswe$'2011'-peakswe$'2012'
staind=which(peakswe$'2011'<peakswe$'2012')
stationid=peakswe[staind,'Station_ID']
peakswe2=peakswe[staind,]
snotelinfo=snotellocs.df[snotellocs.df$Station_ID %in% stationid,]
snotelpeaks=merge(snotelinfo,peakswe2,by='Station_ID')

phvind=which(phvrec$Station_ID %in% snotelpeaks$Station_ID)
phvrecsub=phvrec[phvind,]
snotelpeaks=merge(snotelpeaks,phvrecsub,by='Station_ID')

write.table(snotelpeaks,file='data/snotelQAQC/peakswe_2012gt2011.txt',sep='\t',row.names=F,quote=F)

pdf('data/snotelQAQC/peakswe_2012gt2011_map.pdf')
plot(snotellocs)
coordinates(snotelpeaks)=~Longitude+Latitude
plot(snotelpeaks,add=T,col='red',pch=1)
dev.off()

snotelpeaks=as.data.frame(snotelpeaks)
ggplot(snotelpeaks)+
geom_bar(aes(x=NWbdiff,y=..density..),binwidth=200)
dev.set(2)
ggplot(snotelpeaks)+
geom_bar(aes(x=Lat,y=..density..),binwidth=2)

snotellocs.df=as.data.frame(snotellocs)
dev.new()
ggplot(phvrec)+
geom_bar(aes(x=NWbdiff,y=..density..),binwidth=200)

ggplot(phvrec)+
geom_bar(aes(x=Lat,y=..density..),binwidth=2)


