library('ProjectTemplate')
setwd('~/Documents/snotel-regression_project')
load.project()

recon.version='v3.1'
covrange='idp1'
snotelscale='scale'
cost='r2'
fscaMatch='fsca'
product='phv'
blend='full'
config='bic-v2.2-nosnoteltransform'
for(product in c('phv','phvrcn')){
	for(blend in c('','full')){
		for(snotelscale in c('scale')){
bn=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',snotelscale,'/fullpreds/',cost,'/netcdf/',config,'/')

get_regsurf=function(yr){
	f2r=paste0(bn,blend,'preds-',product,'_1April',yr,'-',fscaMatch,'.tif')
	r=raster(f2r)
# r=shift(r,0.00416666667/2,-0.0041666667/2)
	projection(r)='+init=epsg:4326'
	r[r>200]=NA
	return(r)
}

s=stack()
for(yr in 2010:2012){
r=get_regsurf(yr)
s=addLayer(s,r)
}

sntl=extract(s,snotellocs,sp=T)
sntl.df=as.data.frame(sntl)
#str(sntl.df)
ind=sntl.df$Station_ID %in% unique(snotelrec$Station_ID)
sntl.df=sntl.df[ind,]
colnames(sntl.df)[10:12]=c('Apr12010','Apr12011','Apr12012')
sntl.df=arrange(sntl.df,Station_ID)
# str(sntl.df)
write.csv(sntl.df,paste0(bn,'snotelvalues_',product,blend,'_',recon.version,'_',covrange,'_',snotelscale,'_',cost,'_',fscaMatch,'.txt'),row.names=F,quote=F)
}
	}
}




ss=subset(snotelrec,mth==4 & dy==1 & (yr==2010 | yr==2011 | yr==2012))
sscast=dcast(ss[,c('yr','mth','dy','Station_ID','swe','date')],Station_ID~date,value.var='swe')
snotellocs.df=as.data.frame(snotellocs)
ind=sscast$Station_ID %in% snotellocs$Station_ID
sscast=merge(sscast[ind,],snotellocs.df,by='Station_ID')
sntl=merge(sscast,phvrec,by='Station_ID')
str(sntl)
write.csv(sntl,paste0('diagnostics/snotelobs.txt'),row.names=F,quote=F)

## ---

snotellocs.df=as.data.frame(snotellocs)

ggplot(ucophv)+
geom_raster(aes(x=Long,y=Lat,fill=Long))

latrast=num2ucoRaster(ucophv$Lat)
longrast=num2ucoRaster(ucophv$Long)
coordstack=stack(longrast,latrast)
tmp=data.frame(coordinates(snotellocs),extract(coordstack,snotellocs))
tmp=mutate(tmp,
	latdiff=Latitude-layer.2,
	longdiff=Longitude-layer.1)
tmp$latflag2=ifelse(tmp$latdiff>15/3600/2,1,0)
tmp$longflag2=ifelse(tmp$longdiff>15/3600/2,1,0)
ind=which(tmp$latflag2==1)
tmp[ind,]

tmp=data.frame(coordinates(snotellocs.sub),extract(coords,snotellocs.sub))
tmp=mutate(tmp,
	latdiff=Latitude-layer.2,
	longdiff=Longitude-layer.1)
tmp$latflag2=ifelse(tmp$latdiff>15/3600/2,1,0)
tmp$longflag2=ifelse(tmp$longdiff>15/3600/2,1,0)
ind=which(tmp$latflag2==1)


# tmp=data.frame(coordinates(snotellocs.sub),extract(coords,snotellocs.sub))
coord.df=as.data.frame(coords)
names(coord.df)=c('Long','Lat')
ggplot()+
geom_raster(data=coord.df,aes(x=Long,y=Lat,fill=rnorm(nrow(coord.df))))+
geom_point(data=tmp2,aes(x=Longitude,y=Latitude),shape=1)+
geom_point(data=tmp2,aes(x=layer.1,y=layer.2),shape=3)+
# coord_cartesian(ylim=c(35.255,35.265),xlim=c(-112.065,-112.055))
coord_cartesian(ylim=c(35.3325,35.345),xlim=c(-111.66,-111.64))



### plot PHV variables
range(ucophv$Wd2ocean)
ucophv$nwbdiffcut=cut(ucophv$NWbdiff,breaks=c(min(ucophv$NWbdiff),seq(0,4100,500)),dig.lab=4)
ucophv$wbdiffcut=cut(ucophv$Wbdiff,breaks=seq(0,max(ucophv$Wbdiff),500),dig.lab=4)
ucophv$wdocut=cut(ucophv$Wd2ocean,breaks=c(min(ucophv$Wd2ocean),seq(750,1750,250)),dig.lab=4)
ucophv$latcut=cut_interval(ucophv$Lat,9,dig.lab=4)
ucophv$elevcut=cut_interval(ucophv$Elev,9,dig.lab=4)

snotellocs.df=as.data.frame(snotellocs)
ggplot(ucophv)+
	geom_raster(aes(x=Long,y=Lat,fill=nwbdiffcut))+
	geom_point(data=snotellocs.df,aes(x=Longitude,y=Latitude))+
	scale_fill_brewer(palette='BuPu')

ggsave('data/gis/NWbdiff_withsnotel.png')

ggplot(ucophv)+
	geom_raster(aes(x=Long,y=Lat,fill=wbdiffcut))+
	geom_point(data=snotellocs.df,aes(x=Longitude,y=Latitude))+
	scale_fill_brewer(palette='BuPu')

ggsave('data/gis/Wbdiff_withsnotel.png')

ggplot(ucophv)+
	geom_raster(aes(x=Long,y=Lat,fill=wdocut))+
	geom_point(data=snotellocs.df,aes(x=Longitude,y=Latitude))+
	scale_fill_brewer(palette='BuPu')

ggsave('data/gis/Wd2ocean_withsnotel.png')

ggplot(ucophv)+
	geom_raster(aes(x=Long,y=Lat,fill=latcut))+
	geom_point(data=snotellocs.df,aes(x=Longitude,y=Latitude))+
	scale_fill_brewer(palette='BuPu')

ggsave('data/gis/Lat_withsnotel.png')

ggplot(ucophv)+
	geom_raster(aes(x=Long,y=Lat,fill=elevcut))+
	geom_point(data=snotellocs.df,aes(x=Longitude,y=Latitude))+
	scale_fill_brewer(palette='BuPu')

ggsave('data/gis/Elev_withsnotel.png')


ucophv=mutate(ucophv,poly=Lat^2*Long^2)
range(ucophv$poly)
ucophv$polycut=cut(ucophv$poly,breaks=seq(min(ucophv$poly),max(ucophv$poly),5000))
ggplot(ucophv)+
	geom_raster(aes(x=Long,y=Lat,fill=polycut))+
	geom_point(data=snotellocs.df,aes(x=Longitude,y=Latitude))+
	scale_fill_brewer(palette='BuPu')

