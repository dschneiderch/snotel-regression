library('ProjectTemplate')
setwd('~/Documents/snotel-regression_project')
load.project()
# library(raster) 
# library(rgdal)
# library(plyr)
# library(reshape2)
# library(gtable)
# library(gridExtra)

recon.version='v3.1'
covrange='idp1'
cost='r2'
product='phvrcn'
fscaMatch='wofsca'
dateflag='B'
resid=''
style='real-time'
config=''#bic-v2.5-removed_NWbdiff'#'bic-v2.2-nosnoteltransform'
# yr=2010
# mth=4
dte=as.POSIXct('1900-06-01',tz='MST')
for(snotelscale in c('noscale','scale')){
for(product in c('phv','phvrcn','reconrt')){#
for(fscaMatch in c('wofsca','fsca')){
for(resid in c('full','')){
for(yr in seq(2001,2012)){
	for(mth in strftime(dte,'%m')){
		month=strftime(dte,'%B')
		for(dy in strftime(dte,'%d')){
			predfile=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',snotelscale,'/',dateflag,'/fullpreds/',cost,'/netcdf/',style,'/',config,'/',resid,'preds-',product,'_',yr,'_blend-',fscaMatch,'.nc')
			s=stack(predfile)
			band=grep(paste0(yr,mth,dy),names(s))
			file2eval=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',snotelscale,'/',dateflag,'/fullpreds/',cost,'/netcdf/',style,'/',config,'/',resid,'preds-',product,'_',dy,month,yr,'-',fscaMatch,'.tif')	
			cmd1=paste0('gdal_translate -b ',band, ' -a_nodata 253 ',predfile,' ',file2eval)
			system(cmd1)
		}
	}
}
}
}
}}
# gmted=raster('data/gis/gmted_westernUS.tif')
# dem=projectRaster(gmted,r)
# writeRaster(dem,'data/gis/gmted_uco_recongrid.tif',overwrite=T)

demm=raster('data/gis/gmted_uco_recongrid.tif')
dem=demm*3.28083989501#meters to feet conversion
dem[dem<4000]=3000
dem[dem>=4000 & dem<5000]=4000
dem[dem>=5000 & dem<6000]=5000
dem[dem>=6000 & dem<7000]=6000
dem[dem>=7000 & dem<8000]=7000
dem[dem>=8000 & dem<9000]=8000
dem[dem>=9000 & dem<10000]=9000
dem[dem>=10000 & dem<11000]=10000
dem[dem>=11000 & dem<12000]=11000
dem[dem>=12000 & dem<13000]=12000
dem[dem>=13000]=13000
#
demm[demm<2000]=1000
demm[demm>=2000 & demm <2500]=2000
demm[demm>=2500 & demm<3000]=2500
demm[demm>=3000 & demm<3500]=3000
demm[demm>=3500 & demm<4000]=3500
demm[demm>=4000 & demm<4500]=4000

elev_area=freq(dem)
elev_area=mutate(as.data.frame(elev_area),
	area=count*500*500*0.000247105)#convert meters sq to acres


# ----------  
# r2=raster('/Volumes/shareProjects/WWA/data/SWE_SNODIS/blended/SWEphvrcn20100401.tif')
# r2
# r2=crop(r2,dem)
# setMinMax(r2)

# recon=raster('/Volumes/shareProjects/WWA/data/SWE_SNODIS/vs3_1/01APR2010.tif')
# huc4=readOGR('/Volumes/shareProjects/WWA/data/hydro/TriState','WBD_HU4_subset')
# coh=huc4[grep('Colorado Headwaters',huc4$HU_4_Name),]
# ugn=huc4[grep('Great Divide-Upper Green',huc4$HU_4_Name),]
# r6=raster('/Volumes/shareProjects/WWA/data/hydro/TriState/HUC6_3_wtshd_bands.tif')
# r.prj=projectRaster(r,r6)

huc6=readOGR('data/gis','HUC6_3WWA_wtshd')
units='imp'
#--- Regression Analysis ---
recon.version='v3.1'
covrange='idp1'
snotelscale='scale'
product='phvrcn'
cost='r2'
fscaMatch='wofsca'
resid=''
dateflag='B'

yr=2010
dte=as.POSIXct('1900-04-01',tz='MST')
dy=strftime(dte,'%d')
month=strftime(dte,'%B')
mnth=strftime(dte,'%b')
for(snotelscale in c('scale')){#,'noscale'
for(fscaMatch in c('wofsca')){
for(resid in c('full','')){
for(product in c('phv','phvrcn')){
swestats_all=data.frame()
for (yr in c(2010,2011,2012) ) {#yr=2010#2011#2012
f2r=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',snotelscale,'/',dateflag,'/fullpreds/',cost,'/netcdf/',config,'/',resid,'preds-',product,'_',dy,month,yr,'-',fscaMatch,'.tif')
r=raster(f2r)
# plot(r)
# r=shift(r,0.00416666667/2,-0.0041666667/2)
projection(r)='+init=epsg:4326'
r[r>200]=NA
if(snotelscale=='scale'){
	fscastack=stack(paste0('data/recon_',recon.version,'/reconfscadata_',yr,'_',recon.version,'.nc'))
	band=grep(paste0(yr,mth,dy),names(fscastack))
	fsca=fscastack[[band]]
	r=r*fsca
}

# setMinMax(r)
swestats=ldply(unique(huc6$HU_6_Name1),function(x){
	basinpoly=huc6[grep(x,huc6$HU_6_Name1),]
	r.ext=unlist(extract(r,basinpoly))
	r.ext[r.ext==253]=NA
	if(units=='metric'){
		dem.ext=unlist(extract(demm,basinpoly))
		valdF=data.frame(swe=r.ext,elev=dem.ext)
		swedF=ddply(valdF,.(elev),function(dF){
		summarise(dF,
			avgswe=mean(swe,na.rm=T),#meters
			vol=avgswe*nrow(dF)*500*500/10^9)# cu. km
		})
	} else {
		dem.ext=unlist(extract(dem,basinpoly))
		valdF=data.frame(swe=r.ext,elev=dem.ext)
		swedF=ddply(valdF,.(elev),function(dF){
		summarise(dF,
			avgswe=mean(swe,na.rm=T)*100/2.54,#inches
			vol=avgswe*100/2.54*nrow(dF)*500*500*0.000247105)#acre-ft
		})
	}
	swedF$basin=x
	swedF=swedF[,c('basin','elev','avgswe','vol')]
	return(swedF)
})
swestats$yr=yr
write.table(format(swestats,digits=2,nsmall=0,scientific=F),paste0('reports/',snotelscale,'/wwa_huc6_',resid,'swetable_',units,'_',dy,mnth,yr,'-',product,'-',covrange,'-',cost,'-',fscaMatch,'.txt'),sep='\t',quote=F,row.names=F)
swestats_all=rbind(swestats_all,swestats)
}

if(units=='metric'){
	g=ggplot(swestats_all)+
		geom_bar(aes(x=as.factor(elev),y=vol,fill=as.factor(yr)),stat='identity',position='dodge')+
		facet_wrap(~basin)+
		theme_bw()+
		labs(title=paste(product,covrange,snotelscale,sep=' - '), y='SWE Volume [cu. km]')
} else {
	g=ggplot(swestats_all)+
		geom_bar(aes(x=as.factor(elev),y=avgswe,fill=as.factor(yr)),stat='identity',position='dodge')+
		facet_wrap(~basin)+
		theme_bw()+
	     	coord_cartesian(ylim=c(0,75))+
	     	labs(title=paste(product,covrange,snotelscale,sep=' - '), y='Avg SWE [inches]')
}
	
ggsave(plot=g,file=paste0('reports/',snotelscale,'/wwa_huc6_',resid,'swetable_',dy,mnth,'-',units,'-',product,'-',covrange,'-',cost,'-',fscaMatch,'.pdf'),width=17,height=6)
}
}
}
}





## --- Recon analysis ---
product='recon'
swestats_all=data.frame()
for (yr in c(2010,2011,2012) ) {#yr=2010#2011#2012

f2r=paste0('/Volumes/shareProjects/WWA/data/SWE_SNODIS/vs3_1/01APR',yr,'.tif')
r=raster(f2r)
projection(r)='+init=epsg:4326'
# r[r==253]=NA
# setMinMax(r)

swestats=ldply(unique(huc6$HU_6_Name1),function(x){
	basinpoly=huc6[grep(x,huc6$HU_6_Name1),]
	r.ext=unlist(extract(r,basinpoly))
	r.ext[r.ext==253]=NA
	dem.ext=unlist(extract(dem,basinpoly))
	valdF=data.frame(swe=r.ext,elev=dem.ext)
	swedF=ddply(valdF,.(elev),function(dF){
		summarise(dF,
			avg=mean(swe,na.rm=T)*100/2.54,#meters to inches
			sum=sum(swe,na.rm=T)*100/2.54)
	})
	swedF$basin=x
	swedF=swedF[,c(ncol(swedF),1,2)]
	return(swedF)
})
swestats$yr=yr
# write.table(format(swestats,digits=2,nsmall=0,scientific=F),paste0('reports/wwa_huc6_swetable_1Apr',yr,'-',product,'.txt'),sep='\t',quote=F,row.names=F)
swestats_all=rbind(swestats_all,swestats)
}
ggplot(swestats_all)+
	geom_bar(aes(x=as.factor(elev),y=avg,fill=as.factor(yr)),stat='identity',position='dodge')+
	facet_wrap(~basin)+
	scale_y_continuous(limits=c(0,205),breaks=seq(0,205,20))+
	# coord_fixed(ratio=0.05)+
	labs(title=paste(product,sep='\n'), y='avg swe [inches]')+
	theme_bw()
ggsave(paste0('reports/wwa_huc6_swetable_1Apr-',product,'.pdf'),width=17,height=6)




### --- snotel elevation analysis
apr1snotel=subset(snotelrec,mth==4 & dy==1 & yr>=2010)
apr1snotel=merge(apr1snotel,as.data.frame(snotellocs),by='Station_ID')
coordinates(apr1snotel)=~Longitude+Latitude
proj4string(apr1snotel)=CRS(proj4string(huc6))
huc6name=huc6[,'HU_6_Name']
snotelbasin=over(apr1snotel,huc6name)
snotelelev=extract(dem,apr1snotel)
apr1swe=as.data.frame(apr1snotel)
apr1swe$basin=snotelbasin$HU_6_Name
apr1swe$demelev=snotelelev
apr1swe=apr1swe[!is.na(apr1swe$basin),c('Station_ID','date','basin','yr','swe','demelev')]

snotel_avg_swe=ddply(apr1swe,.(basin,yr,demelev),function(x){
	summarise(x,
		avg=mean(swe,na.rm=T)*100/2.54,
		n=length(swe))
})

write.table(format(snotel_avg_swe,digits=2,nsmall=0,scientific=F),paste0('reports/wwa_huc6_swetable_1Apr-snotel.txt'),sep='\t',quote=F,row.names=F)

ggplot(snotel_avg_swe)+
	geom_bar(aes(x=as.factor(demelev),y=avg,fill=as.factor(yr)),stat='identity',position='dodge')+
	geom_text(data=subset(snotel_avg_swe,yr=unique(yr)[1]),aes(x=as.factor(demelev),y=max(avg),label=paste0('n=',n)),vjust=0)+
	facet_wrap(~basin)+
	scale_y_continuous(limits=c(0,75),breaks=seq(0,60,10))+
	scale_x_discrete(limits=as.character(seq(4000,13000,1000)),labels=as.character(seq(4000,13000,1000)))+
	# coord_fixed(ratio=0.05)+
	labs(title='snotel', y='avg swe [inches]')+
	theme_bw()

ggsave(paste0('reports/wwa_huc6_swetable_1Apr-snotel.pdf'),width=17,height=6) 