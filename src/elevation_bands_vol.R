library('ProjectTemplate')
# setwd('/Volumes/Dominik/Documents/snotel-regression_project')
load.project()
library(doMC)
library(scales)
library(dplyr)
library(tidyr)
# library(raster) 
# library(rgdal)
# library(plyr)
# library(reshape2)
# library(gtable)
# library(gridExtra)



# dem=raster('data/gis/gmted_uco_recongrid.tif')
# dem=dem*3.28083989501#meters to feet conversion
# dem[dem<4000]=3000
# dem[dem>=4000 & dem<5000]=4000
# dem[dem>=5000 & dem<6000]=5000
# dem[dem>=6000 & dem<7000]=6000
# dem[dem>=7000 & dem<8000]=7000
# dem[dem>=8000 & dem<9000]=8000
# dem[dem>=9000 & dem<10000]=9000
# dem[dem>=10000 & dem<11000]=10000
# dem[dem>=11000 & dem<12000]=11000
# dem[dem>=12000 & dem<13000]=12000
# dem[dem>=13000]=13000

dem.m=raster('data/gis/gmted_uco_recongrid.tif')
demm=dem.m
demm[demm<2000]=1000
demm[demm>=2000 & demm <2500]=2000
demm[demm>=2500 & demm<3000]=2500
demm[demm>=3000 & demm<3500]=3000
demm[demm>=3500 & demm<4000]=3500
demm[demm>=4000 & demm<4500]=4000

elev_area=freq(demm)
elev_area=mutate(as.data.frame(elev_area),
	area=count*500*500,
	area.km2=area/1000^2)#

### -- make upper CRB mask
# demm.uco=demm
# values(demm.uco)=demm
# tmp=extract(demm.uco,huc4)
# ucomask=demm
# values(ucomask)=NA
# ucomask[unlist(tmp)]=1
# plot(ucomask)
# writeRaster(ucomask,'data/gis/UpperCRB_mask.tif')


# ----Read in some gis files
ucomask=raster('data/gis/UpperCRB_mask.tif')
demm.uco=mask(demm,ucomask)
dem.m.uco=mask(dem.m,ucomask)
huc4=readOGR('data/gis','UpperCRB')
huc4raster=raster('data/gis/UpperCRB_rasterized.tif')
fordenstack=stack('data/forestdensity/umd_forden.nc')
slope.uco=raster('data/gis/gmted_uco_recongrid_aspect.tif')


units='metric'
#--- Regression Analysis ---
recon.version='v3.1'
covrange='idp1'
scalesnotel='scale'
product='phvrcn'
cost='r2'
fscaMatch='wofsca'
style='real-time'
residblending='unblended'
if(residblending=='blended'){
	resid='full'
} else {
	resid=''
}
dateflag='B'
config=''
graphbase=paste0('graphs/rswe_',recon.version,'/covrange',covrange,'/snotel',scalesnotel,'/',dateflag,'/',fscaMatch,'/',style,'/')


registerDoMC(3)
yr=2010
dte=as.POSIXct('1900-04-01',tz='MST')
dy=strftime(dte,'%d')
month=strftime(dte,'%B')
mnth=strftime(dte,'%b')
mth=strftime(dte,'%m')
for(snotelscale in c('scale')){#,'noscale'
for(fscaMatch in c('wofsca')){
for(resid in c('')){
for(product in c('phv','phvrcn')){
swestats_all=data.frame()
ucoelevzone_all=data.frame()
for (yr in seq(2001,2012) ) {#yr=2010#2011#2012
	predfile=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',snotelscale,'/',dateflag,'/fullpreds/',cost,'/netcdf/',style,'/',config,'/',resid,'preds-',product,'_',yr,'_blend-',fscaMatch,'.nc')
	s=stack(predfile)
	band=grep(paste0(yr,mth,dy),names(s))
	r=s[[band]]
	r[r==253]=NA
	r=mask(r,ucomask)
# plot(r)
	# projection(r)='+init=epsg:4326'
if(snotelscale=='noscale'){
	fscastack=stack(paste0('data/recon_',recon.version,'/reconfscadata_',yr,'_',recon.version,'.nc'))
	band=grep(paste0(yr,mth,dy),names(fscastack))
	fsca=fscastack[[band]]
	r=r*fsca
}

## - forest density
	layerid=paste0('X',yr-2011+12)
	layernum=grep(layerid,names(fordenstack))
	forden=fordenstack[[layernum]]
	# basinforden=as.data.frame(zonal(forden,huc4raster,'mean'))
	# basinforden=merge(as.data.frame(huc4[,c('HUC_4','HU_4_Name')]),basinforden,by.x='HUC_4',by.y='zone')
	# colnames(basinforden)=c('zone','basin','fordenavg')
	# basinforden$fordenavg=floor(basinforden$fordenavg*100)/100

x=unique(huc4$HU_4_Name)[3]
print(yr)
swestats=ldply(unique(huc4$HU_4_Name),.parallel=F,function(x){#don't parallelize, will have segfault
	basinpoly=huc4[grep(x,huc4$HU_4_Name),]
	print(basinpoly$HU_4_Name)
	r.ext=unlist(raster::extract(r,basinpoly))
	# zonal(r,huc4raster)
	if(units=='metric'){
		dem.ext=unlist(raster::extract(dem.m,basinpoly))
		dem.zone=unlist(raster::extract(demm,basinpoly))
		demcount=as.data.frame(table(dem.zone))
		forden.ext=unlist(raster::extract(forden,basinpoly))
		valdF=data.frame(basin=x,yr=yr,swe=r.ext,zone=dem.zone,elev=dem.ext,forden=forden.ext)
		swedF=ddply(valdF,.(basin,yr,zone),function(dF){
		# print(dF[1:40,])
		summarise(dF,
			elevavg=mean(elev,na.rm=T),#meters
			fordenavg=mean(forden,na.rm=T),
			sweavg=mean(swe,na.rm=T),#meters
			vol=sum(swe,na.rm=T)*500*500/10^9)#cubic km
		})
		swedF=mutate(swedF,
			areapct=demcount$Freq/sum(demcount$Freq),
			basinsweavg=sum(sweavg*areapct),
			basinswevol=basinsweavg*nrow(valdF)*500*500/10^9)
	} else {
	print('only metric supported')
	stop()
	}
	return(swedF)
})
swestats_all=rbind(swestats_all,swestats)

## get uco  elev zone stats
elevzone=as.data.frame(zonal(r,demm.uco,'mean'))
colnames(elevzone)=c('zone','sweavg')
elevzone2=as.data.frame(zonal(r,demm.uco,'sum'))
colnames(elevzone2)=c('zone','vol')
elevzone2=mutate(elevzone2,vol=vol*500*500/10^9)
ucoelevzone=data.frame(yr,merge(elevzone,elevzone2,by='zone'))
ucoelevzone_all=rbind(ucoelevzone_all,ucoelevzone)
}
write.table(format(swestats_all,digits=2,nsmall=0,scientific=F),paste0('reports/',style,'/',snotelscale,'/upperCRB_huc4_elevbands_',residblending,'swetable_',units,'_',dy,mnth,'-',product,'-',covrange,'-',cost,'-',fscaMatch,'.txt'),sep='\t',quote=F,row.names=F)
saveRDS(ucoelevzone_all,paste0('reports/',style,'/',snotelscale,'/upperCRB_elevbands_',residblending,'swetable_',units,'_',dy,mnth,'-',product,'-',covrange,'-',cost,'-',fscaMatch,'.rds'))
write.table(format(ucoelevzone_all,digits=2,nsmall=0,scientific=F),paste0('reports/',style,'/',snotelscale,'/upperCRB_elevbands_',residblending,'swetable_',units,'_',dy,mnth,'-',product,'-',covrange,'-',cost,'-',fscaMatch,'.txt'),sep='\t',quote=F,row.names=F)
saveRDS(swestats_all,paste0('reports/',style,'/',snotelscale,'/upperCRB_huc4_elevbands_',residblending,'swetable_',units,'_',dy,mnth,'-',product,'-',covrange,'-',cost,'-',fscaMatch,'.rds'))
}}}}



recon.version='v3.1'
covrange='idp1'
snotelscale='scale'
product='phvrcn'
cost='r2'
fscaMatch='wofsca'
style='real-time'
residblending='unblended'
if(residblending=='blended'){
	resid='full'
} else {
	resid=''
}
dateflag='B'
config=''
units='metric'

registerDoMC(3)
yr=2010
dte=as.POSIXct('1900-04-01',tz='MST')
dy=strftime(dte,'%d')
month=strftime(dte,'%B')
mnth=strftime(dte,'%b')
mth=strftime(dte,'%m')

##--- Get UCRB statistics based on elevation
product='phvrcn'
ucoelevzone_phvrcn=readRDS(paste0('reports/',style,'/',scalesnotel,'/upperCRB_elevbands_',residblending,'swetable_',units,'_',dy,mnth,'-',product,'-',covrange,'-',cost,'-',fscaMatch,'.rds'))
product='phv'
ucoelevzone_phv=readRDS(paste0('reports/',style,'/',scalesnotel,'/upperCRB_elevbands_',residblending,'swetable_',units,'_',dy,mnth,'-',product,'-',covrange,'-',cost,'-',fscaMatch,'.rds'))
ucoelevzone_phvrcn$model='PHVRCN'
ucoelevzone_phv$model='PHV'
ucoelevzone_all=bind_rows(ucoelevzone_phvrcn,ucoelevzone_phv)

format(ucoelevzone_all,scientific=F)
## Average SWE calculation ----
ucoelevzone_all %>%
     mutate(elevrange=ifelse(zone<3000,'low','high')) %>%
     group_by(model,yr,elevrange) %>%
     summarise(swe=mean(sweavg)) %>%
     spread(model,swe) %>%
     mutate(modeldiff=PHVRCN-PHV,
            modeldiffpct=modeldiff/PHV*100) %>%
     group_by(elevrange) %>%
     summarise(avgdiff=mean(modeldiff),
               avgdiffpct=mean(modeldiffpct))

productavg=dcast(ucoelevzone_all,yr+zone~model,value.var='sweavg')
yrswe=subset(productavg,yr==2011)
#2011
colMeans(yrswe[1:6,c('PHV','PHVRCN')])
low=colMeans(yrswe[1:3,c('PHV','PHVRCN')])
low
high=colMeans(yrswe[4:6,c('PHV','PHVRCN')])
high
low/(low+high)*100
high/(low+high)*100

ddply(yrswe,.(yr,zone),function(dF) {
	summarise(dF,
		zone,
		PHV.pct=PHV/sum(yrswe$PHV)*100,
	 	PHVRCN.pct=PHVRCN/sum(yrswe$PHVRCN)*100)
})

#2012
yrswe=subset(productavg,yr==2012)
colSums(yrswe[1:6,c('PHV','PHVRCN')])
low=colMeans(yrswe[1:3,c('PHV','PHVRCN')])
low
high=colMeans(yrswe[4:6,c('PHV','PHVRCN')])
high
low/(low+high)*100
high/(low+high)*100

ddply(yrswe,.(yr,zone),function(dF) {
	summarise(dF,
		zone,
		PHV.pct=PHV/sum(yrswe$PHV)*100,
	 	PHVRCN.pct=PHVRCN/sum(yrswe$PHVRCN)*100)
})



## Volume calculation ----
## yearly model differences
ucoelevzone_all %>%
     group_by(model,yr,zone) %>%
     summarise(swevol=sum(vol)) %>%
     spread(model,swevol) %>%
     mutate(modeldiff=PHVRCN-PHV,
            modeldiffpct=modeldiff/PHV*100) %>%
     group_by(zone) %>%
     summarise(sumphv=sum(PHV),
               sumphvrcn=sum(PHVRCN) ) %>%
     mutate(modeldiff=sumphvrcn-sumphv,
            modeldiffpct=modeldiff/sumphv*100) %>%
     summarise_each('mean')
     
     summarise(zone.phv=zone[which.max(PHV)],
               zone.phvrcn=zone[which.max(PHVRCN)])
     
ucoelevzone_all %>%
     mutate(elevrange=ifelse(zone<3000,'low','high')) %>%
     group_by(model,yr,elevrange) %>%
     summarise(swe=mean(vol)) %>%
     spread(model,swe) %>%
     mutate(modeldiff=PHVRCN-PHV,
            modeldiffpct=modeldiff/PHV*100) %>% 
     group_by(elevrange) %>%
     summarise(avgdiff=mean(modeldiff),
               avgdiffpct=mean(modeldiffpct))

productvol=dcast(ucoelevzone_all,yr+zone~model,value.var='vol')
volsum=subset(productvol,yr==2011)
#2011
colSums(volsum[1:6,c('PHV','PHVRCN')])
low=colSums(volsum[1:3,c('PHV','PHVRCN')])
diff(low)
high=colSums(volsum[4:6,c('PHV','PHVRCN')])
diff(high)
low/(low+high)*100
high/(low+high)*100
 
ddply(volsum,.(yr,zone),function(dF) {
	summarise(dF,
		zone,
		PHV.pct=PHV/sum(volsum$PHV)*100,
 	PHVRCN.pct=PHVRCN/sum(volsum$PHVRCN)*100)
})


#2012
productvol=dcast(ucoelevzone_all,yr+zone~model,value.var='vol')
volsum=subset(productvol,yr==2012)
colSums(volsum[1:6,c('PHV','PHVRCN')])
low=colSums(volsum[1:3,c('PHV','PHVRCN')])
diff(low)
high=colSums(volsum[4:6,c('PHV','PHVRCN')])
diff(high)
low/(low+high)*100
high/(low+high)*100

ddply(volsum,.(yr,zone),function(dF) {
	summarise(dF,
		zone,
		PHV.pct=PHV/sum(volsum$PHV)*100,
 	PHVRCN.pct=PHVRCN/sum(volsum$PHVRCN)*100)
})



## --- Get UCRB information based on huc4 subbasins and elevation
product='phvrcn'
swestats_phvrcn=readRDS(paste0('reports/',style,'/',scalesnotel,'/upperCRB_huc4_elevbands_',residblending,'swetable_',units,'_',dy,mnth,'-',product,'-',covrange,'-',cost,'-',fscaMatch,'.rds'))
product='phv'
swestats_phv=readRDS(paste0('reports/',style,'/',scalesnotel,'/upperCRB_huc4_elevbands_',residblending,'swetable_',units,'_',dy,mnth,'-',product,'-',covrange,'-',cost,'-',fscaMatch,'.rds'))
swestats_phvrcn$model='PHVRCN'
swestats_phv$model='PHV'
swestats_all=bind_rows(swestats_phvrcn,swestats_phv)
str(swestats_all)
mypal=c(brewer.pal(9,name='Set1'),'brown','grey35','black')
## plot facet by product and basin, elevation band on x axis, all years as colors
g=ggplot(swestats_all)+
	geom_bar(aes(x=zone,y=vol,fill=as.factor(yr)),stat='identity',position='dodge')+
	facet_grid(basin~model,scales='free_y')+
	# scale_y_continuous(label=comma)+
	scale_fill_manual(values=mypal)+
	theme_bw()+
	labs(title=paste(covrange,snotelscale,sep=' - '), y=expression(SWE~Volume~group("[",km^3,"]")))

ggsave(plot=g,file=paste0('reports/',style,'/',snotelscale,'/upperCRB_huc4_elevbands_',residblending,'swetable_',dy,mnth,'-',units,'-',covrange,'-',cost,'-',fscaMatch,'.pdf'),width=17,height=6)



## plot facet 
plot_labeller <- function(variable,value){
     if (variable=='basin') {
          return(basin_names[value])
     } else if (variable=='zone') {
          return(zone_names[value])
     } else {
          return(as.character(value))
     }
}
basin_names <- list(
     'Colorado Headwaters'="Colorado\nHeadwaters",
     'Great Divide-Upper Green'='Great Divide\nUpper Green',
     'Gunnison'='Gunnison',
     'Lower Green'='Lower Green',
     'San Juan'='San Juan',
     'Upper Colorado-Dirty Devil'='U.Colorado\nDirty Devil',
     'Upper Colorado-Dolores'='U.Colorado\nDolores',
     'White-Yampa'='White-Yampa'
)
zone_names <- list(
     '4000'='4000+ m' ,
     '3500'='[3500,4000)',
     '3000'='[3000,3500)',
     '2500'='[2500,3000)',
     '2000'='[2000,2500)',
     '1000'='<2000 m'
)

swestats_all=swestats_all %>%
     mutate(volperc=vol/basinswevol)
filter(swestats_all,yr==2011 | yr==2012) %>%
     {
          ggplot(.)+
          geom_path(aes(zone,volperc,colour=model,linetype=model),size=1)+
          scale_colour_manual(name='Model',values=c('blue','red'),drop=F,labels=c('PHV-baseline','PHV-RCN'))+
          scale_linetype_manual(name='Model', values=c(1,2),labels=c('PHV-baseline','PHV-RCN'))+
          labs(y='fraction of HUC4 Volume [/100]',x='Elevation [m]')+
          facet_grid(yr~basin,labeller=plot_labeller)+
          theme_bw(base_size=12)+
          theme(strip.text=element_text(face='bold'),
                panel.grid=element_blank(),
                strip.background=element_blank(),
                axis.line=element_line(),
                legend.direction='vertical',
                legend.position=c(-0.005,.89),
                legend.justification='left',
                legend.key.height=unit(.04,'npc'),
                legend.key=element_rect(colour=NA),
                legend.key.width=unit(.03,'npc'))
          }
ggsave(filename=paste0(graphbase,'/upperCRB_huc4_elevbands_',residblending,'_percentswevol_',dy,mnth,'-',units,'.pdf'),width=300,height=100,units='mm')



## huc4/elevation band swe volume bar graph ----
statmin=dcast(swestats_all,zone+basin~model,value.var='vol',fun.aggregate=min)
statmax=dcast(swestats_all,zone+basin~model,value.var='vol',fun.aggregate=max)

swevolrange=ddply(swestats_all,.(zone,basin,model),function(x){
	summarise(x,
		min=min(vol),
		max=max(vol))
})
swevolrange$model=ifelse(swevolrange$model=='PHV','PHV-baseline','PHV-RCN')
swevolrange$zone=factor(swevolrange$zone,levels=rev(unique(swevolrange$zone)))
#
swestatsavg=dcast(swestats_all,zone+basin~model,value.var='vol',fun.aggregate=mean)
swestatsavg.plot=gather(swestatsavg,model,vol,PHV:PHVRCN)
swestatsavg.plot$zone=factor(swestatsavg.plot$zone,levels=rev(unique(swestatsavg.plot$zone)))
swestatsavg.plot$model=ifelse(swestatsavg.plot$model=='PHV','PHV-baseline','PHV-RCN')
swestatsavg.plot$vol=ifelse(swestatsavg.plot$vol==0,0.001,swestatsavg.plot$vol)

gg=ggplot(swestatsavg.plot)+
	geom_bar(aes(x=as.factor(model),y=vol,fill=as.factor(model)),stat='identity',width=1,position=position_dodge(width=0))+
	scale_fill_manual(name='Model',values=c('blue','red'),drop=F)+
	geom_errorbar(data=swevolrange,aes(x=model,ymin=min,ymax=max),width=.45,size=0.5)+
	coord_trans(y='sqrt')+#,limy=c(0,12))+
	scale_x_discrete(labels=c('PHV\nbaseline','PHV\nRCN'),expand=c(.05,0))+
	labs(x='Model',y=expression(SWE~Volume~group("[",km^3,"]")))+
	facet_grid(zone~basin,scales='free_y',labeller=plot_labeller)+#plot_labeller)+
	theme_minimal(base_size=8)+
	theme(strip.text=element_text(face='bold'),
		panel.grid=element_blank(),
		strip.background=element_blank(),
		axis.line=element_line(),
		legend.direction='vertical',
		legend.position=c(.85,.95),
		legend.justification='left',
		legend.key.size=unit(.5,'lines'))
 # show(gg)
ggsave(gg,file=paste0(graphbase,'/upperCRB_huc4_elevbands_',residblending,'_swevol_',dy,mnth,'-',units,'.pdf'),width=6.6,height=6)
    




ucosum=dcast(swestats_all,yr+basin~model,value.var='vol',fun.agg=sum)
yrlyucosum=ddply(ucosum,.(yr),summarise,phv=sum(PHV),phvrcn=sum(PHVRCN))
colMeans(yrlyucosum)
(47.0-45.4)/47
max(yrlyucosum$phv)-min(yrlyucosum$phv)
max(yrlyucosum$phvrcn)-min(yrlyucosum$phvrcn)

ucosum.m=melt(ucosum,id.vars=.(basin,yr))
ucosumrange=ddply(ucosum.m,.(basin,variable),function(x){
	summarise(x,
		avg=mean(value),
		min=min(value),
		max=max(value))
})
ucosumrange=mutate(ucosumrange,
	diff=max-min)
ucosumrange_diff=dcast(ucosumrange,basin~variable,value.var='diff')
colMeans(ucosumrange_diff[,c(2,3)])

### difference in average swe by basin
voldiff=mutate(dcast(ucosumrange,basin~variable,value.var='avg'),
	diff=PHVRCN-PHV,
	diffpct=diff/PHV*100
	)
voldiff
colMeans(voldiff[2:ncol(voldiff)])
mean(voldiff[-c(1,5),ncol(voldiff)])

format(swestats_all,digits=2,scientific=F)


elevsum=dcast(swestats_all,yr+zone~model,value.var='vol',fun.agg=mean)
#2011
colMeans(elevsum[1:3,c('PHV','PHVRCN')])
colMeans(elevsum[4:6,c('PHV','PHVRCN')])
#2012
colMeans(elevsum[7:9,c('PHV','PHVRCN')])
colMeans(elevsum[10:12,c('PHV','PHVRCN')])


volsum=dcast(swestats_all,yr+zone~model,value.var='vol',fun.agg=mean)
#2011
diff(colMeans(volsum[1:3,c('PHV','PHVRCN')]))
diff(colMeans(volsum[4:6,c('PHV','PHVRCN')]))
#2012
diff(colMeans(volsum[7:9,c('PHV','PHVRCN')]))
diff(colMeans(volsum[10:12,c('PHV','PHVRCN')]))

## huc4/swe avg bar graph ----
statmin=dcast(swestats_all,zone+basin~model,value.var='sweavg',fun.aggregate=min,na.rm=TRUE)
statmax=dcast(swestats_all,zone+basin~model,value.var='sweavg',fun.aggregate=max,na.rm=TRUE)

swevolrange=ddply(swestats_all,.(zone,basin,model),function(x){
     summarise(x,
               min=min(sweavg),
               max=max(sweavg))
})
swevolrange$model=ifelse(swevolrange$model=='PHV','PHV-baseline','PHV-RCN')
swevolrange$zone=factor(swevolrange$zone,levels=rev(unique(swevolrange$zone)))
#
swestatsavg=dcast(swestats_all,zone+basin~model,value.var='sweavg',fun.aggregate=mean)
swestatsavg.plot=gather(swestatsavg,model,sweavg,PHV:PHVRCN)
swestatsavg.plot$zone=factor(swestatsavg.plot$zone,levels=rev(unique(swestatsavg.plot$zone)))
swestatsavg.plot$model=ifelse(swestatsavg.plot$model=='PHV','PHV-baseline','PHV-RCN')
swestatsavg.plot$sweavg=ifelse(swestatsavg.plot$sweavg==0,0.001,swestatsavg.plot$sweavg)

gg=ggplot(swestatsavg.plot)+
     geom_bar(aes(x=as.factor(model),y=sweavg,fill=as.factor(model)),stat='identity',width=1,position=position_dodge(width=0))+
     scale_fill_manual(name='Model',values=c('blue','red'),drop=F)+
     geom_errorbar(data=swevolrange,aes(x=model,ymin=min,ymax=max),width=.3,size=0.5)+
     # coord_trans(y='sqrt')+#,limy=c(0,12))+
     scale_x_discrete(labels=c('PHV\nbaseline','PHV\nRCN'),expand=c(.05,0))+
     labs(x='Model',y="Average SWE Depth [m]")+
     facet_grid(zone~basin,scales='free_y',labeller=plot_labeller)+#plot_labeller)+
     theme_minimal(base_size=8)+
     theme(strip.text=element_text(face='bold'),
           panel.grid=element_blank(),
           strip.background=element_blank(),
           axis.line=element_line(),
           legend.direction='vertical',
           legend.position=c(.85,.95),
           legend.justification='left',
           legend.key.size=unit(.5,'lines'))
# show(gg)
ggsave(gg,file=paste0(graphbase,'/upperCRB_huc4_elevbands_',residblending,'_sweavg_',dy,mnth,'-',units,'.pdf'),width=6.6,height=6)



## Look at Gunnison basin in 2011 to compare to fassnacht ----
gunn=subset(swestats_all,basin=='Gunnison' & yr==2011)
gunn=dcast(gunn,yr+zone~model,value.var='vol')
colSums(gunn)
gunn$MVR=c(0,0.0014,0.0061,0.0070,0.0043,0)
fassnachtpct=ddply(gunn,.(yr), function(dF){
	summarise(dF,
		zone,
		PHV.pct=PHV/sum(PHV)*100,
 		PHVRCN.pct=PHVRCN/sum(PHVRCN)*100,
 		MVR.pct=MVR/sum(MVR)*100)
})
fassnachtpct
colSums(fassnachtpct[1:3,])
colSums(fassnachtpct[4:6,])




# ### --- snotel elevation analysis
# apr1snotel=subset(snotelrec,mth==4 & dy==1 & yr>=2010)
# apr1snotel=merge(apr1snotel,as.data.frame(snotellocs),by='Station_ID')
# coordinates(apr1snotel)=~Longitude+Latitude
# proj4string(apr1snotel)=CRS(proj4string(huc6))
# huc6name=huc6[,'HU_6_Name']
# snotelbasin=over(apr1snotel,huc6name)
# snotelelev=extract(dem,apr1snotel)
# apr1swe=as.data.frame(apr1snotel)
# apr1swe$basin=snotelbasin$HU_6_Name
# apr1swe$demelev=snotelelev
# apr1swe=apr1swe[!is.na(apr1swe$basin),c('Station_ID','date','basin','yr','swe','demelev')]

# snotel_avg_swe=ddply(apr1swe,.(basin,yr,demelev),function(x){
# 	summarise(x,
# 		avg=mean(swe,na.rm=T)*100/2.54,
# 		n=length(swe))
# })

# write.table(format(snotel_avg_swe,digits=2,nsmall=0,scientific=F),paste0('reports/wwa_huc6_swetable_1Apr-snotel.txt'),sep='\t',quote=F,row.names=F)

# ggplot(snotel_avg_swe)+
# 	geom_bar(aes(x=as.factor(demelev),y=avg,fill=as.factor(yr)),stat='identity',position='dodge')+
# 	geom_text(data=subset(snotel_avg_swe,yr=unique(yr)[1]),aes(x=as.factor(demelev),y=max(avg),label=paste0('n=',n)),vjust=0)+
# 	facet_wrap(~basin)+
# 	scale_y_continuous(limits=c(0,75),breaks=seq(0,60,10))+
# 	scale_x_discrete(limits=as.character(seq(4000,13000,1000)),labels=as.character(seq(4000,13000,1000)))+
# 	# coord_fixed(ratio=0.05)+
# 	labs(title='snotel', y='avg swe [inches]')+
# 	theme_bw()

# ggsave(paste0('reports/wwa_huc6_swetable_1Apr-snotel.pdf'),width=17,height=6) 