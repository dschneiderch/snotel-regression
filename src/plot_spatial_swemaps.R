setwd('/Volumes/Dominik/Documents/snotel-regression_project')
library('ProjectTemplate')
load.project()
library(dplyr)

#### plot swe map and swe difference maps
get_regsurf=function(yr,product,sample=NULL){
path=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',scalesnotel,'/',dateflag,'/fullpreds/',cost,'/netcdf/',style,'/',config,'/')
f2r=paste0(resid,'preds-',product,'_',yr,'_blend-',fscaMatch,'.nc')
s=stack(file.path(path,f2r))
layerid=grep(paste0(yr,'0401'),names(s))
# projection(r)='+init=epsg:4326'
r=s[[layerid]]
r[r>200]=NA
if(is.null(sample)){
	return(as.data.frame(r,xy=T))
# str(sampleRegular(r,1e6))
} else {
	return(raster::sampleRegular(r,2e5,xy=T))
}
}

## Get shaded relief for background
ucorelief=raster('data/shaded.relief.grey/shadedrelief_uco.tif')
ucorelief.agg=aggregate(ucorelief,2)
ucorelief.df=data.frame(rasterToPoints(ucorelief.agg))
names(ucorelief.df)=c('long','lat','ter')

recon.version='v3.1'
covrange='idp1'
scalesnotel='scale'
product='phvrcn'
cost='r2'
fscaMatch='wofsca'
style='real-time'
residblending='unblended'
if(residblending=='unblended'){
	resid=''
} else {
	resid='full'
}
config=''#'bic-v2.5-removed_NWbdiff'
dateflag='B'
graphbase=paste0('graphs/rswe_',recon.version,'/covrange',covrange,'/snotel',scalesnotel,'/',dateflag,'/',fscaMatch,'/',style,'/')

sample=NULL
rdfswe=data.frame()
for(product in c('phv','phvrcn')){
	for(yr in c(2011,2012)){
		swe=get_regsurf(yr,product,sample)
		dF=data.frame(swe)#need this if using sampleRegular in get_regsurf so keep it
		colnames(dF)=c('x','y','swe')
		if(product=='phv') {
			dF$product='PHV' 
		} else {
			dF$product='PHVRCN'
		}
		dF$yr=yr
		str(dF)
		rdfswe=rbind(rdfswe,dF)
	}
}
rdfswe$swecut=cut(rdfswe$swe,breaks=c(0,0.1,0.25,0.5,1,2,5))


huc4=readOGR('data/gis','UpperCRB')
huc4.gg=fortify(huc4)

## facet labeller
facet_labeller <- function(variable,value){
  return(product_names[value])
}

product_names=list(
	'PHV'='PHV-baseline',
	'PHVRCN'='PHV-RCN')

## extract legend from plot
glegend<-function(a.ggplot){
  tmp <- ggplot_gtable(ggplot_build(a.ggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# rdfswe=data.frame(swe=runif(9*4),product=rep(c('a','b',2)),yr=rep(c(2010,2011),2))
# rdfswe$swecut=cut(rdfswe$swe,breaks=seq(0,1,.2))

numcolors=length(levels(rdfswe$swecut))+1
mycols=brewer.pal(numcolors,'YlGnBu')[2:numcolors]
g1=ggplot(rdfswe)+
	geom_raster(aes(x,y,fill=swecut))+
	scale_fill_manual(name='SWE [m]',
		values=mycols,
		guide=guide_legend(rev=T),
		drop=F)+
	geom_raster(data=ucorelief.df,aes(long,lat,alpha=ter))+
	scale_alpha_continuous(range=c(0.5,0.09),guide=F)+
	geom_path(data=huc4.gg,aes(x=long,y=lat,group=id),color='grey30',size=.15)+
	coord_fixed(ratio=1,xlim=c(-112.25,-104.125),ylim=c(33,43.75))+
	scale_y_continuous(breaks=seq(33,43,2),labels=seq(33,43,2))+
	labs(x='Longitude',y='Latitude')+
	facet_grid(yr~product,labeller=facet_labeller)+
	theme_minimal()+
	theme(strip.text=element_text(size=12,face='bold'),
		panel.grid=element_blank(),
		legend.key=element_blank()
		)
		# legend.direction='horizontal',
		# legend.text.align=-1,
		# legend.position=c(.5,.6),
		# legend.justification=c(.5,.5),
		# legend.text=element_text(angle=45))
ggsave(plot=g1,file=paste0(graphbase,'spatialmap_',residblending,'_2011_2012_Apr1_wBasins.png'),dpi=450,width=7,height=7)

## --- save the legend by itself to rotate in inkscape
g=ggplot(rdfswe)+
	geom_raster(aes(x,y,fill=swecut))+
	scale_fill_manual(name='SWE [m]',
		values=mycols,
		drop=F,
		guide=guide_legend(rev=F))+
	theme(legend.key=element_blank())
pdf(paste0(graphbase,'spatialmap_legend.pdf'))
print(grid.arrange(glegend(g)))
dev.off()
## --






# lc$type=with(lc,ifelse(value==41 | value==42 | value==43, 'forest',
# 	ifelse(value==52 | value==71,'grassland',
# 	ifelse(value==81 | value==82, 'ag',
# 	ifelse(value==31, 'barren', NA)))))
# lc=format(mutate(lc,
# 	pct= count/5031000), scientific=F)

##
rdfswediff=ddply(rdfswe[,c('x','y','product','swe','yr')],.(yr),function(dF){
	dF2=dcast(dF,x+y~product,value.var='swe')
	dF2$Difference=with(dF2,PHVRCN-PHV)
	dF2$Difference[dF2$PHVRCN==0]=NA
	return(dF2)
})
rdfswediff.m=melt(rdfswediff[,c('x','y','yr','Difference')],.(yr,x,y))
rdfswediff.m$diffcut=cut(rdfswediff.m$value,breaks=c(-5,-1,-0.25,-0.05,0.05,0.25,1,5))

ggplot(rdfswediff.m)+
	geom_raster(aes(x,y,fill=diffcut))+
	scale_fill_brewer(name='SWE Difference [m]',
			palette='PuOr',
			drop=FALSE,
			guide=guide_legend(rev=T))+
	geom_raster(data=ucorelief.df,aes(long,lat,alpha=ter))+
	scale_alpha_continuous(range=c(0.45,0.03),guide=F)+
	geom_path(data=huc4.gg,aes(x=long,y=lat,group=id),color='grey30',size=0.15)+
	coord_fixed(ratio=1,xlim=c(-112.25,-104.125),ylim=c(33,43.75))+
	scale_y_continuous(breaks=seq(33,43,2),labels=seq(33,43,2))+
	labs(x='Longitude',y='Latitude')+
	facet_grid(yr~variable)+
	theme_minimal()+
	theme(strip.text=element_text(size=12,face='bold'),
		panel.grid=element_blank(),
		legend.key=element_blank()
		)

ggsave(paste0(graphbase,'spatialmap_',residblending,'_mdldiff_2011_2012_Apr1_wBasins.png'),dpi=450,width=7,height=7)


##--- save the legend as pdf to rotate in inkscape
gdiff=ggplot(rdfswediff.m)+
	geom_raster(aes(x,y,fill=diffcut))+
	scale_fill_brewer(name='SWE Difference [m]',
			palette='PuOr',
			drop=FALSE,
			guide=guide_legend(rev=F))+
		theme(legend.key=element_blank())	
# show(gdiff)
pdf(paste0(graphbase,'spatialmapdiff_legend.pdf'))
print(grid.arrange(glegend(gdiff)))
dev.off()
# --


# ##---
# rdf=mutate(rdf,
# 	apr2012cut=cut(apr2012,c(0,0.5,1,1.5,2,2.5,3,5,10)),
# 	apr2011cut=cut(apr2011,c(0,0.5,1,1.5,2,2.5,3,5,10)),
# 	apr2010cut=cut(apr2010,c(0,0.5,1,1.5,2,2.5,3,5,10)),
# 	diff1211=apr2012-apr2011,
# 	diff1211cut=cut(diff1211,seq(-1.5,1.5,.25)),
# 	diff1210=apr2012-apr2010,
# 	diff1210cut=cut(diff1210,seq(-1.5,1.5,.25)),
# 	diff1110=apr2011-apr2010,
# 	diff1110cut=cut(diff1110,seq(-1.5,1.5,.25)))

# range(rdf$apr2012,na.rm=T)


# sample.df <- function(df, n) df[sample(nrow(df), n), , drop = FALSE]
# rdfswe=rdf[,c('x','y','apr2011cut','apr2012cut')]
# test=sample.df(rdfswe,4e6)
# rdfswe2=melt(rdfswe,.(x,y))

# ggplot(rdfswe2)+
# 	geom_raster(aes(x,y,fill=value))+
# 	facet_wrap(~variable,scales='free')+
# 	labs(title='april 1')+
# 	scale_fill_brewer(palette='YlGnBu',guide=guide_legend(rev=T),drop=F)+
# 	theme(panel.background=element_rect('black'),
# 		panel.grid=element_blank())
# 	# expand_limits(x=c(0,0),y=c(0,0))

# rdfdiff=rdf[,c('x','y',apr2011)]


# rdf2=melt(rdf,.(x,y))
# rdf2sub=subset(rdf2,grepl('apr2011',as.character(variable)))
# rdf2sub$variable=droplevels(rdf2sub$variable)
# rdf3sub=subset(rdf2sub,!grepl('cut',as.character(variable)))



# ggplot(test)+
# 	geom_raster(aes(x,y,fill=value))+
# 	facet_wrap(~variable)+
# 	labs(title='april 1, 2012')+
# 	scale_fill_brewer(palette='YlGnBu',guide=guide_legend(rev=T))+
# 	theme(panel.background=element_rect('black'))+
# 	expand_limits(x=c(0,0),y=c(0,0))


# ggplot(rdf)+
# 	geom_raster(aes(x,y,fill=diff1211cut,hpad=0))+
# 	labs(title='april 1 difference, 2012-2011')+
# 	scale_fill_brewer(palette='PuOr')+
# 	theme(panel.background=element_rect('black'))

# dev.new()
# ggplot(rdf)+
# 	geom_raster(aes(x,y,fill=diff1210cut))+
# 	labs(title='april 1 difference, 2012-2010')+
# 	scale_fill_brewer(palette='PuOr')+
# 	theme(panel.background=element_rect('black'))

# dev.new()
# ggplot(rdf)+
# 	geom_raster(aes(x,y,fill=diff1110cut),hpad=0,vpad=0)+
# 	labs(title='april 1 difference, 2011-2010')+
# 	scale_fill_brewer(palette='PuOr')+
# 	theme(panel.background=element_rect('black'))



# range(swestats_all$avg)