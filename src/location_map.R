setwd('~/Documents/snotel-regression_project')
library(ProjectTemplate)
# library(raster)
# library(rgdal)
# library(ggplot2)
# library(reshape2)
# library(plyr)
load.project(list(cache_loading=T))

nat.earth <- stack('data/gis/NE1_LR_LC_SR_W_DR/NE1_LR_LC_SR_W_DR.tif')
names(nat.earth)=c('red','green','blue')
ne_lakes <- readOGR('data/gis/ne_50m_lakes',   'ne_50m_lakes')
ne_rivers <- readOGR('data/gis/ne_50m_rivers_lake_centerlines','ne_50m_rivers_lake_centerlines')
huc4=readOGR('data/gis/','UpperCRB')
cities=readOGR('data/gis/ne_10m_populated_places_simple','ne_10m_populated_places_simple')

#  I have a domain I'm interested in, but there's no reason you can't define something else:
quick.subset <- function(x, longlat){
# longlat should be a vector of four values: c(xmin, xmax, ymin, ymax)
  x@data$id <- rownames(x@data)
  x.f = fortify(x, region="id")
  x.join = join(x.f, x@data, by="id")
  x.subset <- subset(x.join, x.join$long > longlat[1] & x.join$long < longlat[2] &
                           x.join$lat > longlat[3] & x.join$lat < longlat[4])
  x.subset
}

domain <- c(-112.5,-104,33,44)
lakes.subset <- quick.subset(ne_lakes, domain)
river.subset <- quick.subset(ne_rivers, domain)
nat.crop <- crop(nat.earth, y=extent(domain))
rast.table <- data.frame(xyFromCell(nat.crop, 1:ncell(nat.crop)),
                         getValues(nat.crop/255))
rast.table$rgb <- with(rast.table, rgb(red,green,blue,1))

sntl.df=as.data.frame(snotellocs)
huc4.df=fortify(huc4)
states=map_data('state')
coi=c('Denver','Salt Lake City','Phoenix','Laramie','Grand Junction','Santa Fe','Lander')
cities=subset(cities,name %in% coi)
cities.df=as.data.frame(cities)

# et voila!
ggplot(data = rast.table, aes(x = x, y = y)) +
  geom_raster(fill = rast.table$rgb) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed(ratio=1,xlim=c(-112.45,-104.122),ylim=c(33,43.85))+
  labs(x='',y='')+#labs(x='Longitude',y='Latitude')+
  theme_bw()+
  theme(axis.text=element_blank(),
    axis.ticks=element_blank())
  # theme(axis.text=element_text(color='black',size=14,face='bold'))
 
 ggsave(file='graphs/map_background_uco.png',dpi=400,width=95,units='mm')

ggplot() +
  geom_path(data=states,aes(x=long,y=lat,group=group),color='grey35')+
  geom_path(data=huc4.df,aes(long,lat,group=group),color='grey20')+
  geom_point(data=sntl.df,aes(x=Longitude,y=Latitude),size=1.5,color='blue')+
  geom_point(data=cities.df,aes(coords.x1,coords.x2),shape=8,size=3)+
  geom_text(data=cities.df,aes(x=coords.x1,y=coords.x2+.2,label=name),size=4,vjust=0,hjust=0.23,fontface='bold')+
  scaleBar(lon=-107,lat=43.5,distanceLon=100,distanceLat=20,distanceLegend=-30,dist.unit="km",orientation=TRUE,arrow.distance=-1000,arrow.length=80)+# geom_polygon(data=lakes.subset, aes(x = long, y = lat, group = group), fill = '#ADD8E6') +
  # geom_path(data=river.subset, aes(x = long, y = lat, group = group), color = 'blue') +
  # geom_path(data=lakes.subset, aes(x = long, y = lat, group = group), color = 'blue') +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed(ratio=1,xlim=c(-112.45,-104.122),ylim=c(33,43.85))+
  labs(x='Longitude',y='Latitude')+
  theme_bw()+
  theme(axis.text=element_text(color='black',size=14,face='bold'))
 
ggsave(file='graphs/map_vector_uco.pdf',width=95,units='mm')
 
### usa inset map
ucodomain=as(extent(nat.crop),'SpatialPolygons')
ucodomain.df=fortify(ucodomain)

statenames=data.frame(name=c('UT','CO','AZ','NM','WY','ID'),long=c(-111,-106,-111,-106,-108,-114),lat=c(39,39,35,35,43,43.5))
ggplot()+
  	geom_path(data=states,aes(long,lat,group=group))+
  	geom_path(data=ucodomain.df,aes(long,lat, group=group),color='red')+
  	geom_text(data=statenames,aes(x=long,y=lat,label=name),size=2)+
  	labs(x='',y='')+
  	theme_minimal()+
  	theme(axis.ticks=element_blank(),
  		axis.text=element_blank(),
  		panel.grid=element_blank())

  ggsave(file='graphs/map_usa.png',height=2.5,width=4) 

