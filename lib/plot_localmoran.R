## TODO plot both blended and unblended estimates, create 1 dataframe to facet.

plot_localmoran=function(dF,snotellocs.usgs,recon.version){
     snotellocs.df=as.data.frame(snotellocs.usgs)
     maxdist=200000
     kd=dnearneigh(snotellocs.usgs,d1=0,d2=maxdist,row.names=snotellocs.usgs$Station_ID)#d2 should equal abut 200km fro all stations to have some neighbors. min 8 with 200km.
     #summary(kd)
    # plot(kd,coords=coordinates(snotellocs.usgs))#
     getneighbors=function(nm){
          ind=grep(nm,attr(kd,'region.id'))
          vec=kd[[ind]]
          data.frame(station=nm,neighbors=attr(kd,'region.id')[vec])
     } 
     snotelneighbors=ldply(as.list(attr(kd,'region.id')),getneighbors)
#
     dist=nbdists(kd,snotellocs.usgs)
     p=1
     idw=lapply(dist,function(x) 1/(x^p))
     kd.idw=nb2listw(kd,glist=idw,style='W')
#      kd.binary_weight=nb2listw(kd200,style='W')#this just weights the nearest neighbors for each station equally. W is row standardized.
     locsmat=as.matrix(as.data.frame(snotellocs.usgs)[,c('x','y')])    
     exdat=data.frame(locsmat,dF)
     localmoran4plot=function(localm){
          dF=cbind(as.data.frame(localm),x=snotellocs.df$x,y=snotellocs.df$y)
          dF=mutate(dF,
                    Pcut=cut(dF$P,c(0,0.01,0.05,0.1,max(dF$P))),
                    shapecut=as.factor(sign(Ii)))
          #dF$Pcut=as.factor(dF$Pcut,levels=rev(levels(dF$Pcut)))
     }
localm.phvrcn=localmoran4plot(localmoran(exdat$phvrcn.pred,kd.idw))
localm.phv=localmoran4plot(localmoran(exdat$phv.pred,kd.idw))

require(maps)
states_map <- map_data("state")
coordinates(states_map)=~long+lat
proj4string(states_map)=CRS('+proj=longlat +datum=WGS84')
states.usgs=as.data.frame(spTransform(states_map,CRS('+init=epsg:5070')))
usgsbox.df=as.data.frame(usgsbox)

mt1=moran.test(exdat$phv.pred,listw=kd.idw)
     p=ggplot(localm.phv)+
          geom_point(aes(x=x,y=y,size=abs(Ii),shape=shapecut,colour=Pcut))+
          geom_path(data=states.usgs,aes(x,y,group=group))+
          annotate('text',x=-9e5,y=1.25e6,
                   label=paste0('Global I = ', format(mt1$estimate[1],digits=3),'\n','p-value = ',format(mt1$p.value,digits=3)))+
          labs(title='phv regression (no blend)')+
          coord_cartesian(xlim=extent(usgsbox)[1:2],ylim=extent(usgsbox)[3:4])+
          theme_minimal()

mt2=moran.test(exdat$phvrcn.pred,listw=kd.idw)
     pr=ggplot(localm.phvrcn)+
          geom_point(aes(x=x,y=y,size=abs(Ii),shape=shapecut,colour=Pcut))+
     geom_path(data=states.usgs,aes(x,y,group=group))+
          annotate('text',x=-9e5,y=1.25e6,
              label=paste0('Global I = ', format(mt2$estimate[1],digits=3),'\n','p-value = ',format(mt2$p.value,digits=3)))+
          labs(title='phvrcn regression (no blend)')+
     coord_cartesian(xlim=extent(usgsbox)[1:2],ylim=extent(usgsbox)[3:4])+
          theme_minimal()

ga=arrangeGrob(p,pr)

xrdate=unique(dF$recondate)
dt=unique(dF$date)
ggsave(ga,width=6,height=10,filename=paste0('diagnostics/rswe_',recon.version,'/localmoran/xval_localmoran_',dt,'_recondate-',xrdate,'_recon_',recon.version,'.pdf'))
 #    moran.plot(exdat$phvrcn.pred,kd.idw)

mt3=moran.test(exdat$phv.fullpred,listw=kd.idw)
mt4=moran.test(exdat$phvrcn.fullpred,listw=kd.idw)

return(
     data.frame(
     date=dt,
     recondate=xrdate,
     phv.globalI=mt1$estimate[1],
     phvrcn.globalI=mt2$estimate[1],
     phv.pvalue=mt1$p.value,
     phvrcn.pvalue=mt2$p.value,
     phv.full.globalI=mt3$estimate[1],
     phvrcn.full.globalI=mt4$estimate[1],
     phv.full.pvalue=mt3$p.value,
     phvrcn.full.pvalue=mt4$p.value)
     )

}