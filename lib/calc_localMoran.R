## TODO plot both blended and unblended estimates, create 1 dataframe to facet.
calc_localMoran=function(dF,snotellocs.usgs,recon.version,covrange){

xrdate=unique(dF$recondate)
dt=unique(dF$date)
print(dt)

if( length(unique(dF$snotel))>1 ){
     snotellocs.df=as.data.frame(snotellocs.usgs)
     maxdist=as.numeric(gsub('km','',covrange))*1000
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
     p=.5#fassnacht found 0.5 to be optimal until mid-march, then up to 1.3 at end-april then back down to .5.
     idw=lapply(dist,function(x) 1/(x^p))
     kd.idw=nb2listw(kd,glist=idw,style='B')#use 'B' if yu want the values from idw as the weights.
#      kd.binary_weight=nb2listw(kd200,style='W')#this just weights the nearest neighbors for each station equally. W is row standardized.
     locsmat=as.matrix(as.data.frame(snotellocs.usgs)[,c('x','y')])    
     exdat=data.frame(locsmat,dF)
     localmoran4plot=function(localm){
          dF=cbind(as.data.frame(localm),x=snotellocs.df$x,y=snotellocs.df$y)
          dF=mutate(dF,
                    Pcut=cut(dF$Pr,c(0,0.01,0.05,0.1,max(dF$Pr))),
                    shapecut=as.factor(sign(Ii)))
          #dF$Pcut=as.factor(dF$Pcut,levels=rev(levels(dF$Pcut)))
     }
localm.snotel=localmoran4plot(localmoran(exdat$snotel,kd.idw))
localm.phv=localmoran4plot(localmoran(exdat$phv.pred,kd.idw))
localm.phvrcn=localmoran4plot(localmoran(exdat$phvrcn.pred,kd.idw))
localm.phv.full=localmoran4plot(localmoran(exdat$phv.fullpred,kd.idw))
localm.phvrcn.full=localmoran4plot(localmoran(exdat$phvrcn.fullpred,kd.idw))


require(maps)
states_map <- map_data("state")
coordinates(states_map)=~long+lat
proj4string(states_map)=CRS('+proj=longlat +datum=WGS84')
states.usgs=as.data.frame(spTransform(states_map,CRS('+init=epsg:5070')))
usgsbox.df=as.data.frame(usgsbox)


     # p=ggplot(localm.phv)+
     #      geom_point(aes(x=x,y=y,size=abs(Ii),shape=shapecut,colour=Pcut))+
     #      geom_path(data=states.usgs,aes(x,y,group=group))+
     #      annotate('text',x=-9e5,y=1.25e6,
     #               label=paste0('Global I = ', format(mt1$estimate[1],digits=3),'\n','p-value = ',format(mt1$p.value,digits=3)))+
     #      labs(title='phv regression (no blend)')+
     #      coord_cartesian(xlim=extent(usgsbox)[1:2],ylim=extent(usgsbox)[3:4])+
     #      theme_minimal()


     # pr=ggplot(localm.phvrcn)+
     #      geom_point(aes(x=x,y=y,size=abs(Ii),shape=shapecut,colour=Pcut))+
     # geom_path(data=states.usgs,aes(x,y,group=group))+
     #      annotate('text',x=-9e5,y=1.25e6,
     #          label=paste0('Global I = ', format(mt2$estimate[1],digits=3),'\n','p-value = ',format(mt2$p.value,digits=3)))+
     #      labs(title='phvrcn regression (no blend)')+
     # coord_cartesian(xlim=extent(usgsbox)[1:2],ylim=extent(usgsbox)[3:4])+
     #      theme_minimal()

#ga=arrangeGrob(p,pr)

#ggsave(ga,width=6,height=10,filename=paste0('diagnostics/rswe_',recon.version,'/localmoran/covrange',covrange,'/xval_localmoran_',dt,'_recondate-',xrdate,'_recon_',recon.version,'.pdf'))
 #    moran.plot(exdat$phvrcn.pred,kd.idw)


return(
     data.frame(
     date=dt,
     recondate=xrdate,
     Station_ID=dF$Station_ID,
     snotel.I=localm.snotel$Ii,
     snotel.pvalue=localm.snotel$Pr,
     phv.I=localm.phv$Ii,
     phvrcn.I=localm.phvrcn$Ii,
     phv.pvalue=localm.phv$Pr,
     phvrcn.pvalue=localm.phvrcn$Pr,
     phv.full=localm.phv.full$Ii,
     phvrcn.full=localm.phvrcn.full$Ii,
     phv.full.pvalue=localm.phv.full$Pr,
     phvrcn.full.pvalue=localm.phvrcn.full$Pr)
     )

} else {
     return(
          data.frame(
          date=dt,
          recondate=xrdate,
          Station_ID=dF$Station_ID,
          snotel.I=NA,
          snotel.pvalue=NA,
          phv.I=NA,
          phvrcn.I=NA,
          phv.pvalue=NA,
          phvrcn.pvalue=NA,
          phv.full=NA,
          phvrcn.full=NA,
          phv.full.pvalue=NA,
          phvrcn.full.pvalue=NA)
          )
}

}