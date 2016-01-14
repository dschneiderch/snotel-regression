## TODO plot both blended and unblended estimates, create 1 dataframe to facet.
calc_GMoran=function(dF,snotellocs.usgs,covrange){

xrdate=unique(dF$recondate)
dt=unique(dF$date)
print(dt)

if( length(unique(dF$snotel))>1 ){
     snotellocs.df=as.data.frame(snotellocs.usgs)
     if(grepl('km',covrange)) {
          maxdist=as.numeric(gsub('km','',covrange))*1000
     } else maxdist=300000
     
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
     kd.idw=nb2listw(kd,glist=idw,style='B')
#      kd.binary_weight=nb2listw(kd200,style='W')#this just weights the nearest neighbors for each station equally. W is row standardized.
     coordnames(snotellocs.usgs)=c('x','y')
     locsmat=as.matrix(as.data.frame(snotellocs.usgs)[,c('x','y')])    
     exdat=data.frame(locsmat,dF)
     localmoran4plot=function(localm){
          dF=cbind(as.data.frame(localm),x=snotellocs.df$x,y=snotellocs.df$y)
          dF=mutate(dF,
                    Pcut=cut(dF$P,c(0,0.01,0.05,0.1,max(dF$P))),
                    shapecut=as.factor(sign(Ii)))
          #dF$Pcut=as.factor(dF$Pcut,levels=rev(levels(dF$Pcut)))
     }
# localm.phv=localmoran4plot(localmoran(exdat$phvresid,kd.idw))
# localm.phvrcn=localmoran4plot(localmoran(exdat$phvrcnresid,kd.idw))

mt1=moran.test(exdat$phvresid,listw=kd.idw)
mt2=moran.test(exdat$phvrcnresid,listw=kd.idw)
mt3=moran.test(exdat$phv.fullpred-exdat$snotel,listw=kd.idw)
mt4=moran.test(exdat$phvrcn.fullpred-exdat$snotel,listw=kd.idw)

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

} else {
     return(
          data.frame(
          date=dt,
          recondate=xrdate,
          phv.globalI=NA,
          phvrcn.globalI=NA,
          phv.pvalue=NA,
          phvrcn.pvalue=NA,
          phv.full.globalI=NA,
          phvrcn.full.globalI=NA,
          phv.full.pvalue=NA,
          phvrcn.full.pvalue=NA)
     )
}

}