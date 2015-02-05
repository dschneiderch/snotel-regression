check_spatialcovariance <- function(dF,locs,colnm,recon.version=NULL){
     if(colnm!='snotel') {
          if(nargs()!=4) {
               print('need to input recon version')
               break
          }
          plottitle=paste0(colnm,recon.version,'_',unique(dF$date))     
     } else { 
          plottitle=paste0(colnm,'_',unique(dF$date))
     }
     
#     colind=grep(glob2rx(colnm),colnames(dF))
#     if(is.na(dF$get(colnm))) print(length(is.na(dF$get(colnm))))
#      #
     locsmat=as.matrix(as.data.frame(locs)[,c('x','y')])
     exdat=data.frame(locsmat,dF)
     coordinates(exdat) <- ~x+y
     vgm=variogram(get(colnm)~1,data=exdat)
     
     pdf(file=paste0('diagnostics/variograms/',plottitle,'_variogram.pdf'))
     print(plot(vgm,main=plottitle))
     dev.off()
     
#      kd200=dnearneigh(locs,d1=0,d2=200000,row.names=locs$Station_ID)
#      kd200.binary_weight=nb2listw(kd200,style='B')
#      dist=nbdists(kd200,locs)
#      idw=lapply(dist,function(x) 1/(x/1000))
#      kd200.binary_weight=nb2listw(kd200,style='B')
#      kd200.idw_weight=nb2listw(kd200,glist=idw,style='B')
#      moran.test(exdat,listw=kd200.binary_weight)
#      
#      localmoran(exdat$phvrcnresid,kd200.idw_weight)
#      moran.plot(locs$phvrcn.resid,kd200.idw_weight)
#      
     
     
}