doCOST=function(mdlpreds,cost,recon.version){
     #pass in predictions and a cost function as character  (cor, rmse etc)     
     
     if(cost=='cor') { #statistics to maximize
          costfun<-function(y,yhat) cor(y,yhat)
     } else if (cost=='r2'){     
          costfun<-function(y,yhat) 1-sum((y-yhat)^2)/sum((y-mean(y,na.rm=T))^2)
     }
     if(exists(costfun)){
          daystats=ddply(mdlpreds,.(yrdoy,recondate), 
                         function(x){
                              summarise(x,  phvcor=costfun(snotel,phv.fullpred),
                                   phvrcncor=costfun(snotel,phvrcn.fullpred),
                                   reconcor=costfun(snotel,recon))
                                   }
                         )
          daystats$date=as.POSIXct(strptime(daystats$yrdoy,'%Y%j'))
          whichrdate=ddply(daystats,.(date),summarise,
                           yr=strftime(unique(date),'%Y'),
                           phvrcn_recondate=recondate[which.max(phvrcncor)],
                           reconcordate=recondate[which.max(reconcor)])
     } else { #statistics to minimize
          if(cost=='mae') costfun<-function(y, yhat) mean(abs(yhat-y),na.rm=T)
          if(cost=='rmse') costfun<-function(y,yhat) sqrt(mean((yhat-y)^2),na.rm=T)
     
          daystats=ddply(mdlpreds,.(yrdoy,recondate),function(x){
                         summarise(x, phvcor=costfun(snotel,phv.fullpred),
                            phvrcncor=costfun(snotel,phvrcn.fullpred),
                            reconcor=costfun(snotel,recon))
                              }
                         )
          daystats$date=as.POSIXct(strptime(daystats$yrdoy,'%Y%j'))
          whichrdate=ddply(daystats,.(date),summarise,
                           yr=strftime(unique(date),'%Y'),
                           phvrcn_recondate=recondate[which.min(phvrcncor)],
                           reconcordate=recondate[which.min(reconcor)])
     }     
     
     return(whichrdate)
}

