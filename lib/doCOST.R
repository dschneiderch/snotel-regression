doCOST=function(mdlpreds,cost,recon.version){
     #pass in predictions and a cost function as character  (cor, rmse etc)     
     
     if(cost=='cor') { #statistics to maximize
          costfun<-function(y,yhat) cor(y,yhat)
          opti=max
          which.opti=which.max
     } else if (cost=='r2') {     
          costfun<-function(y,yhat) 1-sum((y-yhat)^2)/sum((y-mean(y,na.rm=T))^2)
          opti=max
          which.opti=which.max
     } else if(cost=='mae') { #statistics to minimize
          costfun<-function(y, yhat) mean(abs(yhat-y),na.rm=T)
          opti=min 
          which.opti=which.min
     } else if(cost=='rmse') {
          costfun<-function(y,yhat) sqrt(mean((yhat-y)^2),na.rm=T)
          opti=min
          which.opti=which.min
     }
     
     daystats=ddply(mdlpreds,.(yrdoy,recondate),function(x){
          summarise(x, 
                    phv.cor=costfun(snotel,phv.pred),
                    phvrcn.cor=costfun(snotel,phvrcn.pred),
                    phvfull.cor=costfun(snotel,phv.fullpred),
                    phvrcnfull.cor=costfun(snotel,phvrcn.fullpred),
                    recon.cor=costfun(snotel,recon),
                    reconfull.cor=costfun(snotel,recon.fullpred))
     }
     )
     daystats$date=as.POSIXct(strptime(daystats$yrdoy,'%Y%j'))
     
     whichrdate=ddply(daystats,.(date),summarise,
                      yr=strftime(unique(date),'%Y'),
                      phvrcn_recondate=recondate[which.opti(phvrcncor)],
                      recon_costdate=recondate[which.opti(reconcor)],
                      skill_phv=opti(phv.cor),
                      skill_phvrcn=opti(phvrcn.cor),
                      skill_recon=opti(recon.cor),
                      skill_phvfull=opti(phvfull.cor),
                      skill_phvrcnful=opti(phvrcnfull.cor),
                      skill_reconfull=opti(reconfull.cor))
     
     return(whichrdate)
}