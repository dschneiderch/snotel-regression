doPHVfit=function(dF){
     
     mdl=tryCatch({
          ### Base model
          ## start.param=fitdistr(swe2model$snotel,'gamma')
          ## lmmodel=lm(log(snotel) ~ Lat+Elev+Eastness+Northness+Slope+RegionalSlope+RegionalEastness+RegionalNorthness+FtprtW+Wbdiff+NWbdiff+Wd2ocean+1,data=dF)
          #startvals=c(log(mean(dF$snotel)),rep(0,11))
          #print(startvals)
          dF$snotel=dF$snotel+runif(1,0,0.00001)#doesn't like converging with 0s in snotel
#           print(dF[,c('date','snotel')])#
          #
          glm(snotel ~ Lat+Elev+Eastness+Northness+Slope+RegionalSlope+RegionalEastness+RegionalNorthness+FtprtW+Wbdiff+NWbdiff+Wd2ocean+1,family=gaussian(link='log'),data=dF)},#,start=startvals)},
          error= function(e) {
               print(paste0('returning empty phv model ', dF$date[1]))
          }
     )

     if(!is.character(mdl)) {
          #
          ### Reduce with stepaic
          ## PHV only
          mdl.stepaic.phv=stepAIC(mdl, scope=list(upper= ~ Lat + Elev + Eastness + Northness + Slope + RegionalSlope + RegionalEastness + RegionalNorthness + FtprtW + Wbdiff + NWbdiff + Wd2ocean + 1, lower=~1), direction='both',trace=F,k=log(nrow(dF)))
          #   
          
     }    else     { mdl.stepaic.phv=NA }
     
     return(mdl.stepaic.phv)
}