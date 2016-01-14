#takes static model and fits dyn model with optimized recon
doOptRecon <- function(doydF,static_mdl,style){

if(style=='real-time'){
     foptLUT=optLUT[optLUT$yr<doydF$yr,]
} else { 
     foptLUT=optLUT 
}

     cumdata=optdata[optdata$date %in% unique(doydF$date),]

     mdoy=unique(cumdata$doy)

     LUT=arrange(foptLUT[foptLUT$doy %in% mdoy,],date,Station_ID)
     if((nrow(LUT)%%nrow(cumdata))!=0) stop()

     dmat=data.frame(
          date=unique(doydF$date),
          ryr=LUT$yr,
          Station_ID=LUT$Station_ID,
          historical_apcp=LUT$apcp,
          historical_cumdegday=LUT$cumdegday,
          apcp=cumdata$apcp,
          cumdegday=cumdata$cumdegday,
          dist=sqrt((cumdata$apcp-LUT$apcp)^2+(cumdata$cumdegday-LUT$cumdegday)^2))


# dcast(dmat[,c('yr','Station_ID','d')],yr~Station_ID,value.var='d')
wmat=ddply(dmat,.(Station_ID),function(x){
     mutate(x,
          weight=dist/sum(dist),
          mth=strftime(date,'%m'),
          dy=strftime(date,'%d'))
     # return(x)
     })

# dcast(wmat[,c('yr','Station_ID','w')],yr~Station_ID,value.var='w')
# apply(dcast(wmat[,c('yr','Station_ID','w')],yr~Station_ID,value.var='w'),2,sum) #check to make sure years sum to 1

wmat$recondate=NA
ind=which(as.numeric(wmat$mth) >= 3)
wmat$recondate[ind]=paste(wmat$ryr[ind],wmat$mth[ind],wmat$dy[ind],sep='-')
ind=which(as.numeric(wmat$mth) < 3)
wmat$recondate[ind]=paste0(wmat$ryr[ind],'-03-01')
wmat$recondate=as.POSIXct(wmat$recondate,tz='MST')

snotelrecon=arrange(snotelrecon,date,Station_ID)
snotelrecon.doy=snotelrecon[snotelrecon$date %in% wmat$recondate,]

rmap=data.frame(recon=wmat$weight*snotelrecon.doy$recon,Station_ID=wmat$Station_ID,date=wmat$date)
recon=daply(rmap,.(Station_ID),function(dF) sum(dF$recon))

doyfits=CVwrapper('opt',doydF,static_mdl,recon)

wmat$recon=snotelrecon.doy$recon
wmat$reconopt=rep(recon,each=nrow(wmat)/length(unique(wmat$Station_ID)))
return(list(doyfits,wmat))
}