opt_LUT=function(){
# create opt lookup table
topowx=read.table('data/topowx-stationdata/twx.cumdegday.txt',header=T,sep='\t')
topowx$date=as.POSIXct(topowx$date,tz='MST')
topowx=arrange(topowx,date,Station_ID)
#
snotelrec=arrange(snotelrec,date,Station_ID)
#
optLUT=data.frame(topowx[,c('date','Station_ID','cumdegday')],apcp=snotelrec$apcp)
optLUT$wy=water_yr(optLUT$date)
optLUT$doy=strftime(optLUT$date,'%j')
optlist=dlply(optLUT,.(Station_ID),function(dF){
  # print(dF)
  dF$cumdegday=scale(dF$cumdegday)
  dF$apcp=scale(dF$apcp)
  # dF$euc=with(dF,sqrt(cumdegday^2+apcp^2))
  # euc=with(dF,ifelse(atan(apcp/cumdegday)*180/pi<45 & atan(apcp/cumdegday)*180/pi >= -135,-euc,euc))

  stascale=data.frame(cddavg=attr(dF$cumdegday,'scaled:center'),cddstd=attr(dF$cumdegday,'scaled:scale'),
    apcpavg=attr(dF$apcp,'scaled:center'),apcpstd=attr(dF$apcp,'scaled:scale'))
  return(list(dF,stascale))
})
optLUT=ldply(optlist,'[[',1)
optLUT$mth=as.numeric(strftime(optLUT$date,'%m'))
optLUT$yr=as.numeric(strftime(optLUT$date,'%Y'))
# optLUT=optLUT[optLUT$mth==3 | optLUT$mth==4 | optLUT$mth==5 | optLUT$mth==6,]#this needs to correspond to recon dates
#
stascale=ldply(optlist,'[[',2)
return(list(optLUT,stascale))
}