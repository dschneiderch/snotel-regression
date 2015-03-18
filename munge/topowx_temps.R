library('ProjectTemplate')
setwd('~/Documents/snotel-regression_project')
load.project(list(cache_loading=T))


fn.tmin='data/topowx-stationdata/raw/topowx_final_stn_obs_tmin.nc'
fn.tmax='data/topowx-stationdata/raw/topowx_final_stn_obs_tmax.nc'
# fn=fn.tmin
# tvar='tmin'
twx.t=function(fn,tvar){
     #open netcdf file
     nc=nc_open(fn)
     #extract time dimension values
     zeit=nc$dim[[1]]$vals#these are integers since tstart
     #get indices for dates of interest
     tstart=as.POSIXct('1948-01-01',tz='MST')
     twx.dates=tstart+zeit*24*3600
     tneed1=as.POSIXct('1999-10-01',tz='MST')
     tneed2=as.POSIXct('2012-09-30',tz='MST')
     t1ind=which(twx.dates==tneed1)
     t2ind=which(twx.dates==tneed2)
     zeitind=seq(t1ind,t2ind,1)
     #
     #extract station id dimension values
     stndim=nc$dim[[2]]$vals
     #use station_id from snotelrec get find indices and make station lsit
     stationid=unique(snotelrec$Station_ID)
     twx.stnid=unlist(strsplit(stndim,'_'))[seq(2,nrow(stndim)*2,2)]
     snotelind=which(twx.stnid %in% stationid)
     staid=twx.stnid[snotelind]
#      twx.stnid=laply(stndim,function(x) laply(strsplit(x,'_'),'[',2)) #slow
#      snotelind=twx.stnid %in% unique(snotelrec$Station_ID)
#      staid=twx.stnid[snotelind]
     #
     #extract temp variable from netcdf file. only grab indices calculated above
     twx.temp=ncvar_get(nc,tvar,start=c(snotelind[1],t1ind),count=c(length(snotelind),length(zeitind)))
     #rotate. netcdf convention is time as last dimension so make stations cols and time rows     
     snoteltemp=t(twx.temp)
     #
     #create dataframe of temp, station ids and date
     tt=data.frame(snoteltemp)
     colnames(tt)=staid
     tt$date=seq(tneed1,tneed2,24*3600)
     twx.snotel=melt(tt,.(date))
     colnames(twx.snotel)=c('date','Station_ID',tvar)
     return(twx.snotel)
}
     

#extract tmin
twx.tmin=twx.t(fn.tmin,'tmin')
#extract tmax
twx.tmax=twx.t(fn.tmax,'tmax')
#combine variables into 1 dataframe
twx.snotel=cbind(twx.tmin,tmax=twx.tmax$tmax)

#2 other stations not in original nc files
sta=c("05G04S","08S08S")
stadf=ldply(sta,function(x){
     fn=paste0('data/topowx-stationdata/raw/SNOTEL_',x,'.csv')
     df=read.csv(fn,header=T)
     df$DATE=as.POSIXct(df$DATE,tz='MST')
     df=df[df$DATE>=as.POSIXct('1999-10-01',tz='MST') & df$DATE<=as.POSIXct('2012-09-30',tz='MST'),]
     dfout=data.frame(date=df$DATE,Station_ID=x)
     dfout$tmin=with(df,ifelse(TMIN_HOMOG==-9999,TMIN_INFILL,TMIN_HOMOG))
     dfout$tmax=with(df,ifelse(TMAX_HOMOG==-9999,TMAX_INFILL,TMAX_HOMOG))    
     return(dfout)
})
twx.snotel=rbind(twx.snotel,stadf)

#calculate tavg
twx.snotel=mutate(twx.snotel,
               tavg=colMeans(t(matrix(c(tmin,tmax),ncol=2))))
#calculate degdays from tavg
twx.snotel$degday=twx.snotel[,'tavg']
#twx.snotel$degday=ifelse(twx.snotel$degday<0,0,twx.snotel$degday)
#calc water year
twx.snotel$wyr=water_yr(twx.snotel$date)
#calc cumulative degday by water year and station
twx.snotel=ddply(twx.snotel,.(Station_ID,wyr),function(x){
     mutate(x,
            cumdegday=cumsum(degday))
})
twx.snotel$cumdegday=ifelse(twx.snotel$cumdegday<0,0,twx.snotel$cumdegday)

#write for later use.
write.table(twx.snotel,file='data/topowx-stationdata/twx.cumdegday.txt',sep='\t',row.names=F,quote=F)
